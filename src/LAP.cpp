#include "LAP.h"

float LAP::sfVersion = 1.3;
const double LAP::sdSlow = 0.01;
const double LAP::sdShigh = 1.;
bool LAP::sbInitialized = false;
unique_ptr<TTree> LAP::tree = nullptr;
int LAP::siHowMany = 0;
const int LAP::ciNpts = 50;

double LAP::sdLaplaceLow[8]		= {0,0,0,0,0,0,0,0};
double LAP::sdLaplaceHigh[8]	= {0,0,0,0,0,0,0,0};

LAP::LAP() {
	if (g_verbose) cout << "LAP c'tor\n";
	LAP::siHowMany++;
}

LAP::LAP(int ch, int length, shared_ptr<Digitizer> digitizer) : Method(ch, length, digitizer) {
	if (g_verbose) cout << "LAP " << id << " c'tor\n";
	LAP::siHowMany++;
	if ((id >= MAX_CH) || (id < 0)) iFailed |= (1 << method_error);
	auto dScale(log(LAP::sdShigh/LAP::sdSlow)/ciNpts);
	vector<double> vTemp(iEventlength,0.);
	try {
		dTrace.reserve(iEventlength);
		dExp.reserve(ciNpts);
		dS.reserve(ciNpts);
		dXform.reserve(ciNpts);
	} catch (bad_alloc& ba) {
		iFailed |= (1 << alloc_error);
		return;
	}
	for (auto n = 0; n < ciNpts; n++) {
		dS[n] = LAP::sdSlow*exp(dScale*n);
		for (auto i = 0; i < iEventlength; i++) vTemp[i] = (exp(-dS[n]*i*dScaleT));
		dExp.push_back(vTemp);
	}
	iAve = 4; // 9 point average
	dAveScale = 1./(2.*iAve+1.);
}

LAP::~LAP() {
	if (g_verbose) cout << "LAP " << id << " d'tor\n";

	LAP::siHowMany--;
}

void LAP::root_init(TTree* tree_in) {
	if (!LAP::sbInitialized) {
		LAP::tree = unique_ptr<TTree>(tree_in);

		LAP::tree->Branch("LapLow",		LAP::sdLaplaceLow,	"laplo[8]/D");
		LAP::tree->Branch("LapHigh",	LAP::sdLaplaceHigh,	"laphi[8]/D");

		LAP::sbInitialized = true;
	}
}

void LAP::SetEvent(shared_ptr<Event> ev) {
	event = ev;
	if (event->GetAverage() > 0) iAve = 0;
}

void LAP::Analyze() { // not really optimized at all
	auto peak_x = event->Peak_x();
	LAP::sdLaplaceLow[id] = 0;
	LAP::sdLaplaceHigh[id] = 0;

	if (iAve) {
		for (auto t = peak_x; t < iEventlength-iAve; t++) { // performs the 9pt moving average
			dTrace[t] = 0;
			for (auto tt = -iAve; tt <= iAve; tt++) dTrace[t] += event->Trace(t+tt);
			dTrace[t] *= dAveScale;
		}
	}
	for (auto n = 0; n < ciNpts; n++) {
		dXform[n] = 0;
		for (auto t = peak_x; t < iEventlength-iAve; t++) dXform[n] += dExp[n][t - peak_x]*dTrace[t]; // trailing edge of the pulse
	}
	for (auto n = 0; n < ciNpts-1; n++) {
		(n < (ciNpts >> 1) ? LAP::sdLaplaceLow[id] : LAP::sdLaplaceHigh[id]) += (dS[n+1]-dS[n])*(dXform[n+1]+dXform[n]);
	}
	LAP::sdLaplaceHigh[id] *= 0.5*dScaleV;
	LAP::sdLaplaceLow[id] *= 0.5*dScaleV;
}
