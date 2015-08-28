#include "DFT.h"

float DFT::sfVersion = 1.55;
bool DFT::sbInitialized = false;
unique_ptr<TTree> DFT::tree = nullptr;
int DFT::siHowMany = 0;
const double pi = 3.14159265358979;

double DFT::sdOdd[8]	= {0,0,0,0,0,0,0,0};
double DFT::sdEven[8]	= {0,0,0,0,0,0,0,0};

DFT::DFT() {
	if (g_verbose) cout << "DFT c'tor\n";
	DFT::siHowMany++;
}

DFT::DFT(int ch, int length, shared_ptr<Digitizer> digitizer) : Method(ch, length, digitizer) {
	if (g_verbose) cout << "DFT " << id << " c'tor\n";
	DFT::siHowMany++;
	if ((id >= MAX_CH) || (id < 0)) iFailed |= (1 << method_error);
	double omega;
	for (auto n = 0; n < ciOrder; n++) {
		try {
			dCos[n].reserve(iEventlength);
			dSin[n].reserve(iEventlength);
		} catch (bad_alloc& ba) {iFailed |= (1 << alloc_error); return;}
		omega = 2*n*pi/(iEventlength*dScaleT); // GHz
		for (auto t = 0; t < iEventlength; t++) {
			dCos[n].push_back(cos(omega*t));
			dSin[n].push_back(sin(omega*t)); // simpler than using one table and sin(x) = cos(x-pi/2)
	}	}
}

DFT::~DFT() {
	if (g_verbose) cout << "DFT " << id << " d'tor\n";
	DFT::siHowMany--;
	for (auto i = 0; i < ciOrder; i++) {
		dCos[i].clear();
		dSin[i].clear();
	}
}

void DFT::root_init(TTree* tree_in) {
	if (!DFT::sbInitialized) {
		DFT::tree = unique_ptr<TTree>(tree_in);

		DFT::tree->Branch("Even",	DFT::sdEven,	"even[8]/D");
		DFT::tree->Branch("Odd",	DFT::sdOdd,		"odd[8]/D");

		DFT::sbInitialized = true;
	}
}

void DFT::Analyze() {
	auto dReal(0.), dImag(0.);
	auto n(0), t(0);
	for (n = 0; n < ciOrder; n++) {
		dReal = 0;
		dImag = 0;
		for (t = 0; t < iEventlength; t++) { // fourier series
			dReal += event->Trace(t)*dCos[n][t];
			dImag += event->Trace(t)*dSin[n][t];
		}
		dMagnitude[n] = sqrt(dReal*dReal + dImag*dImag);
	}
	DFT::sdEven[id] = DFT::sdOdd[id] = 0;
	for (t = 2; t < ciOrder; t++) (t%2 ? DFT::sdOdd[id] : DFT::sdEven[id]) += (dMagnitude[t]-dMagnitude[(t%2?1:0)])/dMagnitude[(t%2?1:0)];
}


