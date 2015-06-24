#include "DFT.h"

float DFT::sfVersion = 1.5;
bool DFT::sbInitialized = false;
unique_ptr<TTree> DFT::tree = nullptr;
int DFT::siHowMany = 0;
const double pi = 3.14159265358979;

double DFT::sdMagnitude[8][ciOrder]	= { {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, // 8 channels
										{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, // 17 orders
										{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
										{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
										{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
										{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
										{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
										{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}};
double DFT::sdPhase[8][ciOrder]		= {	{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
										{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
										{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
										{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
										{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
										{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
										{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
										{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}};
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
		omega = (n+1)*pi/(iEventlength*dScaleT); // GHz
		for (auto t = 0; t < iEventlength; t++) {
			dCos[n].push_back(cos(omega*t));
			dSin[n].push_back(sin(omega*t)); // simpler than using one table and sin(x) = cos(x-pi/2)
	}	}
	dScalefactor = dScaleT*2./iEventlength; // 2/period, also in GHz
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
		
		DFT::tree->Branch("Amplitude",	DFT::sdMagnitude,	"mag[8][17]/D");
		DFT::tree->Branch("Phase",		DFT::sdPhase,		"phase[8][17]/D");
		DFT::tree->Branch("Even",		DFT::sdEven,		"even[8]/D");
		DFT::tree->Branch("Odd",		DFT::sdOdd,			"odd[8]/D");
		
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
		DFT::sdMagnitude[id][n] = dScalefactor*sqrt(dReal*dReal + dImag*dImag);
		DFT::sdPhase[id][n] = atan2(dImag,dReal);
		DFT::sdEven[id] = sdOdd[id] = 0;
		for (t = 2; t < ciOrder; t++) (t%2 ? sdEven[id] : sdOdd[id]) += (sdMagnitude[id][t]-sdMagnitude[id][(t%2?0:1)])/sdMagnitude[id][(t%2?0:1)];
	}
}


