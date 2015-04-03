#include "DFT.h"
#include <cstdlib>
#include <cmath>
#include <iostream>

float DFT::sfVersion = 1.41;
bool DFT::sbInitialized = false;
unique_ptr<TTree> DFT::tree = nullptr;
int DFT::siHowMany = 0;
const long double pi = 3.14159265358979l;

double DFT::sdMagnitude[8][4]	= { {0,0,0,0}, // 8 channels
									{0,0,0,0}, // 4 orders
									{0,0,0,0},
									{0,0,0,0},
									{0,0,0,0},
									{0,0,0,0},
									{0,0,0,0},
									{0,0,0,0}};
double DFT::sdPhase[8][4]		= {	{0,0,0,0},
									{0,0,0,0},
									{0,0,0,0},
									{0,0,0,0},
									{0,0,0,0},
									{0,0,0,0},
									{0,0,0,0},
									{0,0,0,0}};

DFT::DFT() : ciOrder(3) {
	if (g_verbose) cout << "DFT c'tor\n";
	DFT::siHowMany++;
}

DFT::DFT(int ch, int length, shared_ptr<Digitizer> digitizer) : Method(ch, length, digitizer), ciOrder(3) {
	if (g_verbose) cout << "DFT " << id << " c'tor\n";
	DFT::siHowMany++;
	if ((id >= MAX_CH) || (id < 0)) iFailed |= method_error;
	const int used_orders[] = {3,4,5}; // 4 orders incl 0th
	double omega;
	try {
		dCos = unique_ptr<double[]>(new double[ciOrder*iEventlength]); // simpler than two-dimensional arrays
		dSin = unique_ptr<double[]>(new double[ciOrder*iEventlength]);
	} catch (bad_alloc& ba) {iFailed |= alloc_error; return;}
	for (auto n = 0; n < ciOrder; n++) {
		omega = used_orders[n]*pi/(iEventlength*digitizer->ScaleT()); // GHz
		for (auto t = 0; t < iEventlength; t++) {
			dCos[n*iEventlength+t] = cos(omega*t);
			dSin[n*iEventlength+t] = sin(omega*t); // simpler than using one table and sin(x) = cos(x-pi/2)
	}	}
	dScalefactor = 2./iEventlength; // 2/period
}

DFT::~DFT() {
	if (g_verbose) cout << "DFT " << id << " d'tor\n";
	DFT::siHowMany--;
	dCos.reset();
	dSin.reset();
}

void DFT::root_init(TTree* tree_in) {
	if (!DFT::sbInitialized) {
		DFT::tree = unique_ptr<TTree>(tree_in);
		
		DFT::tree->Branch("Amplitude",	DFT::sdMagnitude,	"mag[8][4]/D");
		DFT::tree->Branch("Phase",		DFT::sdPhase,		"phase[8][4]/D"); // 4 orders
		
		DFT::sbInitialized = true;
	}
}

void DFT::Analyze() {
	auto dReal(0.), dImag(0.);
	auto n(0), t(0), nt(0);

	DFT::sdMagnitude[id][0] = 0; // phase[][0] never accessed so no need to reset it

	for (t = 0; t < iEventlength; t++) DFT::sdMagnitude[id][0] += event->Trace(t);
	DFT::sdMagnitude[id][0] *= 0.5*dScalefactor;
	for (n = 0; n < ciOrder; n++) {
		dReal = 0;
		dImag = 0;
		for (t = 0; t < iEventlength; t++) { // fourier series
			nt = n*iEventlength + t;
			dReal += event->Trace(t)*dCos[nt];
			dImag += event->Trace(t)*dSin[nt];
		}
		DFT::sdMagnitude[id][n+1] = dScalefactor*sqrt(dReal*dReal + dImag*dImag);
		DFT::sdPhase[id][n+1] = atan2(dImag,dReal);
	}
}


