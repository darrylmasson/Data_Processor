#include "LAP.h"
#include <cmath>
#include <iostream>

float LAP::sfVersion = 1.2;
float LAP::sfS = 1.5;
bool LAP::sbInitialized = false;
unique_ptr<TTree> LAP::tree = nullptr;
int LAP::siHowMany = 0;

double LAP::sdLaplace[8] = {0,0,0,0,0,0,0,0};
double LAP::sdLongInt[8] = {0,0,0,0,0,0,0,0};

LAP::LAP() {
	if (g_verbose) cout << "LAP c'tor\n";
	LAP::siHowMany++;
}

LAP::LAP(int ch, int length, shared_ptr<Digitizer> digitizer) : Method(ch, length, digitizer) {
	if (g_verbose) cout << "LAP " << id << " c'tor\n";
	LAP::siHowMany++;
	if ((id >= MAX_CH) || (id < 0)) iFailed |= method_error;
	try {dExp = unique_ptr<double[]>(new double[iEventlength]);}
	catch (bad_alloc& ba) {
		iFailed |= alloc_error;
		return;
	}
	for (auto i = 0; i < iEventlength; i++) dExp[i] = exp(-LAP::sfS*i/dScaleT); // not sure if the time scale is correct here.
}

LAP::~LAP() {
	if (g_verbose) cout << "LAP " << id << " d'tor\n";
	dExp.reset();
	LAP::siHowMany--;
}

void LAP::root_init(TTree* tree_in) {
	if (!LAP::sbInitialized) {
		LAP::tree = unique_ptr<TTree>(tree_in);
		
		LAP::tree->Branch("Laplace", LAP::sdLaplace, "lap[8]/D");
		LAP::tree->Branch("LongInt", LAP::sdLongInt, "longint[8]/D");
		
		LAP::sbInitialized = true;
	}
}

void LAP::Analyze(const shared_ptr<Event> event) {
	auto peak_x = event->Peak_x();
	LAP::sdLaplace[id] = 0;
	LAP::sdLongInt[id] = 0;
	for (auto t = 0; t < iEventlength; t++) {
		sdLongInt[id] += event->Trace(t);
		if (t >= peak_x) LAP::sdLaplace[id] += dExp[t - peak_x]*event->Trace(t); // trailing edge of the pulse
	}
	LAP::sdLongInt[id] = 2*LAP::sdLongInt[id] - (event->Trace(0) + event->Trace(iEventlength-1)); // trapezoid rule
	LAP::sdLongInt[id] = (event->Baseline()*iEventlength - 0.5*LAP::sdLongInt[id])*dScaleV*dScaleT;
	LAP::sdLaplace[id] *= dScaleV;
}
