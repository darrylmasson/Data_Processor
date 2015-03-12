#include "LAP.h"
#include <cmath>
#include <iostream>

float LAP::sf_version = 1.1;
float LAP::sf_s = 1.5;
bool LAP::sb_initialized = false;
unique_ptr<TTree> LAP::tree = nullptr;
int LAP::si_howmany = 0;

double LAP::sd_laplace[8] = {0,0,0,0,0,0,0,0};
double LAP::sd_longint[8] = {0,0,0,0,0,0,0,0};

LAP::LAP(int ch, int len, shared_ptr<Digitizer> digitizer) {
	failed = 0;
	id = ch;
	eventlength = len;
	try {d_EXP = unique_ptr<double[]>(new double[eventlength]);}
	catch (bad_alloc& ba) {
		failed |= alloc_error;
		return;
	}
	d_scaleV = digitizer->ScaleV();
	d_scaleT = digitizer->ScaleT();
	for (auto i = 0; i < eventlength; i++) d_EXP[i] = exp(-LAP::sf_s*d_scaleT*i);
	LAP::si_howmany++;
}

LAP::~LAP() {
	if (g_verbose) cout << " LAP " << id << " d'tor ";
	d_EXP.reset();
	LAP::si_howmany--;
}

void LAP::root_init(TTree* tree_in) {
	if (!LAP::sb_initialized) {
		LAP::tree = unique_ptr<TTree>(tree_in);
		
		LAP::tree->Branch("Laplace", LAP::sd_laplace, "lap[8]/D");
		LAP::tree->Branch("LongInt", LAP::sd_longint, "longint[8]/D");
		
		LAP::sb_initialized = true;
	}
}

void LAP::evaluate(const shared_ptr<Event> event) {
	auto peak_x = event->Peak_x();
	LAP::sd_laplace[id] = 0;
	LAP::sd_longint[id] = 0;
	for (auto t = 0; t < eventlength; t++) {
		sd_longint[id] += event->Trace(t);
		if (t >= peak_x) LAP::sd_laplace[id] += d_EXP[t - peak_x]*event->Trace(t); // trailing edge of the pulse
	}
	LAP::sd_longint[id] = 2*LAP::sd_longint[id] - (event->Trace(0) + event->Trace(eventlength-1)); // trapezoid rule
	LAP::sd_longint[id] = (event->Baseline()*eventlength - 0.5*LAP::sd_longint[id])*d_scaleV*d_scaleT;
	LAP::sd_laplace[id] *= d_scaleV;
}
