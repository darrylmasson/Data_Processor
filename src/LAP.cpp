#include "LAP.h"
#include <cmath>
#include <iostream>

float LAP::si_version = 1.0;
float LAP::sf_s = 1.5;
bool LAP::sb_initialized = false;
unique_ptr<TTree> LAP::tree = nullptr;
int LAP::si_howmany = 0;

double LAP::laplace[8] = {0,0,0,0,0,0,0,0};

LAP::LAP(int ch, int len, shared_ptr<Digitizer> digitizer) {
	failed = 0;
	id = ch;
	eventlength = len;
	d_scale_v = digitizer->ScaleV();
	try {EXP = unique_ptr<double[]>(new double[eventlength]);}
	catch (bad_alloc& ba) {
		failed |= alloc_error;
		return;
	}
	for (int i = 0; i < eventlength; i++) d_EXP[i] = exp(-LAP::sf_s*digitizer->ScaleT()*i);
	LAP::si_howmany++;
}

LAP::~LAP() {
	if (g_verbose) cout << " LAP " << id << " d'tor ";
	d_EXP.reset();
	LAP::si_howmany--;
}

LAP::root_init(TTree* tree_in) {
	if (!LAP::sb_initialized) {
		LAP::tree = unique_ptr<TTree>(tree_in);
		
		LAP::tree->Branch("Laplace", LAP::laplace, "lap[8]/D");
		
		LAP::sb_initialized = true;
	}
}

LAP::Evaluate(const shared_ptr<Event> event) {
	int t(0);
	LAP::laplace[id] = 0;
	for (t = event->peak_x; t < eventlength; t++) {
		LAP::laplace[id] += d_EXP[t - event->peak_x]*event->trace[t];
	}
	LAP::laplace[id] *= d_scale_v;
}