#include "DFT.h"
#include <cstdlib>
#include <cmath>
#include <iostream>

float DFT::si_version = 1.41;
bool DFT::sb_initialized = false;
unique_ptr<TTree> DFT::tree = nullptr;
int DFT::si_howmany = 0;

double DFT::sd_magnitude[8][4]	= {	{0,0,0,0}, // 8 channels
								{0,0,0,0}, // 4 orders
								{0,0,0,0},
								{0,0,0,0},
								{0,0,0,0},
								{0,0,0,0},
								{0,0,0,0},
								{0,0,0,0}};
double DFT::sd_phase[8][4]		= {	{0,0,0,0},
								{0,0,0,0},
								{0,0,0,0},
								{0,0,0,0},
								{0,0,0,0},
								{0,0,0,0},
								{0,0,0,0},
								{0,0,0,0}};

DFT::DFT(const int ch, const int len, const shared_ptr<Digitizer> digitizer) : ci_order(3) {
	eventlength = len;
	failed = 0;
	id = ch;
	if ((id >= MAX_CH) || (id < 0)) failed |= method_error;
	const int used_orders[] = {3,4,5}; // 4 orders incl 0th
	double omega;
	try {
		d_COS = unique_ptr<double[]>(new double[order*eventlength]);
		d_SIN = unique_ptr<double[]>(new double[order*eventlength]);
	} catch (bad_alloc& ba) {failed |= alloc_error; return;}
	for (int n = 0; n < ci_order; n++) {
		omega = used_orders[n]*pi/(eventlength*digitizer->ScaleT()); // GHz
		for (int t = 0; t < eventlength; t++) {
			d_COS[n*eventlength+t] = cos(omega*t);
			d_SIN[n*eventlength+t] = sin(omega*t);
	}	}
	DFT::si_howmany++;
	d_scalefactor = 2./eventlength; // 2/period
}

DFT::~DFT() {
	if (g_verbose) cout << " DFT " << id << " d'tor ";
	DFT::si_howmany--;
	d_COS.reset();
	d_SIN.reset();
}

void DFT::root_init(TTree* tree_in) {
	if (!DFT::sb_initialized) {
		DFT::tree = unique_ptr<TTree>(tree_in);
		
		DFT::tree->Branch("Amplitude",	DFT::sd_magnitude,	"mag[8][4]/D");
		DFT::tree->Branch("Phase",		DFT::sd_phase,		"phase[8][4]/D"); // 4 orders
		
		DFT::sb_initialized = true;
	}
}

void DFT::evaluate(const shared_ptr<Event> event) {
	double d_re(0), d_im(0);
	int n(0), t(0), nt(0);

	DFT::sd_magnitude[id][0] = 0;
	DFT::sd_phase[id][0] = 0;

	for (t = 0; t < eventlength; t++) DFT::sd_magnitude[id][0] += event->trace[t];
	DFT::sd_magnitude[id][0] *= d_scalefactor/2;
	for (n = 0; n < order; n++) {
		d_re = 0;
		d_im = 0;
		for (t = 0; t < eventlength; t++) {
			nt = n*eventlength + t;
			d_re += event->trace[t]*d_COS[nt];
			d_im += event->trace[t]*d_SIN[nt];
		}
		DFT::sd_magnitude[id][n+1] = d_scalefactor*sqrt(d_re*d_re + d_im*d_im);
		DFT::sd_phase[id][n+1] = atan2(d_im,d_re);
	}
}


