#include "DFT.h"
#include <cstdlib>
#include <cmath>
#include <iostream>

float DFT::sf_version = 1.41;
bool DFT::sb_initialized = false;
unique_ptr<TTree> DFT::tree = nullptr;
int DFT::si_howmany = 0;
const long double pi = 3.14159265358979l;

double DFT::sd_magnitude[8][4]	= { {0,0,0,0}, // 8 channels
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
		d_COS = unique_ptr<double[]>(new double[ci_order*eventlength]); // simpler than two-dimensional arrays
		d_SIN = unique_ptr<double[]>(new double[ci_order*eventlength]);
	} catch (bad_alloc& ba) {failed |= alloc_error; return;}
	for (auto n = 0; n < ci_order; n++) {
		omega = used_orders[n]*pi/(eventlength*digitizer->ScaleT()); // GHz
		for (auto t = 0; t < eventlength; t++) {
			d_COS[n*eventlength+t] = cos(omega*t);
			d_SIN[n*eventlength+t] = sin(omega*t); // simpler than using one table and sin(x) = cos(x-pi/2)
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
	auto d_re(0.), d_im(0.);
	auto n(0), t(0), nt(0);

	DFT::sd_magnitude[id][0] = 0;
	DFT::sd_phase[id][0] = 0;

	for (t = 0; t < eventlength; t++) DFT::sd_magnitude[id][0] += event->Trace(t);
	DFT::sd_magnitude[id][0] *= d_scalefactor/2;
	for (n = 0; n < ci_order; n++) {
		d_re = 0;
		d_im = 0;
		for (t = 0; t < eventlength; t++) { // fourier series are pretty straightforward
			nt = n*eventlength + t;
			d_re += event->Trace(t)*d_COS[nt];
			d_im += event->Trace(t)*d_SIN[nt];
		}
		DFT::sd_magnitude[id][n+1] = d_scalefactor*sqrt(d_re*d_re + d_im*d_im);
		DFT::sd_phase[id][n+1] = atan2(d_im,d_re);
	}
}


