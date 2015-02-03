#include "DFT.h"
#include <cstdlib>
#include <cmath>

float DFT::version = 1.40;
bool DFT::initialized = false;
weak_ptr<TTree> DFT::tree = nullptr;
int DFT::howmany = 0;

double DFT::magnitude[8][4]	= {	{0,0,0,0}, // 8 channels
								{0,0,0,0}, // 4 orders
								{0,0,0,0},
								{0,0,0,0},
								{0,0,0,0},
								{0,0,0,0},
								{0,0,0,0},
								{0,0,0,0}};
double DFT::phase[8][4]		= {	{0,0,0,0},
								{0,0,0,0},
								{0,0,0,0},
								{0,0,0,0},
								{0,0,0,0},
								{0,0,0,0},
								{0,0,0,0},
								{0,0,0,0}};

DFT::DFT(const int ch, const int len) : order(3) {
	eventlength = len;
	failed = 0;
	id = ch;
	if ((id >= MAX_CH) || (id < 0)) failed |= method_error;
	const int used_orders[] = {3,4,5}; // 4 orders incl 0th
	double omega;
	try {
		COS = unique_ptr<double[]>(new double[order*eventlength]);
		SIN = unique_ptr<double[]>(new double[order*eventlength]);
	} catch (bad_alloc& ba) {failed |= alloc_error; return;}
	for (int n = 0; n < order; n++) {
		omega = used_orders[n]*pi/eventlength; // GHz
		for (int t = 0; t < eventlength; t++) {
			COS[n*eventlength+t] = cos(omega*t);
			SIN[n*eventlength+t] = sin(omega*t);
	}	}
	DFT::howmany++;
	scalefactor = 2./eventlength; // 2/period
}

DFT::~DFT() {
	DFT::howmany--;
	COS.reset();
	SIN.reset();
	DFT::tree.reset();
}

void DFT::root_init(shared_ptr<TTree> tree_in) {
	if (!DFT::initialized) {
		DFT::tree = tree_in;
		
		DFT::tree->Branch("Amplitude",	DFT::magnitude,	"mag[8][4]/D");
		DFT::tree->Branch("Phase",		DFT::phase,		"phase[8][4]/D"); // 4 orders
		
		DFT::initialized = true;
	}
}

void DFT::evaluate(const weak_ptr<Event> event) {
	double re(0), im(0);
	int n(0), t(0), nt(0);

	DFT::magnitude[id][0] = 0;
	DFT::phase[id][0] = 0;

	for (t = 0; t < eventlength; t++) DFT::magnitude[id][0] += event->trace[t];
	DFT::magnitude[id][0] *= scalefactor/2;
	for (n = 0; n < order; n++) {
		re = 0;
		im = 0;
		for (t = 0; t < eventlength; t++) {
			nt = n*eventlength + t;
			re += event->trace[t]*COS[nt];
			im += event->trace[t]*SIN[nt];
		}
		DFT::magnitude[id][n+1] = scalefactor*sqrt(re*re + im*im);
		DFT::phase[id][n+1] = atan2(im,re);
	}
}


