#include "XSQ.h"
#include <cstdlib>
#include <algorithm>
#include "TVectorT.h"
#include <iostream> // remove

float	XSQ::version = 1.1;
bool	XSQ::initialized = false;
int		XSQ::howmany = 0;

shared_ptr<TTree> XSQ::tree = nullptr;

double XSQ::xsq_n[4]		= {0,0,0,0};
double XSQ::peakheight_n[4]	= {0,0,0,0};
double XSQ::baseline_n[4]	= {0,0,0,0};
double XSQ::offset_n[4]		= {0,0,0,0};
double XSQ::peak_err_n[4]	= {0,0,0,0};
double XSQ::base_err_n[4]	= {0,0,0,0};
double XSQ::offset_err_n[4]	= {0,0,0,0};
double XSQ::prob_n[4]		= {0,0,0,0};

double XSQ::xsq_y[4]		= {0,0,0,0};
double XSQ::peakheight_y[4]	= {0,0,0,0};
double XSQ::baseline_y[4]	= {0,0,0,0};
double XSQ::offset_y[4]		= {0,0,0,0};
double XSQ::peak_err_y[4]	= {0,0,0,0};
double XSQ::base_err_y[4]	= {0,0,0,0};
double XSQ::offset_err_y[4]	= {0,0,0,0};
double XSQ::prob_y[4]		= {0,0,0,0};

int XSQ::fit_status_n[4]	= {0,0,0,0};
int XSQ::fit_status_y[4]	= {0,0,0,0};

XSQ::XSQ(const int ch, const int len, const float gain_in[], const shared_ptr<Digitizer> digitizer) : nPar(4) {
	failed = 0;
	id = ch;
	eventlength = len;
	gain[n] = gain_in[n];
	gain[y] = gain_in[y];
	XSQ::howmany++;
	char filename[64];
	unique_ptr<TFile> std_file = nullptr;
	TVectorT<double>* wave = nullptr;
	int i(0), p(0);
	if ((id > 3) || (id < 0)) failed |= method_error;
	switch (digitizer->ID()) {
		case dt5751 : 
			if (digitizer->Special() == 0) {
				std_length	= 225;
				std_trig	= 32;
				std_peak_x	= 35;
			} else {
				std_length	= 450;
				std_trig	= 64;
				std_peak_x	= 71;
			} break;
		case dt5751des :
			std_length		= 899;
			std_trig		= 128;
			std_peak_x		= 142;
			break;
		case dt5730 :
			std_length		= 225;
			std_trig		= 32;
			std_peak_x		= 35;
			break;
		case v1724 :
		default :
			failed |= dig_error;
			return;
	}
	
	sprintf(filename, "%sconfig/standard_events.root", path);
	std_file.reset(new TFile(filename, "READ"));
	if (!std_file->IsOpen()) {failed |= file_error; return;}
	for (p = 0; p < P; p++) {
		try {std_wave[p] = unique_ptr<double[]>(new double[std_length]);}
		catch (bad_alloc& ba) {failed |= alloc_error; return;}
		wave = (TVectorT<double>*)std_file->Get((p ? "gamma_wave" : "neutron_wave"));
		if (wave == nullptr) {failed |= root_error; return;}
		std_base[p] = 0;
		std_peak[p] = 1000;
		switch(std_length) {
			case 225 : // 500 MSa/s
				for (i = 0; i < std_length; i++) {
					std_wave[p][i] = ((*wave)[2*i] + (*wave)[2*i+1])/2.;
					if (i < digitizer->Baselength()) std_base[p] += std_wave[p][i];
					std_peak[p] = min(std_peak[p], std_wave[p][i]);
				} break;
			case 899 : // 2 GSa/s
				for (i = 0; i < std_length; i++) {
					std_wave[p][i] = (i%2) ? ((*wave)[(i+1)/2] + (*wave)[(i-1)/2])/2. : (*wave)[i/2];
					if (i < digitizer->Baselength()) std_base[p] += std_wave[p][i];
					std_peak[p] = min(std_peak[p], std_wave[p][i]);
				} break;
			case 450 : // 1 GSa/s
				for (i = 0; i < std_length; i++) {
					std_wave[p][i] = (*wave)[i];
					if (i < digitizer->Baselength()) std_base[p] += std_wave[p][i];
					std_peak[p] = min(std_peak[p], std_wave[p][i]);
				} break;
			default : failed |= method_error;
			return;
		}
		std_base[p] /= digitizer->Baselength();
		std_norm[p] = 1./(std_base[p]-std_peak[p]);
	} // p
	std_file->Close();
	std_file = nullptr;
	wave = nullptr;
	try {
		input_wave	= unique_ptr<int[]>(new int[eventlength]);
		x			= unique_ptr<int[]>(new int[eventlength]);
		
		pars		= unique_ptr<double[]>(new double[nPar]);
		
		fit			= unique_ptr<TF1>(new TF1("fit", this, &XSQ::fitter, 0, min(eventlength, std_length), nPar));
	} catch (bad_alloc& ba) {failed |= alloc_error; return;}
	if (fit->IsZombie()) {failed |= root_error; return;}
	for (i = 0; i < eventlength; i++) x[i] = i;
	
	fit->SetParNames("Peakheight_scale","Baseline_offset","Trigger_offset","Particle");

	graph = nullptr;
}

XSQ::~XSQ() {
	if (g_verbose) cout << " XSQ " << id << " d'tor ";
	XSQ::howmany--;
	for (int p = 0; p < P; p++) {
		std_wave[p].reset();
	}
	pars.reset();
	fit.reset();
	graph.reset();
	input_wave.reset();
	x.reset();
//	XSQ::tree = nullptr;
}

void XSQ::root_init(shared_ptr<TTree> tree_in) {
	if (!XSQ::initialized) {
		XSQ::tree = tree_in;
		
		XSQ::tree->Branch("Chisquare_n",	XSQ::xsq_n,			"xsqn[4]/D");
		XSQ::tree->Branch("Peakscale_n",	XSQ::peakheight_n,	"pkscalen[4]/D");
		XSQ::tree->Branch("Base_shift_n",	XSQ::baseline_n,	"baseshiftn[4]/D");
		XSQ::tree->Branch("Offset_n",		XSQ::offset_n,		"offsetn[4]/D");
		XSQ::tree->Branch("Prob_n",			XSQ::prob_n,		"probn[4]/D");
		
		XSQ::tree->Branch("Chisquare_y",	XSQ::xsq_y,			"xsqy[4]/D");
		XSQ::tree->Branch("Peakscale_y",	XSQ::peakheight_y,	"pkscaley[4]/D");
		XSQ::tree->Branch("Base_shift_y",	XSQ::baseline_y,	"baseshifty[4]/D");
		XSQ::tree->Branch("Offset_y",		XSQ::offset_y,		"offsety[4]/D");
		XSQ::tree->Branch("Prob_y",			XSQ::prob_y,		"proby[4]/D");
		
		XSQ::tree->Branch("Peak_err_n",		XSQ::peak_err_n,	"pkerrn[4]/D");
		XSQ::tree->Branch("Base_err_n",		XSQ::base_err_n,	"berrn[4]/D");
		XSQ::tree->Branch("Offset_err_n",	XSQ::offset_err_n,	"offerrn[4]/D");
		XSQ::tree->Branch("Peak_err_y",		XSQ::peak_err_y,	"pkerry[4]/D");
		XSQ::tree->Branch("Base_err_y",		XSQ::base_err_y,	"berry[4]/D");
		XSQ::tree->Branch("Offset_err_y",	XSQ::offset_err_y,	"offerry[4]/D");
		
		XSQ::tree->Branch("Fit_status_n",	XSQ::fit_status_n,	"fitstatusn[4]/i");
		XSQ::tree->Branch("Fit_status_y",	XSQ::fit_status_y,	"fitstatusy[4]/i");
		
		XSQ::initialized = true;
	}
}

double XSQ::fitter(double* x, double* par) {
	int sample = (x[0] - par[2]), p = par[3];
	double val = std_base[p] + par[1];
	if ((sample > -1) && (sample < std_length)) val += par[0]*(std_wave[p][sample]-std_base[p]);
	if (val < 0) val = 0; // saturated events
	return val;
}

void XSQ::evaluate(const shared_ptr<Event> event) {
	for (int i = 0; i < eventlength; i++) input_wave[i] = event->trace[i];
	try {graph.reset(new TGraph(eventlength, x.get(), input_wave.get()));}
	catch (bad_alloc& ba) { // error codes don't work here
		XSQ::xsq_n[id]			= -1;
		XSQ::peakheight_n[id]	= -1;
		XSQ::baseline_n[id]		= -1;
		XSQ::offset_n[id]		= -1;
		XSQ::prob_n[id]			= -1;
		XSQ::peak_err_n[id]		= -1;
		XSQ::base_err_n[id]		= -1;
		XSQ::offset_err_n[id]	= -1;
		
		XSQ::xsq_y[id]			= -1;
		XSQ::peakheight_y[id]	= -1;
		XSQ::baseline_y[id]		= -1;
		XSQ::offset_y[id]		= -1;
		XSQ::prob_y[id]			= -1;
		XSQ::peak_err_y[id]		= -1;
		XSQ::base_err_y[id]		= -1;
		XSQ::offset_err_y[id]	= -1;
		
		XSQ::fit_status_n[id]	= -1;
		XSQ::fit_status_y[id]	= -1;
		
		return;
	}
	
	pars[0] = (event->baseline - event->peak_y)*std_norm[n];
	pars[1] = event->baseline - std_base[n];
	pars[2] = event->trigger - std_trig;
	pars[3] = n;
	fit->SetParameters(pars.get());
	fit->FixParameter(3, n);
	
	XSQ::fit_status_n[id]	= graph->Fit(fit.get(), "Q N R"); // quiet, no-plot, specified range

	XSQ::xsq_n[id]			= fit->GetChisquare();
	XSQ::peakheight_n[id]	= fit->GetParameter(0)*gain[n];
	XSQ::baseline_n[id]		= fit->GetParameter(1);
	XSQ::offset_n[id]		= fit->GetParameter(2);
	XSQ::prob_n[id]			= fit->GetProb();

	XSQ::peak_err_n[id]		= fit->GetParError(0)*gain[n];
	XSQ::base_err_n[id]		= fit->GetParError(1);
	XSQ::offset_err_n[id]	= fit->GetParError(2);

	
	pars[0] = (event->baseline - event->peak_y)*std_norm[y];
	pars[1] = event->baseline - std_base[y];
	pars[2] = event->trigger - std_trig;
	pars[3] = y;
	fit->SetParameters(pars.get());
	fit->FixParameter(3, y);
	
	XSQ::fit_status_y[id]	= graph->Fit(fit.get(), "Q N R");

	XSQ::xsq_y[id]			= fit->GetChisquare();
	XSQ::peakheight_y[id]	= fit->GetParameter(0)*gain[y];
	XSQ::baseline_y[id]		= fit->GetParameter(1);
	XSQ::offset_y[id]		= fit->GetParameter(2);
	XSQ::prob_y[id]			= fit->GetProb();

	XSQ::peak_err_y[id]		= fit->GetParError(0)*gain[y];
	XSQ::base_err_y[id]		= fit->GetParError(1);
	XSQ::offset_err_y[id]	= fit->GetParError(2);
}
