#include "XSQ.h"
#include <cstdlib>
#include <algorithm>
#include "TVectorT.h"
#include <iostream> // remove

float	XSQ::sf_version = 1.1;
bool	XSQ::sb_initialized = false;
int		XSQ::si_howmany = 0;

unique_ptr<TTree> XSQ::tree = nullptr;

double XSQ::sd_xsq_n[4]			= {0,0,0,0};
double XSQ::sd_peakheight_n[4]	= {0,0,0,0};
double XSQ::sd_baseline_n[4]	= {0,0,0,0};
double XSQ::sd_offset_n[4]		= {0,0,0,0};
double XSQ::sd_peak_err_n[4]	= {0,0,0,0};
double XSQ::sd_base_err_n[4]	= {0,0,0,0};
double XSQ::sd_offset_err_n[4]	= {0,0,0,0};
double XSQ::sd_prob_n[4]		= {0,0,0,0};

double XSQ::sd_xsq_y[4]			= {0,0,0,0};
double XSQ::sd_peakheight_y[4]	= {0,0,0,0};
double XSQ::sd_baseline_y[4]	= {0,0,0,0};
double XSQ::sd_offset_y[4]		= {0,0,0,0};
double XSQ::sd_peak_err_y[4]	= {0,0,0,0};
double XSQ::sd_base_err_y[4]	= {0,0,0,0};
double XSQ::sd_offset_err_y[4]	= {0,0,0,0};
double XSQ::sd_prob_y[4]		= {0,0,0,0};

int XSQ::si_fit_status_n[4]		= {0,0,0,0};
int XSQ::si_fit_status_y[4]		= {0,0,0,0};

XSQ::XSQ(const int ch, const int len, const float gain_in[], const shared_ptr<Digitizer> digitizer) : ci_nPar(4) {
	failed = 0;
	id = ch;
	eventlength = len;
	f_gain[n] = gain_in[n];
	f_gain[y] = gain_in[y];
	XSQ::si_howmany++;
	char filename[64];
	unique_ptr<TFile> std_file = nullptr;
	TVectorT<double>* wave = nullptr;
	int i(0), p(0);
	if ((id > 3) || (id < 0)) failed |= method_error;
	switch (digitizer->ID()) {
		case dt5751 : 
			if (digitizer->Special() == 0) {
				i_std_length	= 225;
				i_std_trig		= 32;
				i_std_peak_x	= 35;
			} else {
				i_std_length	= 450;
				i_std_trig		= 64;
				i_std_peak_x	= 71;
			} break;
		case dt5751des :
			i_std_length		= 899;
			i_std_trig			= 128;
			i_std_peak_x		= 142;
			break;
		case dt5730 :
			i_std_length		= 225;
			i_std_trig			= 32;
			i_std_peak_x		= 35;
			break;
		case v1724 :
		default :
			failed |= dig_error;
			return;
	}
	
	sprintf(filename, "%sconfig/standard_events.root", path);
	try {std_file.reset(new TFile(filename, "READ"));}
	catch (bad_alloc& ba) {failed |= alloc_error; return;}
	if (!std_file->IsOpen()) {failed |= file_error; return;}
	for (p = 0; p < P; p++) {
		try {d_std_wave[p] = unique_ptr<double[]>(new double[i_std_length]);}
		catch (bad_alloc& ba) {failed |= alloc_error; return;}
		wave = (TVectorT<double>*)std_file->Get((p ? "gamma_wave" : "neutron_wave"));
		if (wave == nullptr) {failed |= root_error; return;}
		d_std_base[p] = 0;
		d_std_peak[p] = 1000;
		switch(i_std_length) {
			case 225 : // 500 MSa/s
				for (i = 0; i < i_std_length; i++) {
					d_std_wave[p][i] = ((*wave)[2*i] + (*wave)[2*i+1])/2.;
					if (i < digitizer->Baselength()) d_std_base[p] += d_std_wave[p][i];
					d_std_peak[p] = min(d_std_peak[p], d_std_wave[p][i]);
				} break;
			case 899 : // 2 GSa/s
				for (i = 0; i < i_std_length; i++) {
					d_std_wave[p][i] = (i%2) ? ((*wave)[(i+1)/2] + (*wave)[(i-1)/2])/2. : (*wave)[i/2];
					if (i < digitizer->Baselength()) d_std_base[p] += d_std_wave[p][i];
					d_std_peak[p] = min(d_std_peak[p], d_std_wave[p][i]);
				} break;
			case 450 : // 1 GSa/s
				for (i = 0; i < i_std_length; i++) {
					d_std_wave[p][i] = (*wave)[i];
					if (i < digitizer->Baselength()) d_std_base[p] += d_std_wave[p][i];
					d_std_peak[p] = min(d_std_peak[p], d_std_wave[p][i]);
				} break;
			default : failed |= method_error;
			return;
		}
		d_std_base[p] /= digitizer->Baselength();
		d_std_norm[p] = 1./(d_std_base[p]-d_std_peak[p]);
	} // p
	std_file->Close();
	std_file = nullptr;
	wave = nullptr;
	try {
		i_input_wave	= unique_ptr<int[]>(new int[eventlength]);
		i_x				= unique_ptr<int[]>(new int[eventlength]);
		
		d_pars			= unique_ptr<double[]>(new double[ci_nPar]);
		
		fit				= unique_ptr<TF1>(new TF1("fit", this, &XSQ::fitter, 0, min(eventlength, i_std_length), ci_nPar));
	} catch (bad_alloc& ba) {failed |= alloc_error; return;}
	if (fit->IsZombie()) {failed |= root_error; return;}
	for (i = 0; i < eventlength; i++) i_x[i] = i;
	
	fit->SetParNames("Peakheight_scale","Baseline_offset","Trigger_offset","Particle");

	graph = nullptr;
}

XSQ::~XSQ() {
	if (g_verbose) cout << " XSQ " << id << " d'tor ";
	XSQ::si_howmany--;
	for (int p = 0; p < P; p++) {
		d_std_wave[p].reset();
	}
	d_pars.reset();
	fit.reset();
	graph.reset();
	i_input_wave.reset();
	i_x.reset();
}

void XSQ::root_init(TTree* tree_in) {
	if (!XSQ::sb_initialized) {
		XSQ::tree = unique_ptr<TTree>(tree_in);
		
		XSQ::tree->Branch("Chisquare_n",	XSQ::sd_xsq_n,			"xsqn[4]/D");
		XSQ::tree->Branch("Peakscale_n",	XSQ::sd_peakheight_n,	"pkscalen[4]/D");
		XSQ::tree->Branch("Base_shift_n",	XSQ::sd_baseline_n,		"baseshiftn[4]/D");
		XSQ::tree->Branch("Offset_n",		XSQ::sd_offset_n,		"offsetn[4]/D");
		XSQ::tree->Branch("Prob_n",			XSQ::sd_prob_n,			"probn[4]/D");
		
		XSQ::tree->Branch("Chisquare_y",	XSQ::sd_xsq_y,			"xsqy[4]/D");
		XSQ::tree->Branch("Peakscale_y",	XSQ::sd_peakheight_y,	"pkscaley[4]/D");
		XSQ::tree->Branch("Base_shift_y",	XSQ::sd_baseline_y,		"baseshifty[4]/D");
		XSQ::tree->Branch("Offset_y",		XSQ::sd_offset_y,		"offsety[4]/D");
		XSQ::tree->Branch("Prob_y",			XSQ::sd_prob_y,			"proby[4]/D");
		
		XSQ::tree->Branch("Peak_err_n",		XSQ::sd_peak_err_n,		"pkerrn[4]/D");
		XSQ::tree->Branch("Base_err_n",		XSQ::sd_base_err_n,		"berrn[4]/D");
		XSQ::tree->Branch("Offset_err_n",	XSQ::sd_offset_err_n,	"offerrn[4]/D");
		XSQ::tree->Branch("Peak_err_y",		XSQ::sd_peak_err_y,		"pkerry[4]/D");
		XSQ::tree->Branch("Base_err_y",		XSQ::sd_base_err_y,		"berry[4]/D");
		XSQ::tree->Branch("Offset_err_y",	XSQ::sd_offset_err_y,	"offerry[4]/D");
		
		XSQ::tree->Branch("Fit_status_n",	XSQ::si_fit_status_n,	"fitstatusn[4]/i");
		XSQ::tree->Branch("Fit_status_y",	XSQ::si_fit_status_y,	"fitstatusy[4]/i");
		
		XSQ::sb_initialized = true;
	}
}

double XSQ::fitter(double* x, double* par) {
	int i_sample = (x[0] - par[2]), p = par[3];
	double d_val = d_std_base[p] + par[1];
	if ((i_sample > -1) && (i_sample < i_std_length)) d_val += par[0]*(d_std_wave[p][i_sample]-d_std_base[p]);
	if (d_val < 0) d_val = 0; // saturated events
	return d_val;
}

void XSQ::evaluate(const shared_ptr<Event> event) {
	for (int i = 0; i < eventlength; i++) i_input_wave[i] = event->trace[i];
	try {graph.reset(new TGraph(eventlength, i_x.get(), i_input_wave.get()));}
	catch (bad_alloc& ba) { // error codes don't work here
		XSQ::sd_xsq_n[id]			= -1;
		XSQ::sd_peakheight_n[id]	= -1;
		XSQ::sd_baseline_n[id]		= -1;
		XSQ::sd_offset_n[id]		= -1;
		XSQ::sd_prob_n[id]			= -1;
		XSQ::sd_peak_err_n[id]		= -1;
		XSQ::sd_base_err_n[id]		= -1;
		XSQ::sd_offset_err_n[id]	= -1;
		
		XSQ::sd_xsq_y[id]			= -1;
		XSQ::sd_peakheight_y[id]	= -1;
		XSQ::sd_baseline_y[id]		= -1;
		XSQ::sd_offset_y[id]		= -1;
		XSQ::sd_prob_y[id]			= -1;
		XSQ::sd_peak_err_y[id]		= -1;
		XSQ::sd_base_err_y[id]		= -1;
		XSQ::sd_offset_err_y[id]	= -1;
		
		XSQ::si_fit_status_n[id]	= -1;
		XSQ::si_fit_status_y[id]	= -1;
		
		return;
	}
	
	d_pars[0] = (event->baseline - event->peak_y)*d_std_norm[n];
	d_pars[1] = event->baseline - d_std_base[n];
	d_pars[2] = event->trigger - i_std_trig;
	d_pars[3] = n;
	fit->SetParameters(d_pars.get());
	fit->FixParameter(3, n);
	
	XSQ::si_fit_status_n[id]	= graph->Fit(fit.get(), "Q N R"); // quiet, no-plot, specified range

	XSQ::sd_xsq_n[id]			= fit->GetChisquare();
	XSQ::sd_peakheight_n[id]	= fit->GetParameter(0)*f_gain[n];
	XSQ::sd_baseline_n[id]		= fit->GetParameter(1);
	XSQ::sd_offset_n[id]		= fit->GetParameter(2);
	XSQ::sd_prob_n[id]			= fit->GetProb();

	XSQ::sd_peak_err_n[id]		= fit->GetParError(0)*f_gain[n];
	XSQ::sd_base_err_n[id]		= fit->GetParError(1);
	XSQ::sd_offset_err_n[id]	= fit->GetParError(2);

	
	d_pars[0] = (event->baseline - event->peak_y)*d_std_norm[y];
	d_pars[1] = event->baseline - d_std_base[y];
	d_pars[2] = event->trigger - i_std_trig;
	d_pars[3] = y;
	fit->SetParameters(d_pars.get());
	fit->FixParameter(3, y);
	
	XSQ::si_fit_status_y[id]	= graph->Fit(fit.get(), "Q N R");

	XSQ::sd_xsq_y[id]			= fit->GetChisquare();
	XSQ::sd_peakheight_y[id]	= fit->GetParameter(0)*f_gain[y];
	XSQ::sd_baseline_y[id]		= fit->GetParameter(1);
	XSQ::sd_offset_y[id]		= fit->GetParameter(2);
	XSQ::sd_prob_y[id]			= fit->GetProb();

	XSQ::sd_peak_err_y[id]		= fit->GetParError(0)*f_gain[y];
	XSQ::sd_base_err_y[id]		= fit->GetParError(1);
	XSQ::sd_offset_err_y[id]	= fit->GetParError(2);
}
