#include "XSQ.h"
#include <cstdlib>
#include <algorithm>
#include "TVectorT.h"
#include <iostream> // remove

float	XSQ::sfVersion = 1.1;
bool	XSQ::sbInitialized = false;
int		XSQ::siHowMany = 0;

unique_ptr<TTree> XSQ::tree = nullptr;

double XSQ::sdXsq_n[4]			= {0,0,0,0};
double XSQ::sdPeakheight_n[4]	= {0,0,0,0};
double XSQ::sdBaseline_n[4]		= {0,0,0,0};
double XSQ::sdOffset_n[4]		= {0,0,0,0};
double XSQ::sdPeakErr_n[4]		= {0,0,0,0};
double XSQ::sdBaseErr_n[4]		= {0,0,0,0};
double XSQ::sdOffsetErr_n[4]	= {0,0,0,0};
double XSQ::sdProb_n[4]			= {0,0,0,0};

double XSQ::sdXsq_y[4]			= {0,0,0,0};
double XSQ::sdPeakheight_y[4]	= {0,0,0,0};
double XSQ::sdBaseline_y[4]		= {0,0,0,0};
double XSQ::sdOffset_y[4]		= {0,0,0,0};
double XSQ::sdPeakErr_y[4]		= {0,0,0,0};
double XSQ::sdBaseErr_y[4]		= {0,0,0,0};
double XSQ::sdOffsetErr_y[4]	= {0,0,0,0};
double XSQ::sdProb_y[4]			= {0,0,0,0};

int XSQ::siFitStatus_n[4]		= {0,0,0,0};
int XSQ::siFitStatus_y[4]		= {0,0,0,0};

XSQ::XSQ() {
	if (g_verbose) cout << "XSQ c'tor\n";
	XSQ::siHowMany++;
}

XSQ::XSQ(int ch, int length, shared_ptr<Digitizer> digitizer) : ciNPar(4), Method(ch, length, digitizer) {
	if (g_verbose) cout << "XSQ " << id << " c'tor\n";
	XSQ::siHowMany++;
	unique_ptr<TFile> std_file = nullptr;
	TVectorT<double>* pWave = nullptr;
	int i(0), p(0);
	if ((id > 3) || (id < 0)) iFailed |= method_error;
	switch (digitizer->ID()) { // setting the length and stuff for the standard events
		case dt5751 : 
			if (digitizer->Special() == 0) {
				iStdLength	= 225;
				iStdTrig	= 32;
				iStdPeakX	= 35;
			} else {
				iStdLength	= 450;
				iStdTrig	= 64;
				iStdPeakX	= 71;
			} break;
		case dt5751des :
			iStdLength		= 899;
			iStdTrig		= 128;
			iStdPeakX		= 142;
			break;
		case dt5730 :
			iStdLength		= 225;
			iStdTrig		= 32;
			iStdPeakX		= 35;
			break;
		case v1724 :
		case invalid_dig :
		default :
			iFailed |= dig_error;
			return;
	}
	
	try {std_file.reset(new TFile((path + "/config/standard_events.root").c_str(), "READ"));}
	catch (bad_alloc& ba) {iFailed |= alloc_error; return;}
	if (!std_file->IsOpen()) {iFailed |= file_error; return;}
	for (p = 0; p < P; p++) {
		try {dStdWave[p] = unique_ptr<double[]>(new double[iStdLength]);}
		catch (bad_alloc& ba) {iFailed |= alloc_error; return;}
		pWave = (TVectorT<double>*)std_file->Get((p ? "gamma_wave" : "neutron_wave"));
		if (pWave == nullptr) {iFailed |= root_error; return;}
		dStdBase[p] = 0;
		dStdPeak[p] = 1000;
		switch(iStdLength) {
			case 225 : // 500 MSa/s
				for (i = 0; i < iStdLength; i++) { // averages
					dStdWave[p][i] = ((*pWave)[2*i] + (*pWave)[2*i+1])/2.;
					if (i < digitizer->Baselength()) dStdBase[p] += dStdWave[p][i];
					dStdPeak[p] = min(dStdPeak[p], dStdWave[p][i]);
				} break;
			case 899 : // 2 GSa/s
				for (i = 0; i < iStdLength; i++) { // interpolates
					dStdWave[p][i] = (i%2) ? ((*pWave)[(i+1)/2] + (*pWave)[(i-1)/2])/2. : (*pWave)[i/2]; // i%2==1 so i/2 = (i-1)/2
					if (i < digitizer->Baselength()) dStdBase[p] += dStdWave[p][i];
					dStdPeak[p] = min(dStdPeak[p], dStdWave[p][i]);
				} break;
			case 450 : // 1 GSa/s
				for (i = 0; i < iStdLength; i++) {
					dStdWave[p][i] = (*pWave)[i];
					if (i < digitizer->Baselength()) dStdBase[p] += dStdWave[p][i];
					dStdPeak[p] = min(dStdPeak[p], dStdWave[p][i]);
				} break;
			default : iFailed |= method_error;
			return;
		}
		dStdBase[p] /= digitizer->Baselength();
		dStdNorm[p] = 1./(dStdBase[p]-dStdPeak[p]);
	} // p
	std_file->Close();
	std_file = nullptr;
	pWave = nullptr;
	try {
		dInputWave		= unique_ptr<double[]>(new double[iEventlength]);
		dX				= unique_ptr<double[]>(new double[iEventlength]);
		
		dPars			= unique_ptr<double[]>(new double[ciNPar]);
		
		fit				= unique_ptr<TF1>(new TF1("fit", this, &XSQ::fitter, 0, min(iEventlength, iStdLength), ciNPar));
	} catch (bad_alloc& ba) {iFailed |= alloc_error; return;}
	if (fit->IsZombie()) {iFailed |= root_error; return;}
	for (i = 0; i < iEventlength; i++) dX[i] = i;
	
	fit->SetParNames("Peakheight_scale","Baseline_offset","Trigger_offset","Particle");

	graph = nullptr;
}

XSQ::~XSQ() {
	if (g_verbose) cout << "XSQ " << id << " d'tor\n";
	XSQ::siHowMany--;
	for (auto p = 0; p < P; p++) {
		dStdWave[p].reset();
	}
	dPars.reset();
	fit.reset();
	graph.reset();
	dInputWave.reset();
	dX.reset();
}

void XSQ::root_init(TTree* tree_in) {
	if (!XSQ::sbInitialized) {
		XSQ::tree = unique_ptr<TTree>(tree_in);
		
		XSQ::tree->Branch("Chisquare_n",	XSQ::sdXsq_n,			"xsqn[4]/D");
		XSQ::tree->Branch("Peakscale_n",	XSQ::sdPeakheight_n,	"pkscalen[4]/D");
		XSQ::tree->Branch("Base_shift_n",	XSQ::sdBaseline_n,		"baseshiftn[4]/D");
		XSQ::tree->Branch("Offset_n",		XSQ::sdOffset_n,		"offsetn[4]/D");
		XSQ::tree->Branch("Prob_n",			XSQ::sdProb_n,			"probn[4]/D");
		
		XSQ::tree->Branch("Chisquare_y",	XSQ::sdXsq_y,			"xsqy[4]/D");
		XSQ::tree->Branch("Peakscale_y",	XSQ::sdPeakheight_y,	"pkscaley[4]/D");
		XSQ::tree->Branch("Base_shift_y",	XSQ::sdBaseline_y,		"baseshifty[4]/D");
		XSQ::tree->Branch("Offset_y",		XSQ::sdOffset_y,		"offsety[4]/D");
		XSQ::tree->Branch("Prob_y",			XSQ::sdProb_y,			"proby[4]/D");
		
		XSQ::tree->Branch("Peak_err_n",		XSQ::sdPeakErr_n,		"pkerrn[4]/D");
		XSQ::tree->Branch("Base_err_n",		XSQ::sdBaseErr_n,		"berrn[4]/D");
		XSQ::tree->Branch("Offset_err_n",	XSQ::sdOffsetErr_n,		"offerrn[4]/D");
		XSQ::tree->Branch("Peak_err_y",		XSQ::sdPeakErr_y,		"pkerry[4]/D");
		XSQ::tree->Branch("Base_err_y",		XSQ::sdBaseErr_y,		"berry[4]/D");
		XSQ::tree->Branch("Offset_err_y",	XSQ::sdOffsetErr_y,		"offerry[4]/D");
		
		XSQ::tree->Branch("Fit_status_n",	XSQ::siFitStatus_n,		"fitstatusn[4]/i");
		XSQ::tree->Branch("Fit_status_y",	XSQ::siFitStatus_y,		"fitstatusy[4]/i");
		
		XSQ::sbInitialized = true;
	}
}

double XSQ::fitter(double* x, double* par) {
	int iSample = (x[0] - par[2]), p = par[3];
	double dVal = dStdBase[p] + par[1];
	if ((iSample > -1) && (iSample < iStdLength)) dVal += par[0]*(dStdWave[p][iSample]-dStdBase[p]);
	if (dVal < 0) dVal = 0; // saturated events
	return dVal;
}

void XSQ::SetParameters(void* val, int which, shared_ptr<Digitizer> digitizer) {
	switch (which) {
		case n : fGain[n] = *((float*)val); break;
		case y : fGain[y] = *((float*)val); break;
		default: break;
	}
}

void XSQ::Analyze() {
	for (auto i = 0; i < iEventlength; i++) dInputWave[i] = event->Trace(i);
	try {graph.reset(new TGraph(iEventlength, dX.get(), dInputWave.get()));}
	catch (bad_alloc& ba) { // error codes don't work here
		XSQ::sdXsq_n[id]		= -1;
		XSQ::sdPeakheight_n[id]	= -1;
		XSQ::sdBaseline_n[id]	= -1;
		XSQ::sdOffset_n[id]		= -1;
		XSQ::sdProb_n[id]		= -1;
		XSQ::sdPeakErr_n[id]	= -1;
		XSQ::sdBaseErr_n[id]	= -1;
		XSQ::sdOffsetErr_n[id]	= -1;
		
		XSQ::sdXsq_y[id]		= -1;
		XSQ::sdPeakheight_y[id]	= -1;
		XSQ::sdBaseline_y[id]	= -1;
		XSQ::sdOffset_y[id]		= -1;
		XSQ::sdProb_y[id]		= -1;
		XSQ::sdPeakErr_y[id]	= -1;
		XSQ::sdBaseErr_y[id]	= -1;
		XSQ::sdOffsetErr_y[id]	= -1;
		
		XSQ::siFitStatus_n[id]	= -1;
		XSQ::siFitStatus_y[id]	= -1;
		
		return;
	}
	
	dPars[0] = (event->Baseline() - event->Peak_y())*dStdNorm[n];
	dPars[1] = event->Baseline() - dStdBase[n];
	dPars[2] = event->Trigger() - iStdTrig;
	dPars[3] = n;
	fit->SetParameters(dPars.get());
	fit->FixParameter(3, n);
	
	XSQ::siFitStatus_n[id]	= graph->Fit(fit.get(), "Q N R"); // quiet, no-plot, specified range

	XSQ::sdXsq_n[id]		= fit->GetChisquare();
	XSQ::sdPeakheight_n[id]	= fit->GetParameter(0)*fGain[n]; // detectors have different gains
	XSQ::sdBaseline_n[id]	= fit->GetParameter(1);
	XSQ::sdOffset_n[id]		= fit->GetParameter(2);
	XSQ::sdProb_n[id]		= fit->GetProb();

	XSQ::sdPeakErr_n[id]	= fit->GetParError(0)*fGain[n];
	XSQ::sdBaseErr_n[id]	= fit->GetParError(1);
	XSQ::sdOffsetErr_n[id]	= fit->GetParError(2);

	
	dPars[0] = (event->Baseline() - event->Peak_y())*dStdNorm[y];
	dPars[1] = event->Baseline() - dStdBase[y];
	dPars[2] = event->Trigger() - iStdTrig;
	dPars[3] = y;
	fit->SetParameters(dPars.get());
	fit->FixParameter(3, y);
	
	XSQ::siFitStatus_y[id]	= graph->Fit(fit.get(), "Q N R");

	XSQ::sdXsq_y[id]		= fit->GetChisquare();
	XSQ::sdPeakheight_y[id]	= fit->GetParameter(0)*fGain[y];
	XSQ::sdBaseline_y[id]	= fit->GetParameter(1);
	XSQ::sdOffset_y[id]		= fit->GetParameter(2);
	XSQ::sdProb_y[id]		= fit->GetProb();

	XSQ::sdPeakErr_y[id]	= fit->GetParError(0)*fGain[y];
	XSQ::sdBaseErr_y[id]	= fit->GetParError(1);
	XSQ::sdOffsetErr_y[id]	= fit->GetParError(2);
}
