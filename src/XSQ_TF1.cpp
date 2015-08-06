#include "XSQ_TF1.h"
#include "TVectorT.h"

float	XSQ_TF1::sfVersion = 1.56;
bool	XSQ_TF1::sbInitialized = false;
int		XSQ_TF1::siHowMany = 0;

unique_ptr<TTree> XSQ_TF1::tree = nullptr;

double XSQ_TF1::sdXsq_n[4]			= {0,0,0,0};
double XSQ_TF1::sdPeakheight_n[4]	= {0,0,0,0};
double XSQ_TF1::sdBaseline_n[4]		= {0,0,0,0};
double XSQ_TF1::sdOffset_n[4]		= {0,0,0,0};

double XSQ_TF1::sdXsq_y[4]			= {0,0,0,0};
double XSQ_TF1::sdPeakheight_y[4]	= {0,0,0,0};
double XSQ_TF1::sdBaseline_y[4]		= {0,0,0,0};
double XSQ_TF1::sdOffset_y[4]		= {0,0,0,0};

double XSQ_TF1::sdBaseline_a[4]		= {0,0,0,0};
double XSQ_TF1::sdSigma_a[4]		= {0,0,0,0};

double XSQ_TF1::sdPeak_err_n[4]		= {0,0,0,0};
double XSQ_TF1::sdBase_err_n[4]		= {0,0,0,0};
double XSQ_TF1::sdOff_err_n[4]		= {0,0,0,0};
double XSQ_TF1::sdPeak_err_y[4]		= {0,0,0,0};
double XSQ_TF1::sdBase_err_y[4]		= {0,0,0,0};
double XSQ_TF1::sdOff_err_y[4]		= {0,0,0,0};

XSQ_TF1::XSQ_TF1() {
	if (g_verbose) cout << "XSQ_TF1 c'tor\n";
	XSQ_TF1::siHowMany++;
}

XSQ_TF1::XSQ_TF1(int ch, int length, shared_ptr<Digitizer> digitizer) : Method(ch, length, digitizer) {
	if (g_verbose) cout << "XSQ_TF1 " << id << " c'tor\n";
	XSQ_TF1::siHowMany++;
	unique_ptr<TFile> std_file = nullptr;
	TVectorT<double>* pWave = nullptr;
	int i(0), p(0);
	if ((id > 3) || (id < 0)) iFailed |= (1 << method_error);
	switch (digitizer->ID()) { // setting the length and stuff for the standard events
		case dt5751 :
			iPeakCut = 4;
			if (digitizer->Special() == 0) {
				iStdLength	= 225;
				iStdTrig	= 32;
				iResolutionScale = 1;
			} else {
				iStdLength	= 450;
				iStdTrig	= 64;
				iResolutionScale = 1;
			} break;
		case dt5751des :
			iPeakCut		= 4;
			iStdLength		= 899;
			iStdTrig		= 128;
			iResolutionScale = 1;
			break;
		case dt5730 :
			iPeakCut		= 32;
			iStdLength		= 225;
			iStdTrig		= 32;
			iResolutionScale = 1 << 4;
			break;
		case v1724 :
		case invalid_dig :
		default :
			iFailed |= (1 << dig_error);
			iFailed |= (1 << method_error);
			return;
	}
	
	try {std_file.reset(new TFile((sWorkingDir + "/Data_Processor/config/standard_events.root").c_str(), "READ"));}
	catch (bad_alloc& ba) {iFailed |= alloc_error; return;}
	if (!std_file->IsOpen()) {iFailed |= file_error; return;}
	for (p = 0; p < P; p++) {
		try {dStdWave[p] = unique_ptr<double[]>(new double[iStdLength]);}
		catch (bad_alloc& ba) {iFailed |= alloc_error; return;}
		pWave = (TVectorT<double>*)std_file->Get((p ? "gamma_wave_inv" : "neutron_wave_inv"));
		if (pWave == nullptr) {iFailed |= root_error; return;}
		dStdPeak[p] = 1000;
		switch(iStdLength) {
			case 225 : // 500 MSa/s
				for (i = 0; i < iStdLength; i++) { // averages
					dStdWave[p][i] = ((*pWave)[2*i] + (*pWave)[2*i+1])/2.;
					dStdPeak[p] = min(dStdPeak[p], dStdWave[p][i]);
				} break;
			case 899 : // 2 GSa/s
				for (i = 0; i < iStdLength; i++) { // interpolates
					dStdWave[p][i] = (i%2) ? ((*pWave)[(i+1)/2] + (*pWave)[(i-1)/2])/2. : (*pWave)[i/2]; // i%2==1 so i/2 = (i-1)/2
					dStdPeak[p] = min(dStdPeak[p], dStdWave[p][i]);
				} break;
			case 450 : // 1 GSa/s
				for (i = 0; i < iStdLength; i++) {
					dStdWave[p][i] = (*pWave)[i];
					dStdPeak[p] = min(dStdPeak[p], dStdWave[p][i]);
				} break;
			default : iFailed |= (1 << method_error);
			return;
		}
		dStdNorm[p] = -1./dStdPeak[p]*(iResolutionScale); // waves generated on 10-bit digitizer so must be scaled appropriately for other resolutions
	} // p
	std_file->Close();
	std_file = nullptr;
	pWave = nullptr;
	graph = nullptr;
	
	try {
		dInputWave = unique_ptr<double[]>(new double[iEventlength]);
		dX = unique_ptr<double[]>(new double[iEventlength]);
		fit = unique_ptr<TF1>(new TF1("fit",this,&XSQ_TF1::TF1_fit_func,0,iEventlength,ciNPar));
	} catch (bad_alloc& ba) {
		iFailed |= (1 << alloc_error);
		return;
	}
	fit->SetParNames("Peakscale","Baseline","Offset","particle");
	for (i = 0; i < iEventlength; i++) dX[i] = i;
}

XSQ_TF1::~XSQ_TF1() {
	if (g_verbose) cout << "XSQ_TF1 " << id << " d'tor\n";
	XSQ_TF1::siHowMany--;
	dInputWave.reset();
	dX.reset();
	fit.reset();
	graph.reset();
	for (auto p = 0; p < P; p++) {
		dStdWave[p].reset();
	}
}

void XSQ_TF1::root_init(TTree* tree_in) {
	if (!XSQ_TF1::sbInitialized) {
		XSQ_TF1::tree = unique_ptr<TTree>(tree_in);
		
		XSQ_TF1::tree->Branch("Chisquare_n",	XSQ_TF1::sdXsq_n,			"xsqn[4]/D");
		XSQ_TF1::tree->Branch("Peakscale_n",	XSQ_TF1::sdPeakheight_n,	"pkscalen[4]/D");
		XSQ_TF1::tree->Branch("Base_shift_n",	XSQ_TF1::sdBaseline_n,		"baseshiftn[4]/D");
		XSQ_TF1::tree->Branch("Offset_n",		XSQ_TF1::sdOffset_n,		"offsetn[4]/D");
		
		XSQ_TF1::tree->Branch("Chisquare_y",	XSQ_TF1::sdXsq_y,			"xsqy[4]/D");
		XSQ_TF1::tree->Branch("Peakscale_y",	XSQ_TF1::sdPeakheight_y,	"pkscaley[4]/D");
		XSQ_TF1::tree->Branch("Base_shift_y",	XSQ_TF1::sdBaseline_y,		"baseshifty[4]/D");
		XSQ_TF1::tree->Branch("Offset_y",		XSQ_TF1::sdOffset_y,		"offsety[4]/D");
		
		XSQ_TF1::tree->Branch("Base",			XSQ_TF1::sdBaseline_a,		"baseline[4]/D");
		XSQ_TF1::tree->Branch("Sigma",			XSQ_TF1::sdSigma_a,			"sigma[4]/D");
		
		XSQ_TF1::tree->Branch("Peak_err_n",		XSQ_TF1::sdPeak_err_n,		"pkerrn[4]/D");
		XSQ_TF1::tree->Branch("Base_err_n",		XSQ_TF1::sdBase_err_n,		"baerrn[4]/D");
		XSQ_TF1::tree->Branch("Off_err_n",		XSQ_TF1::sdOff_err_n,		"oferrn[4]/D");
		XSQ_TF1::tree->Branch("Peak_err_y",		XSQ_TF1::sdPeak_err_y,		"pkerry[4]/D");
		XSQ_TF1::tree->Branch("Base_err_y",		XSQ_TF1::sdBase_err_y,		"baerry[4]/D");
		XSQ_TF1::tree->Branch("Off_err_y",		XSQ_TF1::sdOff_err_y,		"oferry[4]/D");

		XSQ_TF1::sbInitialized = true;
	}
}

void XSQ_TF1::SetDefaultParameters() { // for convenience
	XSQ_TF1::sdXsq_n[id]		= -1;
	XSQ_TF1::sdPeakheight_n[id]	= (event->Baseline() - event->Peak_y())*dStdNorm[n]*fGain[n];
	XSQ_TF1::sdBaseline_n[id]	= event->Baseline();
	XSQ_TF1::sdOffset_n[id]		= event->Trigger() - iStdTrig;

	XSQ_TF1::sdXsq_y[id]		= -1;
	XSQ_TF1::sdPeakheight_y[id]	= (event->Baseline() - event->Peak_y())*dStdNorm[y]*fGain[y];
	XSQ_TF1::sdBaseline_y[id]	= event->Baseline();
	XSQ_TF1::sdOffset_y[id]		= event->Trigger() - iStdTrig;
	
	XSQ_TF1::sdBaseline_a[id]	= event->Baseline();
	XSQ_TF1::sdSigma_a[id]		= event->BaseSigma();

	XSQ_TF1::sdPeak_err_n[id]	= -1;
	XSQ_TF1::sdPeak_err_y[id]	= -1;
	XSQ_TF1::sdBase_err_n[id]	= -1;
	XSQ_TF1::sdBase_err_y[id]	= -1;
	XSQ_TF1::sdOff_err_n[id]	= -1;
	XSQ_TF1::sdOff_err_y[id]	= -1;

	return;
}

void XSQ_TF1::SetParameters(void* val, int which, shared_ptr<Digitizer> digitizer) {
	switch (which) {
		case n : fGain[n] = *((float*)val); break;
		case y : fGain[y] = *((float*)val); break;
		default: break;
	}
}

double XSQ_TF1::TF1_fit_func(double* x, double* par) {
	int iSample(x[0]-par[2]), p(par[3]);
	double dVal(par[1]);
	if ((iSample > -1) && (iSample < iStdLength)) dVal += par[0]*dStdWave[p][iSample];
	if (dVal < 0) dVal = 0; // saturated events
	return dVal;
}

void XSQ_TF1::Analyze() {
	int i(0);
	double dBaseline(event->Baseline()), dPeakY(event->Peak_y());
	SetDefaultParameters();
	if (dBaseline - dPeakY < iPeakCut) { // cut noise to save processing time
		return;
	}
	for (i = 0; i < iEventlength; i++) dInputWave[i] = event->Trace(i);
	try {graph.reset(new TGraph(iEventlength, dX.get(), dInputWave.get()));}
	catch (bad_alloc& ba) {
		SetDefaultParameters();
		return;
	}
	fit->SetParameter(0, (event->Baseline() - event->Peak_y())*dStdNorm[n]);
	fit->SetParameter(1, event->Baseline());
	fit->SetParameter(2, event->Trigger() - iStdTrig);
	fit->SetParameter(3, n);
	fit->FixParameter(3, n);
	graph->Fit(fit.get(), "Q N R"); // quiet, no-plot, range
	
	XSQ_TF1::sdXsq_n[id]		= fit->GetChisquare();
	XSQ_TF1::sdPeakheight_n[id]	= fit->GetParameter(0)/fGain[n]/iResolutionScale;
	XSQ_TF1::sdBaseline_n[id]	= fit->GetParameter(1);
	XSQ_TF1::sdOffset_n[id]		= fit->GetParameter(2);
	
	XSQ_TF1::sdPeak_err_n[id]	= fit->GetParError(0)/fGain[n]/iResolutionScale;
	XSQ_TF1::sdBase_err_n[id]	= fit->GetParError(1);
	XSQ_TF1::sdOff_err_n[id]	= fit->GetParError(2);
	
	
	fit->SetParameter(0, (event->Baseline() - event->Peak_y())*dStdNorm[y]);
	fit->SetParameter(1, event->Baseline());
	fit->SetParameter(2, event->Trigger() - iStdTrig);
	fit->SetParameter(3, y);
	fit->FixParameter(3, y);
	graph->Fit(fit.get(), "Q N R");
	
	XSQ_TF1::sdXsq_y[id]		= fit->GetChisquare();
	XSQ_TF1::sdPeakheight_y[id]	= fit->GetParameter(0)/fGain[y]/iResolutionScale;
	XSQ_TF1::sdBaseline_y[id]	= fit->GetParameter(1);
	XSQ_TF1::sdOffset_y[id]		= fit->GetParameter(2);
	
	XSQ_TF1::sdPeak_err_y[id]	= fit->GetParError(0)/fGain[y]/iResolutionScale;
	XSQ_TF1::sdBase_err_y[id]	= fit->GetParError(1);
	XSQ_TF1::sdOff_err_y[id]	= fit->GetParError(2);

}
