#include "XSQ.h"
#include <cstdlib>
#include <algorithm>
#include "TVectorT.h"

float	XSQ::sfVersion = 1.5;
bool	XSQ::sbInitialized = false;
int		XSQ::siHowMany = 0;

unique_ptr<TTree> XSQ::tree = nullptr;

double XSQ::sdXsq_n[4]			= {0,0,0,0};
double XSQ::sdPeakheight_n[4]	= {0,0,0,0};
double XSQ::sdBaseline_n[4]		= {0,0,0,0};
double XSQ::sdOffset_n[4]		= {0,0,0,0};

double XSQ::sdXsq_y[4]			= {0,0,0,0};
double XSQ::sdPeakheight_y[4]	= {0,0,0,0};
double XSQ::sdBaseline_y[4]		= {0,0,0,0};
double XSQ::sdOffset_y[4]		= {0,0,0,0};

double XSQ::sdPeak_err_n[4]		= {0,0,0,0};
double XSQ::sdBase_err_n[4]		= {0,0,0,0};
double XSQ::sdOff_err_n[4]		= {0,0,0,0};
double XSQ::sdPeak_err_y[4]		= {0,0,0,0};
double XSQ::sdBase_err_y[4]		= {0,0,0,0};
double XSQ::sdOff_err_y[4]		= {0,0,0,0};

short XSQ::ssIterations[2][4]		= {{0,0,0,0},{0,0,0,0}};

XSQ::XSQ() : ciVariant(-1) {
	if (g_verbose) cout << "XSQ c'tor\n";
	XSQ::siHowMany++;
}

XSQ::XSQ(int ch, int length, shared_ptr<Digitizer> digitizer, int variant) : Method(ch, length, digitizer), ciVariant(variant) {
	if (g_verbose) cout << "XSQ " << id << " c'tor\n";
	XSQ::siHowMany++;
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
		fit = unique_ptr<TF1>(new TF1("fit",this,&XSQ::TF1_fit_func,0,min(iStdLength, iEventlength),TF1_pars));
	} catch (bad_alloc& ba) {
		iFailed |= (1 << alloc_error);
		return;
	}
	fit->SetParNames("Peakscale","Baseline","Offset","particle");
	for (i = 0; i < iEventlength; i++) dX[i] = i;
	if ((ciVariant != VAR_TF1) && (ciVariant != VAR_NEW)) {iFailed |= (1 << method_error); return;}
	iIterations = 5;
	dConvergence = 0.001; // 0.1% change
}

XSQ::~XSQ() {
	if (g_verbose) cout << "XSQ " << id << " d'tor\n";
	XSQ::siHowMany--;
	dInputWave.reset();
	dX.reset();
	fit.reset();
	graph.reset();
	for (auto p = 0; p < P; p++) {
		dStdWave[p].reset();
	}
}

void XSQ::root_init(TTree* tree_in) {
	if (!XSQ::sbInitialized) {
		XSQ::tree = unique_ptr<TTree>(tree_in);
		
		XSQ::tree->Branch("Chisquare_n",	XSQ::sdXsq_n,			"xsqn[4]/D");
		XSQ::tree->Branch("Peakscale_n",	XSQ::sdPeakheight_n,	"pkscalen[4]/D");
		XSQ::tree->Branch("Base_shift_n",	XSQ::sdBaseline_n,		"baseshiftn[4]/D");
		XSQ::tree->Branch("Offset_n",		XSQ::sdOffset_n,		"offsetn[4]/D");
		
		XSQ::tree->Branch("Chisquare_y",	XSQ::sdXsq_y,			"xsqy[4]/D");
		XSQ::tree->Branch("Peakscale_y",	XSQ::sdPeakheight_y,	"pkscaley[4]/D");
		XSQ::tree->Branch("Base_shift_y",	XSQ::sdBaseline_y,		"baseshifty[4]/D");
		XSQ::tree->Branch("Offset_y",		XSQ::sdOffset_y,		"offsety[4]/D");
		
		XSQ::tree->Branch("Peak_err_n",		XSQ::sdPeak_err_n,		"pkerrn[4]/D");
		XSQ::tree->Branch("Base_err_n",		XSQ::sdBase_err_n,		"baerrn[4]/D");
		XSQ::tree->Branch("Off_err_n",		XSQ::sdOff_err_n,		"oferrn[4]/D");
		XSQ::tree->Branch("Peak_err_y",		XSQ::sdPeak_err_y,		"pkerry[4]/D");
		XSQ::tree->Branch("Base_err_y",		XSQ::sdBase_err_y,		"baerry[4]/D");
		XSQ::tree->Branch("Off_err_y",		XSQ::sdOff_err_y,		"oferry[4]/D");
		
		XSQ::tree->Branch("Iterations",		XSQ::ssIterations,		"iters[2][4]/S");

		XSQ::sbInitialized = true;
	}
}

void XSQ::SetDefaultParameters() { // for convenience
	XSQ::sdXsq_n[id]		= -1;
	XSQ::sdPeakheight_n[id]	= (event->Baseline() - event->Peak_y())*dStdNorm[n]*fGain[n];
	XSQ::sdBaseline_n[id]	= event->Baseline();
	XSQ::sdOffset_n[id]		= event->Trigger() - iStdTrig;

	XSQ::sdXsq_y[id]		= -1;
	XSQ::sdPeakheight_y[id]	= (event->Baseline() - event->Peak_y())*dStdNorm[y]*fGain[y];
	XSQ::sdBaseline_y[id]	= event->Baseline();
	XSQ::sdOffset_y[id]		= event->Trigger() - iStdTrig;

	XSQ::sdPeak_err_n[id]	= -1;
	XSQ::sdPeak_err_y[id]	= -1;
	XSQ::sdBase_err_n[id]	= -1;
	XSQ::sdBase_err_y[id]	= -1;
	XSQ::sdOff_err_n[id]	= -1;
	XSQ::sdOff_err_y[id]	= -1;
	
	XSQ::ssIterations[0][id] = -1;
	XSQ::ssIterations[1][id] = -1;

	return;
}

void XSQ::SetParameters(void* val, int which, shared_ptr<Digitizer> digitizer) {
	switch (which) {
		case n : fGain[n] = *((float*)val); break;
		case y : fGain[y] = *((float*)val); break;
		default: break;
	}
}

double XSQ::TF1_fit_func(double* x, double* par) {
	int iSample(x[0]-par[2]), p(par[3]);
	double dVal(par[1]);
	if ((iSample > -1) && (iSample < iStdLength)) dVal += par[0]*dStdWave[p][iSample];
	if (dVal < 0) dVal = 0; // saturated events
	return dVal;
}

double XSQ::FindChisquare(int p, double dPeak, double dBase, int iOff) {
	double dVal(0), dTrace(0), dDiff(0), dChiSquare(0);
	int iSample(0);
	for (int i = 0; i < iEventlength; i++) {
		dTrace = event->Trace(i);
		iSample = i - iOff;
		if (dTrace == 0) continue;
		dVal = dBase;
		if ((iSample > -1) && (iSample < iStdLength)) dVal += dPeak*dStdWave[p][iSample];
		dDiff = dTrace-dVal;
		dChiSquare += (dDiff*dDiff)/(dTrace*dTrace);
	}
	return dChiSquare;
}

void XSQ::Analyze() {
	int i(0), iter(iIterations);
	double dEpsilonP(0.0005), dEpsilonB(0.0001), dPeak(0), dBase(0); // looking in a small region
	double dDet(0), a(0), b(0), c(0), d(0), e(0), f(0), dDiff(0);
	double d000(0), dp00(0), d0p0(0), d00p(0), dpp0(0), dp0p(0), d0pp(0), dm00(0), d0m0(0), d00m(0);
	double dChiSquareLast(0), dChiSquare(0), dEpsilonO(1), iOff(0);
/*	if (event->Baseline() - event->Peak_y() < iPeakCut) { // cut noise to save processing time
		SetDefaultParameters();
		return;
	} */
	if (ciVariant) { // Newton's method, could use some optimization
		for (int p = 0; p < P; p++) {
			ssIterations[p][id] = 0;
			iter = iIterations; // in case extra time was used
			dPeak = (event->Baseline() - event->Peak_y())*dStdNorm[p];
			dBase = event->Baseline();
			iOff = event->Trigger() - iStdTrig;
			dChiSquare = FindChisquare(p, dPeak, dBase, iOff);
			for (i = 0; i < iter; i++) { // x_n+1 = x_n - (main step size)(inverse Hessian matrix)(gradient of function at x_n)
				d000 = dChiSquareLast = dChiSquare;
				dp00 = FindChisquare(p, dPeak + dEpsilonP, dBase, iOff);
				d0p0 = FindChisquare(p, dPeak, dBase + dEpsilonB, iOff);
				d00p = FindChisquare(p, dPeak, dBase, iOff + dEpsilonO);
				dpp0 = FindChisquare(p, dPeak + dEpsilonP, dBase + dEpsilonB, iOff);
				dp0p = FindChisquare(p, dPeak + dEpsilonP, dBase, iOff + dEpsilonO); 
				d0pp = FindChisquare(p, dPeak, dBase + dEpsilonB, iOff + dEpsilonO);
				dm00 = FindChisquare(p, dPeak - dEpsilonP, dBase, iOff);
				d0m0 = FindChisquare(p, dPeak, dBase - dEpsilonB, iOff);
				d00m = FindChisquare(p, dPeak, dBase, iOff - dEpsilonO);
				
				a = (dp00 - 2*d000 + dm00)/(dEpsilonP*dEpsilonP);
				b = (dpp0 - dp00 - d0p0 + d000)/(dEpsilonB*dEpsilonP);
				c = (dp0p - dp00 - d00p + d000)/(dEpsilonP); // iEpsilonO == 1
				d = (d0p0 - 2*d000 + d0m0)/(dEpsilonB*dEpsilonB);
				e = (d0pp - d0p0 - d00p + d000)/(dEpsilonB);
				f = (d00p - 2*d000 + d00m);
				
				dDet = a*(d*f-e*e)-b*(b*f-c*e)+c*(b*e-c*d);
				if (dDet == 0) { // not sure what to do about this, really
					SetDefaultParameters();
					return;
				}
				dDet = 1./dDet; // save a few clock cycles
				
				dHessianInv[0][0] = (d*f-e*e)*dDet;
				dHessianInv[0][1] = dHessianInv[1][0] = (c*e-b*f)*dDet; // | a b c |     1  | x u v |  x = df-ee, u = ce-bf
				dHessianInv[0][2] = dHessianInv[2][0] = (b*e-c*d)*dDet; // | b d e | -> --- | u y w |, v = be-cd, w = bc-ae
				dHessianInv[1][1] = (a*f-c*c)*dDet;						// | c e f |    det | v w z |  y = af-cc, z = ad-bb
				dHessianInv[1][2] = dHessianInv[2][1] = (b*c-a*e)*dDet;
				dHessianInv[2][2] = (a*d-b*b)*dDet;
				
				dGradient[0] = (dp00-dm00)/(2*dEpsilonP);
				dGradient[1] = (d0p0-d0m0)/(2*dEpsilonB);
				dGradient[2] = (d00p-d00m)/2;
				
				dStep[0] = dGradient[0]*dHessianInv[0][0] + dGradient[1]*dHessianInv[0][1] + dGradient[2]*dHessianInv[0][2];
				dStep[1] = dGradient[0]*dHessianInv[1][0] + dGradient[1]*dHessianInv[1][1] + dGradient[2]*dHessianInv[1][2];
				dStep[2] = dGradient[0]*dHessianInv[2][0] + dGradient[1]*dHessianInv[2][1] + dGradient[2]*dHessianInv[2][2];
				
				dPeak -= dStep[0];
				dBase -= dStep[1];
				iOff -= dStep[2];
				
				dChiSquare = FindChisquare(p, dPeak, dBase, iOff);
				dDiff = dChiSquare - dChiSquareLast;
				if (dDiff < 0) dDiff *= -1;
				if (dDiff/dChiSquareLast < dConvergence) { // convergence check, 0.1%
					if (p) {
						XSQ::sdXsq_y[id]		= dChiSquare;
						XSQ::sdPeakheight_y[id]	= dPeak*fGain[y]/iResolutionScale;
						XSQ::sdBaseline_y[id]	= dBase;
						XSQ::sdOffset_y[id]		= iOff;
						
						XSQ::sdPeak_err_y[id]	= dStep[0] > 0 ? dStep[0] : -dStep[0];
						XSQ::sdBase_err_y[id]	= dStep[1] > 0 ? dStep[1] : -dStep[1];
						XSQ::sdOff_err_y[id]	= dStep[2] > 0 ? dStep[2] : -dStep[2];
						
					} else {
						XSQ::sdXsq_n[id]		= dChiSquare;
						XSQ::sdPeakheight_n[id]	= dPeak*fGain[n]/iResolutionScale;
						XSQ::sdBaseline_n[id]	= dBase;
						XSQ::sdOffset_n[id]		= iOff;
						
						XSQ::sdPeak_err_n[id]	= dStep[0] > 0 ? dStep[0] : -dStep[0];
						XSQ::sdBase_err_n[id]	= dStep[1] > 0 ? dStep[1] : -dStep[1];
						XSQ::sdOff_err_n[id]	= dStep[2] > 0 ? dStep[2] : -dStep[2];
					}
					XSQ::ssIterations[p][id] = i;
					break; // break from iterations
				} else if (i == iIterations-1) { // more time to converge, only once
					iter += iIterations; 
				}
			} // iterations loop
		} // p loop
	} else { // TF1
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
		
		XSQ::sdXsq_n[id]		= fit->GetChisquare();
		XSQ::sdPeakheight_n[id]	= fit->GetParameter(0)*fGain[n]/iResolutionScale;
		XSQ::sdBaseline_n[id]	= fit->GetParameter(1);
		XSQ::sdOffset_n[id]		= fit->GetParameter(2);
		
		XSQ::sdPeak_err_n[id]	= fit->GetParError(0)*fGain[n]/iResolutionScale;
		XSQ::sdBase_err_n[id]	= fit->GetParError(1);
		XSQ::sdOff_err_n[id]	= fit->GetParError(2);
		
		
		fit->SetParameter(0, (event->Baseline() - event->Peak_y())*dStdNorm[y]);
		fit->SetParameter(1, event->Baseline());
		fit->SetParameter(2, event->Trigger() - iStdTrig);
		fit->SetParameter(3, y);
		fit->FixParameter(3, y);
		graph->Fit(fit.get(), "Q N R");
		
		XSQ::sdXsq_y[id]		= fit->GetChisquare();
		XSQ::sdPeakheight_y[id]	= fit->GetParameter(0)*fGain[y]/iResolutionScale;
		XSQ::sdBaseline_y[id]	= fit->GetParameter(1);
		XSQ::sdOffset_y[id]		= fit->GetParameter(2);
		
		XSQ::sdPeak_err_y[id]	= fit->GetParError(0)*fGain[y]/iResolutionScale;
		XSQ::sdBase_err_y[id]	= fit->GetParError(1);
		XSQ::sdOff_err_y[id]	= fit->GetParError(2);
	}
}
