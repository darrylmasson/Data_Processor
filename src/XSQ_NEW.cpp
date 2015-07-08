#include "XSQ_NEW.h"
#include "TVectorT.h"

float	XSQ_NEW::sfVersion = 1.06;
bool	XSQ_NEW::sbInitialized = false;
int		XSQ_NEW::siHowMany = 0;

unique_ptr<TTree> XSQ_NEW::tree = nullptr;

double XSQ_NEW::sdXsq[2][4]			= {{0,0,0,0},{0,0,0,0}};
double XSQ_NEW::sdPeakheight[2][4]	= {{0,0,0,0},{0,0,0,0}};
double XSQ_NEW::sdBaseline[2][4]	= {{0,0,0,0},{0,0,0,0}};
double XSQ_NEW::sdOffset[2][4]		= {{0,0,0,0},{0,0,0,0}};

double XSQ_NEW::sdPeak_err[2][4]	= {{0,0,0,0},{0,0,0,0}};
double XSQ_NEW::sdBase_err[2][4]	= {{0,0,0,0},{0,0,0,0}};
double XSQ_NEW::sdOff_err[2][4]		= {{0,0,0,0},{0,0,0,0}};

short XSQ_NEW::ssIterations[2][4]	= {{0,0,0,0},{0,0,0,0}};
double XSQ_NEW::sdConvergence[2][4]	= {{0,0,0,0},{0,0,0,0}};

XSQ_NEW::XSQ_NEW() {
	if (g_verbose) cout << "XSQ_NEW c'tor\n";
	XSQ_NEW::siHowMany++;
}

XSQ_NEW::XSQ_NEW(int ch, int length, shared_ptr<Digitizer> digitizer) : Method(ch, length, digitizer) {
	if (g_verbose) cout << "XSQ_NEW " << id << " c'tor\n";
	XSQ_NEW::siHowMany++;
	unique_ptr<TFile> std_file = nullptr;
	TVectorT<double>* pWave = nullptr;
	int i(0), p(0);
	if ((id > 4) || (id < 0)) iFailed |= (1 << method_error);
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

	iIterations = 5;
	dConvergence = 0.001; // 0.1% change
}

XSQ_NEW::~XSQ_NEW() {
	if (g_verbose) cout << "XSQ_NEW " << id << " d'tor\n";
	XSQ_NEW::siHowMany--;
	for (auto p = 0; p < P; p++) {
		dStdWave[p].reset();
	}
}

void XSQ_NEW::root_init(TTree* tree_in) {
	if (!XSQ_NEW::sbInitialized) {
		XSQ_NEW::tree = unique_ptr<TTree>(tree_in);
		
		XSQ_NEW::tree->Branch("Chisquare",	XSQ_NEW::sdXsq,			"xsq[2][4]/D");
		XSQ_NEW::tree->Branch("Peakscale",	XSQ_NEW::sdPeakheight,	"pkscale[2][4]/D");
		XSQ_NEW::tree->Branch("Base_shift",	XSQ_NEW::sdBaseline,	"baseshift[2][4]/D");
		XSQ_NEW::tree->Branch("Offset",		XSQ_NEW::sdOffset,		"offset[2][4]/D");
	
		XSQ_NEW::tree->Branch("Peak_err",	XSQ_NEW::sdPeak_err,	"pkerr[2][4]/D");
		XSQ_NEW::tree->Branch("Base_err",	XSQ_NEW::sdBase_err,	"baerr[2][4]/D");
		XSQ_NEW::tree->Branch("Off_err",	XSQ_NEW::sdOff_err,		"oferr[2][4]/D");
		
		XSQ_NEW::tree->Branch("Iterations",	XSQ_NEW::ssIterations,	"iters[2][4]/S");
		XSQ_NEW::tree->Branch("Convergence",XSQ_NEW::sdConvergence, "conv[2][4]/D");

		XSQ_NEW::sbInitialized = true;
	}
}

void XSQ_NEW::SetDefaultParameters() { // for convenience
	for (int p = 0; p < P; p++) {
		XSQ_NEW::sdXsq[p][id]			= -1;
		XSQ_NEW::sdPeakheight[p][id]	= (event->Baseline() - event->Peak_y())*dStdNorm[n]*fGain[n];
		XSQ_NEW::sdBaseline[p][id]		= event->Baseline();
		XSQ_NEW::sdOffset[p][id]		= event->Trigger() - iStdTrig;

		XSQ_NEW::sdPeak_err[p][id]		= -1;
		XSQ_NEW::sdBase_err[p][id]		= -1;
		XSQ_NEW::sdOff_err[p][id]		= -1;
		
		XSQ_NEW::sdConvergence[p][id]	= -1;
	}
	return;
}

void XSQ_NEW::SetParameters(void* val, int which, shared_ptr<Digitizer> digitizer) {
	switch (which) {
		case n : fGain[n] = *((float*)val); break;
		case y : fGain[y] = *((float*)val); break;
		default: break;
	}
}

double XSQ_NEW::FindChisquare(int p, double dPeak, double dBase, int iOff) {
	double dVal(0), dTrace(0), dDiff(0), dChiSquare(0);
	int iSample(0);
	for (int i = 0; i < iEventlength; i++) {
		dTrace = event->Trace(i);
		if (dTrace == 0) continue;
		iSample = i - iOff;
		dVal = dBase;
		if ((iSample > -1) && (iSample < iStdLength)) dVal += dPeak*dStdWave[p][iSample];
		dDiff = dTrace-dVal;
		dChiSquare += (dDiff*dDiff)/(dTrace*dTrace);
	}
	return dChiSquare;
}

void XSQ_NEW::Analyze() {
	int iter(iIterations);
	double dEpsilonP(0.0005), dEpsilonB(0.0001), dPeak(event->Peak_y()), dBase(event->Baseline()); // looking in a small region
	double dDet(0), a(0), b(0), c(0), d(0), e(0), f(0), dDiff(0);
	double d000(0), dp00(0), d0p0(0), d00p(0), dpp0(0), dp0p(0), d0pp(0), dm00(0), d0m0(0), d00m(0);
	double dChiSquareLast(0), dChiSquare(0), dEpsilonO(1), iOff(0);
	SetDefaultParameters();
	XSQ_NEW::ssIterations[0][id] = XSQ_NEW::ssIterations[1][id] = 0;
	if (dBase - dPeak < iPeakCut) { // cut noise to save processing time
		return;
	}
	for (int p = 0; p < P; p++) {
		iter = iIterations; // in case extra time was used
		dPeak = (event->Baseline() - event->Peak_y())*dStdNorm[p];
		dBase = event->Baseline();
		iOff = event->Trigger() - iStdTrig;
		dChiSquare = FindChisquare(p, dPeak, dBase, iOff);
		for (XSQ_NEW::ssIterations[p][id] = 1; XSQ_NEW::ssIterations[p][id] <= iter; XSQ_NEW::ssIterations[p][id]++) {
			// x_n+1 = x_n - (main step size)(inverse Hessian matrix)(gradient of function at x_n)
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
				break;
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
				break; // break from iterations
			} else if (XSQ_NEW::ssIterations[p][id] == iIterations) { // more time to converge, only once
				iter += iIterations; 
			}
		} // iterations loop
		XSQ_NEW::sdXsq[p][id]			= dChiSquare;
		XSQ_NEW::sdPeakheight[p][id]	= dPeak*fGain[p]/iResolutionScale;
		XSQ_NEW::sdBaseline[p][id]		= dBase;
		XSQ_NEW::sdOffset[p][id]		= iOff;
		
		XSQ_NEW::sdPeak_err[p][id]		= dStep[0] > 0 ? dStep[0] : -dStep[0];
		XSQ_NEW::sdBase_err[p][id]		= dStep[1] > 0 ? dStep[1] : -dStep[1];
		XSQ_NEW::sdOff_err[p][id]		= dStep[2] > 0 ? dStep[2] : -dStep[2];
		
		XSQ_NEW::sdConvergence[p][id]	= dDiff/dChiSquareLast;
	} // p loop
}
