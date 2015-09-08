#include "Method.h"
#include "TFile.h"
#include "TVectorT.h"

const double pi = acos(-1.0);

float Method::sfVersion = 1.0;

Method::Method(int length, int fast, int slow, int samples, float gain[2], double scaleT, double scaleV, shared_ptr<Event> ev) : dSlow(0.01), dShigh(1.0) {
	if (g_verbose) cout << "Method c'tor\n";
	iFailed = 0;
	iEventlength = length;
	iFastTime = fast;
	iSlowTime = slow;
	iPGASamples = samples;
	fGain[n] = gain[n];
	fGain[y] = gain[y];
	dScaleT = scaleT;
	dScaleV = scaleV;
	event = ev;
	iPGAAverage = event->GetAverage() > 0 ? 0 : 5; // no sense averaging averages
	iLAPAverage = event->GetAverage() > 0 ? 0 : 4;

	unique_ptr<TFile> std_file = nullptr;
	TVectorT<double>* pWave = nullptr;

	double omega;
	for (auto n = 0; n < ciDFTOrder; n++) {
		try {
			dCos[n].reserve(iEventlength);
			dSin[n].reserve(iEventlength);
		} catch (bad_alloc& ba) {cout << error_message[alloc_error] << "DFT lookup\n"; return;}
		omega = 2*n*pi/(iEventlength*dScaleT); // GHz
		for (auto t = 0; t < iEventlength; t++) {
			dCos[n].push_back(cos(omega*t));
			dSin[n].push_back(sin(omega*t)); // simpler than using one table and sin(x) = cos(x-pi/2)
	}	}

	double dScale(log(dShigh/dSlow)/ciLAPNpts);
	try {
		dTrace.reserve(iEventlength);
		dS.reserve(ciLAPNpts);
		dXform.reserve(ciLAPNpts);
	} catch (bad_alloc& ba) {
		cout << error_message[alloc_error] << "LAP lookup\n";
		return;
	}
	for (auto n = 0; n < ciLAPNpts; n++) {
		try {
			dExp[n].reserve(iEventlength);
		} catch (bad_alloc& ba) {
			cout << error_message[alloc_error] << "LAP lookup\n";
			return;
		}
		dS[n] = dSlow*exp(dScale*n);
		for (auto i = 0; i < iEventlength; i++) dExp[n].push_back(exp(-dS[n]*i*dScaleT));
	}

	int i(0), p(0);
	if (dScaleT == 1) { // ns/Sa
		iStdLength	= 450;
		iStdTrig	= 64;
	} else if (dScaleT == 2) {
		iStdLength	= 225;
		iStdTrig	= 32;
	} else if (dScaleT == 0.5) {
		iStdLength		= 899;
		iStdTrig		= 128;
	} else {
		cout << error_message[dig_error];
		cout << error_message[method_error] << "Std Events\n";
		return;
	}
	if (dScaleV == 1./1024) { // volts/bin
		iPeakCut = 4;
		iResolutionScale = 1;
	} else if (dScaleV == 2./(1 << 14)) { // TODO fix for simulated resolutions
		iPeakCut = 32;
		iResolutionScale = 1 << 3;
	} else {
		cout << error_message[dig_error];
		cout << error_message[method_error] << "Std Events\n";
		return;
	}

	try {std_file.reset(new TFile((sWorkingDir + "/Data_Processor/config/standard_events.root").c_str(), "READ"));}
	catch (bad_alloc& ba) {cout << error_message[alloc_error] << "Std Events\n"; return;}
	if (!std_file->IsOpen()) {cout << error_message[file_error] << "Std Events\n"; return;}
	for (p = 0; p < P; p++) {
		try {dStdWave[p].reserve(iStdLength);}
		catch (bad_alloc& ba) {cout << error_message[alloc_error] << "Std Events\n"; return;}
		pWave = (TVectorT<double>*)std_file->Get((p ? "gamma_wave_inv" : "neutron_wave_inv"));
		if (pWave == nullptr) {cout << error_message[root_error] << "Std Events\n"; return;}
		dStdPeak[p] = 1000;
		switch(iStdLength) {
			case 225 : // 500 MSa/s
				for (i = 0; i < iStdLength; i++) { // averages
					dStdWave[p].push_back(((*pWave)[2*i] + (*pWave)[2*i+1])/2.);
					dStdPeak[p] = min(dStdPeak[p], dStdWave[p][i]);
				} break;
			case 899 : // 2 GSa/s
				for (i = 0; i < iStdLength; i++) { // interpolates
					dStdWave[p].push_back((i%2) ? ((*pWave)[(i+1)/2] + (*pWave)[(i-1)/2])/2. : (*pWave)[i/2]); // i%2==1 so i/2 = (i-1)/2
					dStdPeak[p] = min(dStdPeak[p], dStdWave[p][i]);
				} break;
			case 450 : // 1 GSa/s
				for (i = 0; i < iStdLength; i++) {
					dStdWave[p].push_back((*pWave)[i]);
					dStdPeak[p] = min(dStdPeak[p], dStdWave[p][i]);
				} break;
			default : cout << error_message[method_error];
			return;
		}
		dStdNorm[p] = -1./dStdPeak[p]*(iResolutionScale); // waves generated on 10-bit digitizer so must be scaled appropriately for other resolutions
	} // p
	std_file->Close();
	std_file = nullptr;
	pWave = nullptr;
	graph = nullptr;

	try {
		dX.reserve(iEventlength);
		fit = unique_ptr<TF1>(new TF1("fit",this,&Method::TF1_fit_func,0,iEventlength,4));
	} catch (bad_alloc& ba) {
		cout << error_message[alloc_error] << "Fitter\n";
		return;
	}
	fit->SetParNames("Peakscale","Baseline","Offset","particle");
	for (i = 0; i < iEventlength; i++) dX.push_back(i);
}

Method::~Method() {
	if (g_verbose) cout << "Method d'tor\n";
	fit.reset();
	graph.reset();
	event.reset();
}

double Method::TF1_fit_func(double* x, double* par) {
	int iSample(x[0]-par[2]), p(par[3]);
	double dVal(par[1]);
	if ((iSample > -1) && (iSample < iStdLength)) dVal += par[0]*dStdWave[p][iSample];
	if (dVal < 0) dVal = 0; // saturated events
	return dVal;
}

void Method::SetDefaultValues() {

	*bTruncated = ((event->Peak.itStart + iSlowTime) >= event->itEnd);

	//CCM
	*dBaseline		= *(event->dBaseline) * dScaleV;
	*dBaseSigma		= *(event->dBaseSigma) * dScaleV;
	*dBasePost		= *(event->dBasePost) * dScaleV;
	*dBasePostSigma	= *(event->dBasePostSigma) * dScaleV;
	*dBasePeakP		= (*(event->itBasePkP) - *(event->dBaseline))* dScaleV;
	*dBasePeakN		= (*(event->dBaseline) - *(event->itBasePkN))* dScaleV;
	*dPeak1			= (*(event->dBaseline) - *(event->Peak.itPeak))* dScaleV;
	*dPeak2			= *dPeak1;
	*dSlowInt		= 0;
	*dFastInt		= 0;

	//PGA
	*dSample		= 0;

	//DFT
	*dOdd			= 0;
	*dEven			= 0;

	//LAP
	*dLaplaceHigh	= 0;
	*dLaplaceLow	= 0;

	//TF1
	*dXsq_n			= -1;
	*dPeakheight_n	= (*(event->dBaseline) - *(event->dPeak0))*dStdNorm[n]*fGain[n];
	*dBaseline_n	= *(event->dBaseline);
	*dOffset_n		= *(event->sTrigger) - iStdTrig;

	*dXsq_y			= -1;
	*dPeakheight_y	= (*(event->dBaseline) - *(event->dPeak0))*dStdNorm[y]*fGain[y];
	*dBaseline_y	= *(event->dBaseline);
	*dOffset_y		= *(event->sTrigger) - iStdTrig;

	*dPeak_err_n	= -1;
	*dPeak_err_y	= -1;
	*dBase_err_n	= -1;
	*dBase_err_y	= -1;
	*dOff_err_n		= -1;
	*dOff_err_y		= -1;

	for (auto it = dTrace.begin(); it < dTrace.end(); it++) *it = 0;

	return;
}

void Method::Analyze() {
	auto i(0), iPGA_average(5);
	auto dTemp(0.), dReal(0.), dImag(0.);
	auto m(0), t(0), iFast(0), iSlow(0);

	SetDefaultValues();

	//CCM
	iFast = (iFastTime < event->itEnd - 1 - event->Peak.itStart ? iFastTime : event->itEnd - 1 - event->Peak.itStart);
	iSlow = (iSlowTime < event->itEnd - 1 - event->Peak.itStart ? iSlowTime : event->itEnd - 1 - event->Peak.itStart); // local integration limits for fast and slow

	if (((event->Peak.itPeak + 2) < event->itEnd) && (event->Peak.itPeak - 1 > event->itBegin) && !*(event->bSaturated)) { // peak averaging
		for (auto it = event->Peak.itPeak-1; it < event->Peak.itPeak; it++) dTemp += *it;
		*dPeak1 = (*(event->dBaseline) - (0.333*dTemp))*dScaleV; // averaged with adjacent samples
		dTemp += *(event->Peak.itPeak-2) + *(event->Peak.itPeak + 2);
		*dPeak2 = (*(event->dBaseline) - (0.2*dTemp))*dScaleV; // averaged with two adjacent samples
	}

	for (auto it = event->Peak.itStart; it <= event->Peak.itStart + iFast; it++) *dFastInt += *it; // integrator
	*dSlowInt = *dFastInt;
	for (auto it = event->Peak.itStart + iFast; it <= event->Peak.itStart + iSlow; it++) *dSlowInt += *it;

	*dFastInt *= 2.; // trapezoid rule: integral = f(0) + 2*f(1) + ... + 2*f(n-1) + f(n)
	*dSlowInt *= 2.; // faster to do sum f(i), double, and subtract endpoints
	*dFastInt -= (*event->Peak.itStart + *(event->Peak.itStart + iFast));
	*dSlowInt -= (*event->Peak.itStart + *(event->Peak.itStart + iSlow));

	*dSlowInt = ((*event->dBaseline * (iSlow)) - 0.5 * (*dSlowInt)) * dScaleV * dScaleT; // baseline subtraction
	*dFastInt = ((*event->dBaseline * (iFast)) - 0.5 * (*dFastInt)) * dScaleV * dScaleT;

	//PGA
	if ((event->Peak.itPeak + iPGASamples + iPGA_average) < event->itEnd) {
		*dSample = 0;
		for (i = -iPGA_average; i <= iPGA_average; i++) *dSample += *(event->Peak.itPeak + iPGASamples + i); // average to reduce statistical fluctuations
		*dSample /= (2.*iPGA_average + 1);
	} else *dSample = -1;

	//DFT
	for (m = 0; m < ciDFTOrder; m++) {
		dReal = 0;
		dImag = 0;
		for (auto it = event->itBegin; it < event->itEnd; it++) { // fourier series
			dReal += (*it)*dCos[m][t];
			dImag += (*it)*dSin[m][t];
		}
		dMagnitude[m] = sqrt(dReal*dReal + dImag*dImag);
	}
	*dEven = *dOdd = 0;
	for (t = 2; t < ciDFTOrder; t++) (t%2 ? *dOdd : *dEven) += (dMagnitude[t]-dMagnitude[(t%2?1:0)])/dMagnitude[(t%2?1:0)];

	//LAP
	*dLaplaceLow = 0;
	*dLaplaceHigh = 0;

	if (iLAPAverage) {
		auto itA = dTrace.begin();
		for (auto itD = event->itSatEnd; itD < event->itEnd-iLAPAverage; itD++, itA++) { // performs the 9pt moving average
			*itA = 0;
			for (auto tt = -iLAPAverage; tt <= iLAPAverage; tt++) *itA += *(itD+tt);
			*itA /= (2.*iLAPAverage + 1.);
		}
	}
	for (m = 0; m < ciLAPNpts; m++) {
		dXform[m] = 0;
		t = 0;
		for (auto it = dTrace.begin(); it < dTrace.end(); it++, t++) dXform[m] += (*it)*dExp[m][t]; // trailing edge of the pulse
	}
	for (m = 0; m < ciLAPNpts-1; m++) {
		(m < (ciLAPNpts >> 1) ? *dLaplaceLow : *dLaplaceHigh) += (dS[m+1]-dS[m])*(dXform[m+1]+dXform[m]);
	}
	*dLaplaceHigh *= 0.5*dScaleV;
	*dLaplaceLow *= 0.5*dScaleV;

	//NGM
	double integral(*event->dIntegral); // gets around g++ optimization errors
	if (integral == 0) { // doesn't process if no trigger to save time
		return;
	}
	try {graph.reset(new TGraph(iEventlength, dX.data(), event->vTrace.data()));}
	catch (bad_alloc& ba) {
		return;
	}
	fit->SetParameter(0, (*event->dBaseline-*event->Peak.itPeak)*dStdNorm[n]);
	fit->SetParameter(1, *event->dBaseline);
	fit->SetParameter(2, *event->sTrigger - iStdTrig);
	fit->FixParameter(3, n);
	graph->Fit(fit.get(), "Q N"); // quiet, no-plot

	*dXsq_n			= fit->GetChisquare();
	*dPeakheight_n	= fit->GetParameter(0)/fGain[n]/iResolutionScale;
	*dBaseline_n	= fit->GetParameter(1);
	*dOffset_n		= fit->GetParameter(2);

	*dPeak_err_n	= fit->GetParError(0)/fGain[n]/iResolutionScale;
	*dBase_err_n	= fit->GetParError(1);
	*dOff_err_n		= fit->GetParError(2);


	fit->SetParameter(0, (*event->dBaseline - *event->Peak.itPeak)*dStdNorm[y]);
	fit->SetParameter(1, *event->dBaseline);
	fit->SetParameter(2, *event->sTrigger - iStdTrig);
	fit->FixParameter(3, y);
	graph->Fit(fit.get(), "Q N");

	*dXsq_y			= fit->GetChisquare();
	*dPeakheight_y	= fit->GetParameter(0)/fGain[y]/iResolutionScale;
	*dBaseline_y	= fit->GetParameter(1);
	*dOffset_y		= fit->GetParameter(2);

	*dPeak_err_y	= fit->GetParError(0)/fGain[y]/iResolutionScale;
	*dBase_err_y	= fit->GetParError(1);
	*dOff_err_y		= fit->GetParError(2);

}

void Method::SetAddresses(vector<void*> add) {
	int i(0);
	bTruncated = (bool*)add[i++];

	dBaseline = (double*)add[i++];
	dBaseSigma = (double*)add[i++];
	dBasePost = (double*)add[i++];
	dBasePostSigma = (double*)add[i++];
	dBasePeakN = (double*)add[i++];
	dBasePeakP = (double*)add[i++];
	dSlowInt = (double*)add[i++];
	dFastInt = (double*)add[i++];
	dPeak1 = (double*)add[i++];
	dPeak2 = (double*)add[i++];

	//PGA
	dSample = (double*)add[i++];

	//DFT
	dOdd = (double*)add[i++];
	dEven = (double*)add[i++];

	//LAP
	dLaplaceHigh = (double*)add[i++];
	dLaplaceLow = (double*)add[i++];

	//NGM
	dXsq_n = (double*)add[i++];
	dXsq_y = (double*)add[i++];
	dPeakheight_n = (double*)add[i++];
	dPeakheight_y = (double*)add[i++];
	dBaseline_n = (double*)add[i++];
	dBaseline_y = (double*)add[i++];
	dOffset_n = (double*)add[i++];
	dOffset_y = (double*)add[i++];
	dPeak_err_n = (double*)add[i++];
	dPeak_err_y = (double*)add[i++];
	dBase_err_n = (double*)add[i++];
	dBase_err_y = (double*)add[i++];
	dOff_err_n = (double*)add[i++];
	dOff_err_y = (double*)add[i++];
}