#include "Method.h"
#include "TFile.h"
#include "TGraphErrors.h"

using std::cout;

const auto pi = acos(-1.0);

float Method::sfVersion = 1.28;

Method::Method(const int length, const int fast, const int slow, const int samples, const array<float,2>& gain, const double scaleT, const double scaleV, shared_ptr<Event> ev, std::string sConfDir) : dSlow(0.01), dShigh(1.0) , sConfigDir(sConfDir) {
	if (g_verbose > 1) cout << "Method c'tor\n";
	iFailed = 0;
	iEventlength = length;
	iFastTime = fast;
	iSlowTime = slow;
	iPGASamples = samples;
	fGain[n] = gain[n];
	fGain[y] = gain[y];
	dNsPerSample = scaleT;
	dVoltsPerBin = scaleV;
	event = ev;
	iPGAAverage = event->GetAverage() > 0 ? 0 : 5; // no sense averaging averages
	iLAPAverage = event->GetAverage() > 0 ? 0 : 4;

	unique_ptr<TFile> std_file = nullptr;
	TGraphErrors* pWave = nullptr;

	double omega;
	try {
		dCos.fill(vector<double>(iEventlength));
		dSin.fill(vector<double>(iEventlength));
		dTrace.reserve(iEventlength);
		dExp.fill(vector<double>(iEventlength));
		dX.reserve(iEventlength);
	} catch (std::bad_alloc& ba) {
		cout << error_message[alloc_error] << "Method\n";
		iFailed = 1;
		return;
	}
	for (auto n = 0; n < ciDFTOrder; n++) {
		omega = 2.*n*pi/(iEventlength*dNsPerSample); // GHz or /ns
		for (auto t = 0; t < iEventlength; t++) {
			dCos[n][t] = cos(omega*t);
			dSin[n][t] = sin(omega*t); // simpler than using one table and sin(x) = cos(x-pi/2)
	}	}

	double dScale(log(dShigh/dSlow)/ciLAPNpts);

	for (auto n = 0; n < ciLAPNpts; n++) {
		dS[n] = dSlow*exp(dScale*n);
		for (auto i = 0; i < iEventlength; i++) dExp[n][i] = exp(-dS[n]*i*dNsPerSample);
	}

	if (dNsPerSample == 1) { // ns/Sa
		iStdLength	= 600;
		iStdTrig	= 50;
	} else if (dNsPerSample == 2) {
		iStdLength	= 300;
		iStdTrig	= 25;
	} else if (dNsPerSample == 0.5) {
		iStdLength		= 1199;
		iStdTrig		= 100;
	} else {
		cout << error_message[dig_error];
		cout << error_message[method_error] << "Std Events\n";
		iFailed = 1;
		return;
	}
	if (dVoltsPerBin == 1./1024) { // volts/bin
		iResolutionScale = 1;
	} else if (dVoltsPerBin == 2./(1 << 14)) {
		iResolutionScale = 1 << 3;
	} else {
		cout << error_message[dig_error];
		cout << error_message[method_error] << "Std Events\n";
		iFailed = 1;
		return;
	}

	try {std_file.reset(new TFile((sConfigDir + "/config/standard_events.root").c_str(), "READ"));}
	catch (std::bad_alloc& ba) {
		cout << error_message[alloc_error] << "Std Events file\n";
		iFailed = 1;
		return;
	}
	if (std_file->IsZombie()) {
		cout << error_message[file_error] << "Std Events\n";
		iFailed = 1;
		return;
	}
	for (int p = 0; p < P; p++) {
		try {dStdWave[p].reserve(iStdLength);}
		catch (std::bad_alloc& ba) {
			cout << error_message[alloc_error] << "Std Events\n";
			iFailed = 1;
			return;
		}
		pWave = (TGraphErrors*)std_file->Get((p ? "gamma_low" : "neutron_low"));
		if (pWave == nullptr) {
			cout << error_message[root_error] << "Std Events\n";
			iFailed = 1;
			return;
		}
		dStdPeak[p] = 0;
		switch(iStdLength) {
			case 300 : // 500 MSa/s
				for (int i = 0; i < iStdLength; i++) { // averages
					dStdWave[p][i] = (pWave->GetY()[2*i] + pWave->GetY()[2*i+1])/2.;
					dStdPeak[p] = std::min(dStdPeak[p], dStdWave[p][i]);
				} break;
			case 1199 : // 2 GSa/s
				for (int i = 0; i < iStdLength; i++) { // interpolates
					dStdWave[p][i] = (i%2) ? (pWave->GetY()[(i+1)/2] + pWave->GetY()[(i-1)/2])/2. : pWave->GetY()[i/2]; // i%2==1 so i/2 = (i-1)/2
					dStdPeak[p] = std::min(dStdPeak[p], dStdWave[p][i]);
				} break;
			case 600 : // 1 GSa/s
				for (int i = 0; i < iStdLength; i++) {
					dStdWave[p][i] = pWave->GetY()[i];
					dStdPeak[p] = std::min(dStdPeak[p], dStdWave[p][i]);
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
		fit_n	= unique_ptr<TF1>(new TF1("fitn",this,&Method::TF1_fit_func,0,iEventlength,4));
		fit_y	= unique_ptr<TF1>(new TF1("fity",this,&Method::TF1_fit_func,0,iEventlength,4));
		fit_n_f	= unique_ptr<TF1>(new TF1("fitnf",this,&Method::TF1_fit_func,0,iEventlength,4));
		fit_y_f	= unique_ptr<TF1>(new TF1("fityf",this,&Method::TF1_fit_func,0,iEventlength,4));
	} catch (std::bad_alloc& ba) {
		cout << error_message[alloc_error] << "Fitter\n";
		iFailed = 1;
		return;
	}
	fit_n->SetParNames("Peakscale","Baseline","Offset","particle");
	fit_y->SetParNames("Peakscale","Baseline","Offset","particle");
	fit_n_f->SetParNames("Peakscale","Baseline","Offset","particle");
	fit_y_f->SetParNames("Peakscale","Baseline","Offset","particle");

	fit_n->FixParameter(3, n);
	fit_y->FixParameter(3, y);
	fit_n_f->FixParameter(3, n);
	fit_y_f->FixParameter(3, y);

	for (auto i = 0; i < iEventlength; i++) dX[i] = i;
}

Method::~Method() {
	if (g_verbose > 1) cout << "Method d'tor\n";
	fit_n.reset();
	fit_y.reset();
	fit_n_f.reset();
	fit_y_f.reset();
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

	*bTruncated = (event->vPeaks.size() > 0) && ((event->vPeaks.front().itStart + iSlowTime) >= event->itEnd);

	//CCM
	*dBaseline		= (*(event->dBaseline)-dZero) * dVoltsPerBin;
	*dBaseSigma		= *(event->dBaseSigma) * dVoltsPerBin;
	*dBasePost		= (*(event->dBasePost)-dZero) * dVoltsPerBin;
	*dBasePostSigma	= *(event->dBasePostSigma) * dVoltsPerBin;
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
	*dPeakheight_n	= (event->dPeak0->front()/dVoltsPerBin)*dStdNorm[n]*fGain[n];
	*dBaseline_n	= *(event->dBaseline);
	*dOffset_n		= *(event->sTrigger) - iStdTrig;

	*dXsq_y			= -1;
	*dPeakheight_y	= (event->dPeak0->front()/dVoltsPerBin)*dStdNorm[y]*fGain[y];
	*dBaseline_y	= *(event->dBaseline);
	*dOffset_y		= *(event->sTrigger) - iStdTrig;

	*dPeak_err_n	= -1;
	*dPeak_err_y	= -1;
	*dBase_err_n	= -1;
	*dBase_err_y	= -1;
	*dOff_err_n		= -1;
	*dOff_err_y		= -1;

	*dXsq_n_f		= -1;
	*dPeakheight_n_f= (event->dPeak0->front()/dVoltsPerBin)*dStdNorm[n]*fGain[n];
	*dOffset_n_f	= *(event->sTrigger) - iStdTrig;

	*dXsq_y_f		= -1;
	*dPeakheight_y_f= (event->dPeak0->front()/dVoltsPerBin)*dStdNorm[y]*fGain[y];
	*dOffset_y_f	= *(event->sTrigger) - iStdTrig;

	*dPeak_err_n_f	= -1;
	*dPeak_err_y_f	= -1;
	*dOff_err_n_f	= -1;
	*dOff_err_y_f	= -1;

	return;
}

void Method::Analyze() {
	auto i(0);
	auto dReal(0.), dImag(0.);
	auto m(0), t(0), iFast(0), iSlow(0);

	SetDefaultValues();
	if (event->vPeaks.size() == 0) return; // noise, nothing to do here
	//CCM
	iFast = (event->vPeaks.front().itStart + iFastTime - 1 < event->itEnd ? iFastTime : event->itEnd - event->vPeaks.front().itStart - 1);
	iSlow = (event->vPeaks.front().itStart + iSlowTime - 1 < event->itEnd ? iSlowTime : event->itEnd - event->vPeaks.front().itStart - 1); // local integration limits

	for (auto it = event->vPeaks.front().itStart; it <= event->vPeaks.front().itStart + iFast; it++) *dFastInt += *it;
	*dSlowInt = *dFastInt;
	for (auto it = event->vPeaks.front().itStart + iFast+1; it <= event->vPeaks.front().itStart + iSlow; it++) *dSlowInt += *it;

	*dFastInt *= 2.; // trapezoid rule: integral ~ f(0) + 2*f(1) + ... + 2*f(n-1) + f(n)
	*dSlowInt *= 2.; // faster to do sum f(i), double, and subtract endpoints
	*dFastInt -= (*event->vPeaks.front().itStart + *(event->vPeaks.front().itStart + iFast));
	*dSlowInt -= (*event->vPeaks.front().itStart + *(event->vPeaks.front().itStart + iSlow));

	*dSlowInt = ((*event->dBaseline * (iSlow)) - 0.5 * (*dSlowInt)) * dVoltsPerBin * dNsPerSample; // baseline subtraction
	*dFastInt = ((*event->dBaseline * (iFast)) - 0.5 * (*dFastInt)) * dVoltsPerBin * dNsPerSample;

#ifndef CCM_ONLY

	//PGA
	if ((event->vPeaks.front().itPeak + iPGASamples + iPGAAverage) < event->itEnd) {
		for (i = -iPGAAverage; i <= iPGAAverage; i++) *dSample += *(event->vPeaks.front().itPeak + iPGASamples + i); // average to reduce statistical fluctuations
		*dSample /= (2.*iPGAAverage + 1);
	} else *dSample = -1;

	//DFT
	for (m = 0; m < ciDFTOrder; m++) {
		dReal = 0;
		dImag = 0;
		t = 0;
		for (auto it = event->itBegin; it < event->itEnd; it++, t++) { // fourier series
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
		t = 0;
		for (auto itD = event->itSatEnd; itD < event->itEnd-iLAPAverage; itD++, t++) { // performs the 9pt moving average over the trailing edge
			dTrace[t] = 0;
			for (auto tt = -iLAPAverage; tt <= iLAPAverage; tt++) dTrace[t] += *(itD+tt);
			dTrace[t] /= (2.*iLAPAverage + 1.);
		}
		for (; t < iEventlength; t++) dTrace[t] = *(event->dBasePost);
	} else {
		t = 0;
		for (auto itD = event->itSatEnd; itD < event->itEnd; itD++, t++) dTrace[t] = *itD;
		for (; t < iEventlength; t++) dTrace[t] = *(event->dBasePost);
	}
	for (m = 0; m < ciLAPNpts; m++) {
		dXform[m] = 0;
		t = 0;
	//	for (auto it = dTrace.data(); it < dTrace.data()+iEventlength; it++, t++) dXform[m] += (*it)*dExp[m][t];
		for (auto& it : dTrace) dXform[m] += it*dExp[m][t++];
	}
	for (m = 0; m < ciLAPNpts-1; m++) {
		(m < (ciLAPNpts >> 1) ? *dLaplaceLow : *dLaplaceHigh) += (dS[m+1]-dS[m])*(dXform[m+1]+dXform[m]);
	}
	*dLaplaceHigh *= 0.5*dVoltsPerBin;
	*dLaplaceLow *= 0.5*dVoltsPerBin;

	//NGM
	try {
		graph.reset(new TGraph(event->Length(), dX.data(), event->itBegin));
	} catch (std::bad_alloc& ba) {
		return;
	}
	fit_n->SetParameter(0, *dPeakheight_n);
	fit_n->SetParameter(1, *(event->dBaseline));
	fit_n->SetParameter(2, *dOffset_n);
	graph->Fit(fit_n.get(), "Q N"); // quiet, no-plot

	*dXsq_n			= fit_n->GetChisquare();
	*dPeakheight_n	= fit_n->GetParameter(0)/fGain[n]/iResolutionScale;
	*dBaseline_n	= fit_n->GetParameter(1);
	*dOffset_n		= fit_n->GetParameter(2);

	*dPeak_err_n	= fit_n->GetParError(0)/fGain[n]/iResolutionScale;
	*dBase_err_n	= fit_n->GetParError(1);
	*dOff_err_n		= fit_n->GetParError(2);

	fit_y->SetParameter(0, *dPeakheight_y);
	fit_y->SetParameter(1, *(event->dBaseline));
	fit_y->SetParameter(2, *dOffset_y);
	graph->Fit(fit_y.get(), "Q N");

	*dXsq_y			= fit_y->GetChisquare();
	*dPeakheight_y	= fit_y->GetParameter(0)/fGain[y]/iResolutionScale;
	*dBaseline_y	= fit_y->GetParameter(1);
	*dOffset_y		= fit_y->GetParameter(2);

	*dPeak_err_y	= fit_y->GetParError(0)/fGain[y]/iResolutionScale;
	*dBase_err_y	= fit_y->GetParError(1);
	*dOff_err_y		= fit_y->GetParError(2);

	fit_n_f->SetParameter(0, *dPeakheight_n_f); // Fixed baseline fit
	fit_n_f->FixParameter(1, *(event->dBaseline));
	fit_n_f->SetParameter(2, *dOffset_n_f);
	graph->Fit(fit_n_f.get(), "Q N"); // quiet, no-plot

	*dXsq_n_f		= fit_n_f->GetChisquare();
	*dPeakheight_n_f= fit_n_f->GetParameter(0)/fGain[n]/iResolutionScale;
	*dOffset_n_f	= fit_n_f->GetParameter(2);

	*dPeak_err_n_f	= fit_n_f->GetParError(0)/fGain[n]/iResolutionScale;
	*dOff_err_n_f	= fit_n_f->GetParError(2);

	fit_y_f->SetParameter(0, *dPeakheight_y_f);
	fit_y_f->FixParameter(1, *(event->dBaseline));
	fit_y_f->SetParameter(2, *dOffset_y_f);
	graph->Fit(fit_y_f.get(), "Q N"); // quiet, no-plot

	*dXsq_y_f		= fit_y_f->GetChisquare();
	*dPeakheight_y_f= fit_y_f->GetParameter(0)/fGain[y]/iResolutionScale;
	*dOffset_y_f	= fit_y_f->GetParameter(2);

	*dPeak_err_y_f	= fit_y_f->GetParError(0)/fGain[y]/iResolutionScale;
	*dOff_err_y_f	= fit_y_f->GetParError(2);

#endif
}

void Method::SetAddresses(const vector<void*>& add) {
	int i(0);
	bTruncated =		(bool*)add[i++];

	dBaseline =			(double*)add[i++];
	dBaseSigma =		(double*)add[i++];
	dBasePost =			(double*)add[i++];
	dBasePostSigma =	(double*)add[i++];
	dFastInt =			(double*)add[i++];
	dSlowInt =			(double*)add[i++];

	//PGA
	dSample =			(double*)add[i++];

	//DFT
	dOdd =				(double*)add[i++];
	dEven =				(double*)add[i++];

	//LAP
	dLaplaceHigh =		(double*)add[i++];
	dLaplaceLow =		(double*)add[i++];

	//NGM
	dXsq_n =			(double*)add[i++];
	dXsq_y =			(double*)add[i++];
	dXsq_n_f =			(double*)add[i++];
	dXsq_y_f =			(double*)add[i++];
	dPeakheight_n =		(double*)add[i++];
	dPeakheight_y =		(double*)add[i++];
	dPeakheight_n_f =	(double*)add[i++];
	dPeakheight_y_f =	(double*)add[i++];
	dBaseline_n =		(double*)add[i++];
	dBaseline_y =		(double*)add[i++];
	dOffset_n =			(double*)add[i++];
	dOffset_y =			(double*)add[i++];
	dOffset_n_f =		(double*)add[i++];
	dOffset_y_f =		(double*)add[i++];
	dPeak_err_n =		(double*)add[i++];
	dPeak_err_y =		(double*)add[i++];
	dPeak_err_n_f =		(double*)add[i++];
	dPeak_err_y_f =		(double*)add[i++];
	dBase_err_n =		(double*)add[i++];
	dBase_err_y =		(double*)add[i++];
	dOff_err_n =		(double*)add[i++];
	dOff_err_y =		(double*)add[i++];
	dOff_err_n_f =		(double*)add[i++];
	dOff_err_y_f =		(double*)add[i++];
}
