#include "Event.h"

float Event::sfVersion = 1.31;

Event::Event() : usTrace(nullptr), itBegin(nullptr), itEnd(nullptr), ciChan(-1) {
	if (g_verbose > 1) cout << "Event c'tor\n";
}

Event::Event(int eventlength, int baselength, int average, int threshold, int chan, unsigned short* usStart, double* dStart) : usTrace(usStart), itBegin(dStart), itEnd(dStart+eventlength-average), ciChan(chan) {
	if (g_verbose > 1) cout << "Event " << chan << " c'tor\n";
	if ((ciChan >= MAX_CH) || (ciChan < 0)) {
		cout << error_message[method_error] << "Channel\n";
		iFailed = 1;
		return;
	}
	iFailed = 0;
	iEventlength = eventlength;
	iBaselength = baselength;
	iAverage = average;
	iThreshold = threshold;
	try {
		vPeakCandidates.reserve(16);
		vPeaks.reserve(16);
	}
	catch (exception& e) {
		cout << error_message[alloc_error] << "Vector\n";
		iFailed = 1;
		return;
	}
}

Event::~Event() {
	if (g_verbose > 1) cout << "Event " << ciChan << " d'tor\n";
}

void Event::Analyze() {
	auto it = itBegin, itt = itBegin;
	for (auto itU = usTrace; itU < usTrace+iEventlength && it < itEnd; itU++, it++) *it = *itU; // copies from unsigned short buffer into double buffer
	*dBaseline		= 0.;
	*dBaseSigma		= 0.;
	*dBasePost		= 0.;
	*dBasePostSigma	= 0.;
	*dIntegral		= 0.;
	dPeak0->clear();
	dPeak2->clear();
	*dBasePeakP		= 0.;
	*dBasePeakN		= 0.;

	*sDecay			= 0;
	sRise->clear();
	sPeakX->clear();
	*sTrigger		= 0;
	sHWHM->clear();
	*sSatDur		= 0;
	*sNumPeaks		= 0;

	*bSaturated		= false;
	*bFullWaveform	= true;

	itBasePkP = itBasePkN = itSatEnd = itBegin;

	for (it = itBegin, itt = itEnd-1; it-itBegin < iBaselength; it++, itt--) {
		*dBaseline += *it; // baseline at start of event
		*dBaseSigma += (*it)*(*it);
		if (*itBasePkP < *it) itBasePkP = it;
		if (*itBasePkN > *it) itBasePkN = it;
		*dBasePost += *itt; // baseline at end of event
		*dBasePostSigma += (*itt)*(*itt);
	}
	*dBaseline /= iBaselength;
	*dBasePost /= iBaselength;
	*dBasePeakN = (*dBaseline - *itBasePkN)*dScaleV;
	*dBasePeakP = (*itBasePkP - *dBaseline)*dScaleV;

	*dBaseSigma = sqrt(*dBaseSigma/iBaselength-(*dBaseline)*(*dBaseline)); // sigma^2 = <x^2> - <x>^2
	*dBasePostSigma = sqrt(*dBasePostSigma/iBaselength-(*dBasePost)*(*dBasePost));
	Peakfinder();
	if (vPeaks.size() == 0) {
		sPeakX->push_back(0);
		sRise->push_back(0);
		sHWHM->push_back(0);
		dPeak0->push_back(0);
		dPeak2->push_back(0);
		return;
	} else if (*(vPeaks.front().itPeak) == 0) { // saturated event
		*bSaturated = true;
		for (it = vPeaks.front().itPeak; it < itEnd; it++) if (*it != 0) break;
		itSatEnd = it-1;
		*sSatDur = (itSatEnd-vPeaks.front().itPeak)*dScaleT;
	} else {
		itSatEnd = vPeaks.front().itPeak;
	}
	if (iAverage) {
		Average();
		Peakfinder();
	}
	if (vPeaks.size() == 0) {
		sPeakX->push_back(0);
		sRise->push_back(0);
		sHWHM->push_back(0);
		dPeak0->push_back(0);
		dPeak2->push_back(0);
		return;
	}
	for (it = vPeaks.front().itPeak; it > itBegin; it--) if ((*it < iThreshold) && (*(it-1) >= iThreshold)) break;
	*sTrigger = (it - itBegin)*dScaleT;
	for (it = vPeaks.front().itPeak; it < itEnd; it++) if (*it > *dBaseline-3*(*dBaseSigma)) break;
	*bFullWaveform = (it != itEnd);
	vPeaks.front().itEnd = (it == itEnd ? itEnd-1 : it);
	for (auto& p : vPeaks) { //(auto it = vPeaks.begin(); it < vPeaks.end(); it++) {
		sPeakX->push_back((p.itPeak-itBegin)*dScaleT);
		sRise->push_back((p.itPeak-p.itStart)*dScaleT);
		sHWHM->push_back(p.HWHM()*dScaleT);
		dPeak0->push_back((*(p.itStart)-*(p.itPeak))*dScaleV);
		if (iAverage) dPeak2->push_back(dPeak0->back()); // if the waveform is averaged there's no need to do peak averaging
		else {
			dPeak2->push_back(0);
			if ((p.itPeak - itBegin > 1) && (itEnd - p.itPeak > 2)) {
				for (itt = p.itPeak - 2; itt <= p.itPeak + 2; itt++) dPeak2->back() += *itt;
				dPeak2->back() = (*(p.itStart) - 0.2*dPeak2->back())*dScaleV;
			} else {
				dPeak2->back() = dPeak0->back();
			}
		}
	}

	*sNumPeaks = vPeaks.size();
	*sDecay = (vPeaks.front().itEnd - vPeaks.front().itPeak)*dScaleT;
	for (it = vPeaks.front().itStart; it <= vPeaks.front().itEnd; it++) *dIntegral += *it; // integrator
	*dIntegral = (*dIntegral)*2 - (*vPeaks.front().itStart + *vPeaks.front().itEnd);
	*dIntegral = ((*dBaseline)*(vPeaks.front().itEnd-vPeaks.front().itStart) - 0.5*(*dIntegral))*dScaleT*dScaleV;
}

inline void Event::Average() {
	auto itU = usTrace + (itSatEnd-itBegin) + (itSatEnd-itBegin < iAverage ? iAverage - (itSatEnd - itBegin): 0);
	for (auto it = itSatEnd + (itSatEnd-itBegin < iAverage ? iAverage - (itSatEnd - itBegin): 0); (it < itEnd) && (itU < usTrace + iEventlength - iAverage); it++, itU++) {
		*it = 0.;
		for (auto itt = itU-iAverage; itt <= itU+iAverage; itt++) *it += *itt;
		*it /= (2.*iAverage+1.);
	}
}

void Event::Peakfinder() {
	int iPeakCutBL(8/(1000.*dScaleV)), iPeakCutMin(16/(1000.*dScaleV)); // ignores Peaks less than 8 or 16 mV depending on if the pulse reaches baseline or not
	vPeakCandidates.clear();
	vPeaks.clear();
	Peak_t peak(itBegin);
	auto it = itBegin, itt = itBegin;

	for (it = itBegin; it < itEnd; it += iBaselength/2) { // peak-finding logic
		peak.itPeak = it;
		for (itt = it; itt < it + iBaselength; itt++) { // finds minimum value in window
			if (itt == itEnd) break;
			if (*peak.itPeak > *itt) {
				peak.itPeak = itt;
			}
		}
		if (peak.itPeak - it < iBaselength/4) {
			continue; // in the first quarter of the window, falling edge or nothing.
		}
		else if (peak.itPeak - it < 3*iBaselength/4) { // not in the end of the window
			if (*(peak.itPeak) < iThreshold) {
				vPeakCandidates.push_back(peak); // cuts noise
			} else {
				continue;
			}
		} else { // in last quarter of window
			continue;
		}
	}

	if (vPeakCandidates.size() == 0) { // no Peaks found, use default values
		return;
	} else { // peaks found, select candidates of sufficient height
		auto itPrev = itBegin, itMin = itBegin;
		for (auto& p : vPeakCandidates) { // evaluates heights of peak candidates
			itMin = p.itPeak;
			for (it = p.itPeak; it > itPrev; it--) { // looks for start of peak or minimum value
				if (*itMin < *it) itMin = it;
				if (*it > (*dBaseline)-3*(*dBaseSigma)) break;
			}
			if (it == itPrev) { // didn't reach baseline
				if ((*itMin - *(p.itPeak) > iPeakCutMin)){
					p.itStart = itMin;
					vPeaks.push_back(p);
				}
			} else { // did reach baseline
				if (((*dBaseline) - *(p.itPeak) > iPeakCutBL)) {
					p.itStart = it;
					vPeaks.push_back(p);
				}
			} // choosing it or itMin
		} // iter loop
	}
	stable_sort(vPeaks.begin(), vPeaks.end(), [&] (const Peak_t& lhs, const Peak_t& rhs){return (*lhs.itStart-*lhs.itPeak) > (*rhs.itStart-*rhs.itPeak);}); // sorts descending via lambda

/*	for (unsigned i = 0; i < vPeaks.size(); i++) { // sorts peaks by descending size
		auto BigPeak = vPeaks.begin() + i;
		Peak_t temp = nullptr;
		for (auto iter = vPeaks.begin() + i; iter < vPeaks.end(); iter++) {
			if (*iter > *BigPeak) BigPeak = iter;
		}
		if (BigPeak == vPeaks.begin() + i) continue;
		else {
			temp = *(vPeaks.begin() + i);
			*(vPeaks.begin() + i) = *BigPeak;
			*BigPeak = temp;
		}
	} */
	return;
}

void Event::SetAddresses(vector<void*> add) {
	int i(0);
	bFullWaveform	= (bool*)add[i++];
	bSaturated		= (bool*)add[i++];

	sDecay			= (short*)add[i++];
	sRise			= (vector<double>*)add[i++];
	sPeakX			= (vector<double>*)add[i++];
	sTrigger		= (short*)add[i++];
	sHWHM			= (vector<double>*)add[i++];
	sSatDur			= (short*)add[i++];
	sNumPeaks		= (short*)add[i++];

	dBaseline		= (double*)add[i++];
	dBaseSigma		= (double*)add[i++];
	dBasePost		= (double*)add[i++];
	dBasePostSigma	= (double*)add[i++];
	dPeak0			= (vector<double>*)add[i++];
	dPeak2			= (vector<double>*)add[i++];
	dIntegral		= (double*)add[i++];
	dBasePeakP		= (double*)add[i++];
	dBasePeakN		= (double*)add[i++];
}
