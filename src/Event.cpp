#include "Event.h"

float Event::sfVersion = 1.15;

Event::Event() : usTrace(nullptr), itBegin(nullptr), itEnd(nullptr), ciChan(-1) {
	if (g_verbose) cout << "Event c'tor\n";
}

Event::Event(int eventlength, int baselength, int average, int threshold, int chan, unsigned short* usStart, double* dStart) : usTrace(usStart), itBegin(dStart), itEnd(dStart+eventlength-average), ciChan(chan) {
	if (g_verbose) cout << "Event " << chan << " c'tor\n";
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
}

Event::~Event() {
	if (g_verbose) cout << "Event " << ciChan << " d'tor\n";
}

void Event::Analyze() {
	auto it = itBegin, itt = itBegin;
	for (auto itU = usTrace; itU < usTrace+iEventlength && it < itEnd; itU++, it++) *it = *itU; // copies from unsigned short buffer into double buffer
	*dBaseline		= 0.;
	*dBaseSigma		= 0.;
	*dBasePost		= 0.;
	*dBasePostSigma	= 0.;
	*dIntegral		= 0.;
	*dPeak0			= 0.;
	*dPeak0s		= 0.;
	*dBasePeakP		= 0.;
	*dBasePeakN		= 0.;

	*sDecay			= 0;
	*sRise			= 0;
	*sPeakX			= 0;
	*sPeakXs		= 0;
	*sTrigger		= 0;
	*sPeaks			= 0;

	*bSaturated		= false;
	*bFullWaveform	= true;

	itBasePkP = itBasePkN = itSatEnd = Peak.itPeak = Peak.itStart = itBegin;
	Peak.itEnd = itEnd;
	PeakS = itBegin;

	auto dTemp(0.);

	for (it = itBegin, itt = itEnd-1; it-itBegin < iBaselength; it++, itt--) {
		*dBaseline += *it; // baseline at start of event
		if (*itBasePkP < *it) itBasePkP = it;
		if (*itBasePkN > *it) itBasePkN = it;
		*dBasePost += *itt; // baseline at end of event
	}
	*dBaseline /= iBaselength;
	*dBasePost /= iBaselength;
	*dBasePeakN = (*dBaseline - *itBasePkN)*dScaleV;
	*dBasePeakP = (*itBasePkP - *dBaseline)*dScaleV;

	for (it = itBegin, itt = itEnd-1; it-itBegin < iBaselength; it++, itt--) { // RMS devations of baselines
		dTemp = *it - *dBaseline;
		*dBaseSigma += dTemp*dTemp;
		dTemp = *itt - *dBasePost;
		*dBasePostSigma += dTemp*dTemp;
	}
	*dBaseSigma = sqrt(*dBaseSigma/iBaselength);
	*dBasePostSigma = sqrt(*dBasePostSigma/iBaselength);

	*sPeaks = Peakfinder();
	if (*Peak.itPeak == 0) { // saturated event
		*bSaturated = true;
		for (it = Peak.itPeak; it < itEnd; it++) if (*it != 0) break;
		itSatEnd = it-1;
	} else {
		itSatEnd = Peak.itPeak;
	}
	if (iAverage) {
		Average();
		*sPeaks = Peakfinder();
	}
	if (*sPeaks == 0) return;

	for (it = Peak.itPeak; it > itBegin; it--) if ((*it < iThreshold) && (*(it-1) >= iThreshold)) break;
	*sTrigger = (it - itBegin)*dScaleT;
	for (it = Peak.itPeak; it < itEnd; it++) if (*it > *dBaseline-3*(*dBaseSigma)) break;
	*bFullWaveform = (it != itEnd);
	Peak.itEnd = (it == itEnd ? itEnd-1 : it);

	*sPeakX = (Peak.itPeak - itBegin)*dScaleT;
	*sPeakXs = (PeakS.itPeak - itBegin)*dScaleT;
	*sRise = (Peak.itPeak - Peak.itStart)*dScaleT;
	*sDecay = (Peak.itEnd - Peak.itPeak)*dScaleT;

	*dPeak0 = (*Peak.itStart - *Peak.itPeak)*dScaleV;
	*dPeak0s = (*PeakS.itStart - *PeakS.itPeak)*dScaleV;

	for (it = Peak.itStart; it <= Peak.itEnd; it++) *dIntegral += *it; // integrator
	*dIntegral = (*dIntegral)*2 - (*Peak.itStart + *Peak.itEnd);
	*dIntegral = ((*dBaseline)*(Peak.itEnd-Peak.itStart) - 0.5*(*dIntegral))*dScaleT*dScaleV;
}

inline void Event::Average() {
	auto itU = usTrace + (itSatEnd-itBegin);
	for (auto it = itSatEnd; it < itEnd; it++, itU++) {
		*it = 0.;
		for (auto itt = itU-iAverage; itt <= itU+iAverage; itt++) *it += *itt;
		*it /= (2.*iAverage+1.);
	}
}

int Event::Peakfinder() { // returns number of peaks found
	auto iPeakCutBL(8), iPeakCutMin(16); // Peaks less than this height don't get tagged, TODO add for DT5730
	auto iPrimaryTrigger(110); // primary peak in front of this point, TODO fix for other digitizers
	vector<Peak_t> vPeakCandidates, vFoundPeaks, vPrimaryPeaks;
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
		return 0;
	} else { // peaks found, select candidates of sufficient height
		auto itPrev = itBegin, itMin = itBegin;
		for (auto iter = vPeakCandidates.begin(); iter < vPeakCandidates.end(); iter++) { // evaluates heights of peak candidates
			itMin = (*iter).itPeak;
			for (it = (*iter).itPeak; it > itPrev; it--) { // looks for start of peak or minimum value
				if (*itMin < *it) itMin = it;
				if (*it > (*dBaseline)-3*(*dBaseSigma)) break;
			}
			if (it == itPrev) { // didn't reach baseline
				if (*itMin - *((*iter).itPeak) > iPeakCutMin) {
					(*iter).itStart = itMin;
					vFoundPeaks.push_back(*iter);
				}
			} else { // did reach baseline
				if ((*dBaseline) - *((*iter).itPeak) > iPeakCutBL) {
					(*iter).itStart = it;
					vFoundPeaks.push_back(*iter);
				}
			} // choosing it or itMin
		} // iter loop
	}

	if (vFoundPeaks.size() == 0) { // no peaks
		return 0;
	} else if (vFoundPeaks.size() == 1) { // only one peak, probably the majority of cases
		Peak = vFoundPeaks.front();
		return 1;
	} else { // two or more Peaks, primary is the tallest peak in the trigger region or the first
		for (auto iter = vFoundPeaks.begin(); iter < vFoundPeaks.end(); iter++) { // gathers Peaks in the trigger region (first 110ns of waveform)
			if ((*iter).itPeak-itBegin < iPrimaryTrigger) {
				vPrimaryPeaks.push_back(*iter);
			}
		}
		if (vPrimaryPeaks.size() == 0) { // no primaries, take first Peak
			Peak = vFoundPeaks.front();
		} else if (vPrimaryPeaks.size() == 1) { // one primary
			Peak = vPrimaryPeaks.front();
		} else { // multiple primaries, pick tallest
			Peak = vPrimaryPeaks.front();
			for (auto iter = vPrimaryPeaks.begin(); iter < vPrimaryPeaks.end(); iter++) {
				if (*iter > Peak) Peak = *iter;
			}
		}
		PeakS = itBegin;
		for (auto iter = vFoundPeaks.begin(); iter < vFoundPeaks.end(); iter++) {
			if ((*iter > PeakS) && !(*iter == Peak)) PeakS = *iter;
		}
	}
	return vFoundPeaks.size();
}

void Event::SetAddresses(vector<void*> add) {
	int i(0);
	bFullWaveform	= (bool*)add[i++];
	bSaturated		= (bool*)add[i++];

	sDecay			= (short*)add[i++];
	sRise			= (short*)add[i++];
	sPeakX			= (short*)add[i++];
	sPeakXs			= (short*)add[i++];
	sTrigger		= (short*)add[i++];
	sPeaks			= (short*)add[i++];

	dBaseline		= (double*)add[i++];
	dBaseSigma		= (double*)add[i++];
	dBasePost		= (double*)add[i++];
	dBasePostSigma	= (double*)add[i++];
	dPeak0			= (double*)add[i++];
	dPeak0s			= (double*)add[i++];
	dIntegral		= (double*)add[i++];
	dBasePeakP		= (double*)add[i++];
	dBasePeakN		= (double*)add[i++];
}
