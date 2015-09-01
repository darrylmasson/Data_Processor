#include "Event.h"

int Event::siLength = 0;
int Event::siHowMany = 0;

Event::Event() : ciSamples(0) {
	Event::siHowMany++;
	if (g_verbose) cout << "Event c'tor\n";
}

Event::Event(int eventlength, int baselength, int average, unsigned short* start, unsigned short* end, int threshold, int chan) : ciChan(chan) {
	if (g_verbose) cout << "Event " << chan << " c'tor\n";
	iFailed = 0;
	iEventlength = eventlength;
	iBaselength = baselength;
	iAverage = average;
	try{
		uspTrace.assign(start,end);
		dTrace.assign(uspTrace.size()-2*iAverage);
	} catch (bad_alloc& ba) {
		iFailed |= (1 << alloc_error);
		return;
	}
	itBegin = vTrace.begin();
	itEnd = vTrace.end();
	iThreshold = threshold;
}

Event::~Event() {
	if (g_verbose) cout << "Event " << ciChan << " d'tor\n";
}

void Event::Analyze() {
	PreAnalyze();
	dBaseline = 0.;
	dBaseSigma = 0.;
	dBasePost = 0.;
	dBasePostSigma = 0.;
	dIntegral = 0.;
	bSaturated = false;
	bPileUp = false;
	bFullWaveform = true;
	itPeakY = itBasePkP = itBasePkP = itPeakPos = itTrigger = itSatEnd = itPulseStart = itPulseEnd = itBegin;
	auto it = itBegin, itt = itBegin, itMin = itBegin;
	auto dTemp(0.);
	vector<vector<double>::iterator> vPeakCandidates(10), vAllTriggers(10);
	vector<double>::iterator this_peak;
	for (it, itt = itEnd-1; it < itEnd; it++, itt--) {
		if (*itPeakPos < *it) itPeakPos = it;
		if (it-itBegin < iBaselength) {
			dBaseline += *it; // baseline at start of event
			if (*itBasePkP < *it) itBasePkP = it;
			if (*itBasePkN > *it) itBasePkN = it;
			dBasePost += *itt; // baseline at end of event
		}
		if ((it-itBegin > 0) && (*it < iThreshold) && (*(it-1) >= iThreshold)) vAllTriggers.push_back(it);
	}
	dBaseline /= iBaselength;
	dBasePost /= iBaselength;

	for (it = itBegin, itt = itEnd-1; it-itBegin < iBaselength; it++, itt--) { // RMS devations of baselines
		dTemp = *it - dBaseline;
		dBaseSigma += dTemp*dTemp;
		dTemp = *itt - dBasePost;
		dBasePostSigma += dTemp*dTemp;
	}
	dBaseSigma = sqrt(dBaseSigma/iBaselength);
	dBasePostSigma = sqrt(dBasePostSigma/iBaselength);

	for (it = itBegin; it < itEnd; it += iBaselength>>1) { // peak-finding logic
		this_peak = it;
		for (itt = it; itt < it + iBaselength; itt++) {
			if (itt == itEnd) break;
			if (*this_peak > *itt) this_peak = itt;
		}
		if (this_peak - it < iBaselength>>2) continue; // in the first quarter of the window, falling edge or nothing.
		else if (this_peak - it < 3*iBaselength/4) { // not in the end of the window
			if (*this_peak < dBaseline-3*dBaseSigma) vPeakCandidates.push_back(this_peak); // cuts noise
		} else { // in last quarter of window
			continue;
		}
	}

	if (vPeakCandidates.size() == 0) { // no peaks found, use default values
		return;
	} else {
		auto iter = vPeakCandidates.begin();
		itMin = *iter;
		for (it = *iter; it >= itBegin; it--) {
			if (*itMin < *it) itMin = it;

		}
	}

	if (*itPeakY == 0) { // saturated event
		bSaturated = true;
		for (it = itPeakY; it < itEnd; it++) if (*it != 0) break;
		itSatEnd = it-1;
	}

	for (it = itPulseStart; it <= itPulseEnd; it++) dIntegral += *it; // integrator
	dIntegral = dIntegral*2 - (*itPulseStart + *itPulseEnd);
	dIntegral = dBaseline*(itPulseEnd-itPulseStart) - 0.5*dIntegral; // at this point, bins*clock cycles, not V*ns
}

inline void Event::PreAnalyze() {
	for (auto itU = uspTrace.begin(), auto itD = itBegin; itU != uspTrace.end(); itU++, itD++) *itD = *itU;
	// insert averaging and other shenanigans here
}