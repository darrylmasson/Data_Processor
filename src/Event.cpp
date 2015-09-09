#include "Event.h"

float Event::sfVersion = 1.0;

Event::Event() : ciChan(-1) {
	if (g_verbose) cout << "Event c'tor\n";
}

Event::Event(int eventlength, int baselength, int average, unsigned short* start, int threshold, int chan) : ciChan(chan) {
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
	uspTrace = start;
	try {
		usvTrace.assign(start,start + eventlength);
		vTrace.assign((int)usvTrace.size()-2*iAverage, 0);
	} catch (bad_alloc& ba) {
		cout << error_message[alloc_error] << "Event\n";
		iFailed = 1;
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
	*dBaseline = 0.;
	*dBaseSigma = 0.;
	*dBasePost = 0.;
	*dBasePostSigma = 0.;
	*dIntegral = 0.;
	*dPeak0 = 0.;
	*dBasePeakP = 0.;
	*dBasePeakN = 0.;

	*sDecay = 0;
	*sRise = 0;
	*sPeakX = 0;
	*sTrigger = 0;

	*bSaturated = false;
	*bPileUp = false;
	*bFullWaveform = true;

	itBasePkP = itBasePkN = itSatEnd = Peak.itPeak = Peak.itStart = itBegin;
	Peak.itEnd = itEnd;
	auto it = itBegin, itt = itBegin;
	auto dTemp(0.);
	auto iPeakCut(8); // Peaks less than this height don't get tagged
	vector<Peak_t> vPeakCandidates, vFoundPeaks, vPrimaryPeaks;

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
	
	cout << *dBaseline << " ";

	for (it = itBegin, itt = itEnd-1; it-itBegin < iBaselength; it++, itt--) { // RMS devations of baselines
		dTemp = *it - *dBaseline;
		*dBaseSigma += dTemp*dTemp;
		dTemp = *itt - *dBasePost;
		*dBasePostSigma += dTemp*dTemp;
	}
	*dBaseSigma = sqrt(*dBaseSigma/iBaselength);
	*dBasePostSigma = sqrt(*dBasePostSigma/iBaselength);

	for (it = itBegin; it < itEnd; it += iBaselength/2) { // peak-finding logic
		Peak.itPeak = it;
		for (itt = it; itt < it + iBaselength; itt++) { // finds minimum value in window
			if (itt == itEnd) break;
//			cout << *itt << " " << *Peak.itPeak << " ";
			if (*Peak.itPeak > *itt) {
				Peak.itPeak = itt;// cout << Peak.itPeak - itBegin << " ";
			}
		}
		if (Peak.itPeak - it < iBaselength/4) {
//			cout << "92 ";
			continue; // in the first quarter of the window, falling edge or nothing.
		}
		else if (Peak.itPeak - it < 3*iBaselength/4) { // not in the end of the window
			if (*(Peak.itPeak) < iThreshold) { cout << "96 ";
				vPeakCandidates.push_back(Peak); // cuts noise
			} else {
				cout << "99 ";			
				continue;
			}
		} else { // in last quarter of window
			cout << "103 ";
			continue;
		} cout << '\n';
	}
//	cout << vPeakCandidates.size() << " " << vFoundPeaks.size() << " " << vPrimaryPeaks.size() << " ";
	if (vPeakCandidates.size() == 0) { // no Peaks found, use default values
		return;
	} else { // peaks found, select candidates of sufficient height
		auto itPrev = itBegin, itMin = itBegin;
		for (auto iter = vPeakCandidates.begin(); iter < vPeakCandidates.end(); iter++) { // evaluates heights of peak candidates
			itMin = (*iter).itPeak; cout << "105 "; cout << itBegin - (*iter).itPeak << " ";
			for (it = (*iter).itPeak; it > itPrev; it--) { // looks for start of peak or minimum value
				if (*itMin < *it) itMin = it;
				if (*it > (*dBaseline)-3*(*dBaseSigma)) break;
				cout << it-(*iter).itPeak << " ";
			} cout << "109 "; cout << *itMin << " " << *dBaseline << *((*iter).itPeak) << " ";
			if (it == itPrev) { // didn't reach baseline
				cout << "111 ";
				if (*itMin - *((*iter).itPeak) > iPeakCut) {
					(*iter).itStart = itMin;
					vFoundPeaks.push_back(*iter);
				} cout << "115 ";
			} else { // did reach baseline
				cout << "117 ";
				if ((*dBaseline) - *((*iter).itPeak) > iPeakCut) { cout << "118 ";
					(*iter).itStart = it;
					vFoundPeaks.push_back(*iter);
				} cout << "121 ";
			} // choosing it or itMin
			itPrev = (*iter).itPeak; cout << "123 ";
		} // iter loop
	}
//	cout << "Location ";
	if (vFoundPeaks.size() == 0) { // no peaks
		return;
	} else if (vFoundPeaks.size() == 1) { // only one peak, probably the majority of cases
		Peak = vFoundPeaks.front();
	} else { // two or more Peaks, primary is the tallest peak in the trigger region or the first
		*bPileUp = true;
		for (auto iter = vFoundPeaks.begin(); iter < vFoundPeaks.end(); iter++) { // gathers Peaks in the trigger region (first third of waveform)
			if ((*iter).itPeak-itBegin < iEventlength/3) {
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
				if (*((*iter).itPeak) - *((*iter).itStart) < *(Peak.itPeak) - *(Peak.itStart)) Peak = *iter;
			}
		}
	}
//	cout << "Basics ";
	for (it = Peak.itPeak; it > itBegin; it--) if ((*it < iThreshold) && (*(it-1) >= iThreshold)) break;
	*sTrigger = (it - itBegin)*dScaleT;
	for (it = Peak.itPeak; it < itEnd; it++) if (*it > *dBaseline-3*(*dBaseSigma)) break;
	*bFullWaveform = (it != itEnd);
	Peak.itEnd = (it == itEnd ? itEnd-1 : it);

	*sPeakX = (Peak.itPeak - itBegin)*dScaleT;
	*sRise = (Peak.itPeak - Peak.itStart)*dScaleT;
	*sDecay = (Peak.itEnd - Peak.itPeak)*dScaleT;

	if (*Peak.itPeak == 0) { // saturated event
		*bSaturated = true;
		for (it = Peak.itPeak; it < Peak.itEnd; it++) if (*it != 0) break;
		itSatEnd = it-1;
	} else {
		itSatEnd = Peak.itPeak;
	}

	*dPeak0 = (*dBaseline - *Peak.itPeak)*dScaleV;
//	cout << "Integral\n";
	for (it = Peak.itStart; it <= Peak.itEnd; it++) *dIntegral += *it; // integrator
	*dIntegral = (*dIntegral)*2 - (*Peak.itStart + *Peak.itEnd);
	*dIntegral = ((*dBaseline)*(Peak.itEnd-Peak.itStart) - 0.5*(*dIntegral))*dScaleInt;
}

inline void Event::PreAnalyze() {
	auto itD = itBegin;
	usvTrace.assign(uspTrace,uspTrace+iEventlength);
	for (auto itU = usvTrace.begin(); itU < usvTrace.end() && itD < itEnd; itU++, itD++) {
		*itD = *itU;
	//	cout << *itD << " " << *itU << " ";
	}
	//cout << '\n';
	// insert averaging and other shenanigans here
}

void Event::SetAddresses(vector<void*> add) {
	int i(0);
	bFullWaveform	= (bool*)add[i++];
	bSaturated		= (bool*)add[i++];
	bPileUp			= (bool*)add[i++];

	sDecay			= (short*)add[i++];
	sRise			= (short*)add[i++];
	sPeakX			= (short*)add[i++];
	sTrigger		= (short*)add[i++];

	dBaseline		= (double*)add[i++];
	dBaseSigma		= (double*)add[i++];
	dBasePost		= (double*)add[i++];
	dBasePostSigma	= (double*)add[i++];
	dPeak0			= (double*)add[i++];
	dIntegral		= (double*)add[i++];
	dBasePeakP		= (double*)add[i++];
	dBasePeakN		= (double*)add[i++];
}
