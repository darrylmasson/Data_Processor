#include "Event_ave.h"
#include <cmath>
#include <algorithm>
#include <iostream>

Event_ave::Event_ave(int len, std::shared_ptr<Digitizer> dig, int dc_offset, int threshold_in, int average_in) {
	if (Event::siHowMany == 1) Event::siLength = len;
	digitizer = dig;
	iSpecial = digitizer->Special();
	iBaselength = digitizer->Baselength();
	dTrace = nullptr;
	threshold = threshold_in;
	iAverage = average_in;
	dScale = 1./iAverage;
	dZero = digitizer->Resolution()*(1. - (double)dc_offset/65535.); // conversion from wavedump documentation
	iFailed = 0;
	if ((iSpecial == 0) && (Event::siHowMany == 1)) Event::siLength >>= 1;
	if ((iAverage > 0) && (Event::siHowMany == 1)) Event::siLength -= iAverage;
	iEventlength = Event::siLength;
	try { dTrace = std::unique_ptr<double[]>(new double[iEventlength]);}
	catch (std::bad_alloc& ba) {failed |= alloc_error; return;}
}

Event_ave::~Event_ave() {
	if (g_verbose) std::cout << " event_ave " << Event::siHowMany << " d'tor "; // d'tor for Event automatically called, so no decrement
	dTrace.reset();
	digitizer.reset();
}

void Event_ave::Set(unsigned short* in) {
	int i(0), j(0);
	if (iSpecial > 0) for (i = 0; i < iEventlength; i++) dTrace[i] = in[i] >> iSpecial; // special resolution
	if (iSpecial == 0) for (i = 0; i < iEventlength; i++) dTrace[i] = (in[2*i] + in[2*i+1]) >> 1; // special samplerate
	for (i = 0; i < iEventlength; i++) { // waveform averaging
		dTrace[i] = 0;
		for (j = 0; j < iAverage; j++) dTrace[i] += in[i+j];
		dTrace[i] *= dScale;
	}
	dBaseline = 0;
	dBaseSigma = 0;
	dPeakY = 15000; // some arbitratily high number as -1 doesn't work for floats
	usPeakX = 0;
	dBasePkP = 0;
	dBasePkN = 15000;
	dPeakPos = 0;
	usTrigger = 0;
	dBasePost = 0;
	dBasePostSigma = 0;
	double dTemp(0);
	for (i = 0; i < iEventlength; i++) {
		dPeakPos = std::max(dPeakPos, dTrace[i]);
		if (i < iBaselength) {
			dBaseline += dTrace[i]; // baseline at start
			dBasePkP = std::max(dTrace[i],dBasePkP);
			dBasePkN = std::min(dTrace[i],dBasePkN);
			dBasePost += dTrace[iEventlength-iBaselength+i]; // baseline at end
		}
		if (dPeakY > dTrace[i]) { // peakfinder
			dPeakY = dTrace[i];
			usPeakX = i;
		}
		if ((usTrigger == 0) && (dTrace[i] < threshold)) usTrigger = i; // trigger
	}
	dBaseline /= iBaselength;
	dBasePost /= iBaselength;
	for (i = 0; i < iBaselength; i++) { // RMS baseline deviation
		dTemp = dTrace[i] - dBaseline;
		dBaseSigma += dTemp*dTemp;
		dTemp = dTrace[iEventlength-iBaselength+i] - dBasePost;
		dBasePostSigma += dTemp*dTemp;
	}
	dBaseSigma = sqrt(dBaseSigma/iBaselength);
	dBasePostSigma = sqrt(dBasePostSigma/iBaselength);
}