#include "Event.h"
#include <cmath>
#include <algorithm>
#include <iostream>

int Event::siLength = 0;
int Event::siHowMany = 0;

Event::Event(int len, std::shared_ptr<Digitizer> dig, int dc_offset, int threshold_in) {
	Event::siHowMany++;
	if (Event::siHowMany == 1) Event::siLength = len;
	digitizer = dig;
	iSpecial = digitizer->Special();
	iBaselength = digitizer->Baselength();
	usTrace = nullptr;
	iThreshold = threshold_in;
	dZero = digitizer->Resolution()*(1. - (double)dc_offset/65535.); // conversion from wavedump documentation
	iFailed = 0;
	if ((iSpecial == 0) && (Event::siHowMany == 1)) Event::siLength >>= 1;
	iEventlength = Event::siLength;
}

Event::~Event() {
	if (g_verbose) std::cout << " event " << --Event::siHowMany << " d'tor "; // no ID tag on Event classes
	usTrace = nullptr;
	digitizer.reset();
}

void Event::Set(unsigned short* in) {
	usTrace = in;
	int i(0);
	if (iSpecial > 0) for (i = 0; i < iEventlength; i++) usTrace[i] >>= iSpecial; // special resolution
	if (iSpecial == 0) for (i = 0; i < iEventlength; i++) usTrace[i] = (usTrace[2*i] + usTrace[2*i+1]) >> 1; // special samplerate
	dBaseline = 0;
	dBaseSigma = 0;
	usPeakY = -1;
	usPeakX = 0;
	usBasePkP = 0;
	usBasePkN = -1;
	usPeakPos = 0;
	usTrigger = 0;
	dBasePost = 0;
	dBasePostSigma = 0;
	double dTemp(0);
	for (i = 0; i < iEventlength; i++) {
		usPeakPos = std::max(usPeakPos, usTrace[i]);
		if (i < iBaselength) {
			dBaseline += usTrace[i]; // baseline at start of event
			usBasePkP = std::max(usTrace[i],usBasePkP);
			usBasePkN = std::min(usTrace[i],usBasePkN);
			dBasePost += usTrace[iEventlength-iBaselength+i]; // baseline at end of event
		}
		if (usPeakY > usTrace[i]) { // finding primary peak
			usPeakY = usTrace[i];
			usPeakX = i;
		}
		if ((usTrigger == 0) && (usTrace[i] < iThreshold)) usTrigger = i; // finding trigger
	}
	dBaseline /= iBaselength;
	dBasePost /= iBaselength;
	for (i = 0; i < iBaselength; i++) { // RMS devations of baselines
		dTemp = usTrace[i] - dBaseline;
		dBaseSigma += dTemp*dTemp;
		dTemp = usTrace[iEventlength-iBaselength+i] - dBasePost;
		dBasePostSigma += dTemp*dTemp;
	}
	dBaseSigma = sqrt(dBaseSigma/iBaselength);
	dBasePostSigma = sqrt(dBasePostSigma/iBaselength);
}