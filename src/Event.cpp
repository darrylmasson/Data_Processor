#include "Event.h"
#include <cmath>
#include <algorithm>
#include <iostream>

int Event::siLength = 0;
int Event::siHowMany = 0;

Event::Event() {
	Event::siHowMany++;
	if (g_verbose) cout << "Event c'tor\n";
}

Event::Event(int len, shared_ptr<Digitizer> digitizer) : ciSamples(len) {
	if (g_verbose) cout << "Event " << Event::siHowMany << " c'tor\n";
	Event::siHowMany++;
	iEventlength = ciSamples;
	iSpecial = digitizer->Special();
	iBaselength = digitizer->Baselength();
	dBaseScale = 1./iBaselength;
	uspTrace = nullptr;
	iFailed = 0;
	if (iSpecial == 0) iEventlength >>= 1;
	Event::siLength = iEventlength;
}

Event::~Event() {
	Event::siHowMany--;
	if (g_verbose) cout << "Event " << Event::siHowMany << " d'tor\n"; // no ID tag on Event classes
	uspTrace = nullptr;
}

void Event::Analyze() {
	int i(0);
	if (iSpecial > 0) for (i = 0; i < iEventlength; i++) uspTrace[i] >>= iSpecial; // special resolution
	if (iSpecial == 0) for (i = 0; i < iEventlength; i++) uspTrace[i] = (uspTrace[2*i] + uspTrace[2*i+1]) >> 1; // special samplerate
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
		usPeakPos = max(usPeakPos, uspTrace[i]);
		if (i < iBaselength) {
			dBaseline += uspTrace[i]; // baseline at start of event
			usBasePkP = max(uspTrace[i],usBasePkP);
			usBasePkN = min(uspTrace[i],usBasePkN);
			dBasePost += uspTrace[iEventlength-iBaselength+i]; // baseline at end of event
		}
		if (usPeakY > uspTrace[i]) { // finding primary peak
			usPeakY = uspTrace[i];
			usPeakX = i;
		}
		if ((usTrigger == 0) && (uspTrace[i] < iThreshold)) usTrigger = i; // finding trigger
	}
	dBaseline *= dBaseScale;
	dBasePost *= dBaseScale;
	for (i = 0; i < iBaselength; i++) { // RMS devations of baselines
		dTemp = uspTrace[i] - dBaseline;
		dBaseSigma += dTemp*dTemp;
		dTemp = uspTrace[iEventlength-iBaselength+i] - dBasePost;
		dBasePostSigma += dTemp*dTemp;
	}
	dBaseSigma = sqrt(dBaseSigma*dBaseScale);
	dBasePostSigma = sqrt(dBasePostSigma*dBaseScale);
}