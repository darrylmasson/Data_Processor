#include "CCM.h"
#include <cstdlib>
#include <cmath>
#include <iostream>

float CCM::sfVersion = 2.81;
bool CCM::sbInitialized = false;
unique_ptr<TTree> CCM::tree = nullptr;
int CCM::siHowMany = 0;

bool CCM::sbFullWave[8]			= {0,0,0,0,0,0,0,0};
bool CCM::sbSaturated[8]		= {0,0,0,0,0,0,0,0}; // up to 8 channels on digitizers
bool CCM::sbTruncated[8]		= {0,0,0,0,0,0,0,0};

short CCM::ssDecay[8]			= {0,0,0,0,0,0,0,0};
short CCM::ssFastStop[8]		= {0,0,0,0,0,0,0,0};
short CCM::ssRise[8]			= {0,0,0,0,0,0,0,0};
short CCM::ssPeakX[8]			= {0,0,0,0,0,0,0,0};
short CCM::ssSlowStop[8]		= {0,0,0,0,0,0,0,0};

double CCM::sdBaseline[8]		= {0,0,0,0,0,0,0,0};
double CCM::sdBaseSigma[8]		= {0,0,0,0,0,0,0,0};
double CCM::sdBasePost[8]		= {0,0,0,0,0,0,0,0};
double CCM::sdBasePostSigma[8]	= {0,0,0,0,0,0,0,0};
double CCM::sdBasePeakP[8]		= {0,0,0,0,0,0,0,0};
double CCM::sdBasePeakN[8]		= {0,0,0,0,0,0,0,0};
double CCM::sdFullInt[8]		= {0,0,0,0,0,0,0,0};
double CCM::sdFastInt[8]		= {0,0,0,0,0,0,0,0};
double CCM::sdSlowInt[8]		= {0,0,0,0,0,0,0,0};
double CCM::sdPeak0[8]			= {0,0,0,0,0,0,0,0};
double CCM::sdPeak1[8]			= {0,0,0,0,0,0,0,0};
double CCM::sdPeak2[8]			= {0,0,0,0,0,0,0,0};
double CCM::sdPeakP[8]			= {0,0,0,0,0,0,0,0};

double CCM::sdGradient[8]		= {0,0,0,0,0,0,0,0};

CCM::CCM() {
	CCM::siHowMany++;
	if (g_verbose) cout << "CCM c'tor\n";
}

CCM::CCM(int ch, int length, shared_ptr<Digitizer> digitizer) : Method(ch, length, digitizer) {
	CCM::siHowMany++;
	if (g_verbose) cout << "CCM " << id << " c'tor\n";
	if ((id >= MAX_CH) || (id < 0)) iFailed |= method_error;
}

CCM::~CCM() {
	if (g_verbose) cout << "CCM " << id << " d'tor\n";
	CCM::siHowMany--;
}

void CCM::SetParameters(void* val, int which, shared_ptr<Digitizer> digitizer) {
	auto special = digitizer->Special();
	int value = *((int*)val);
	switch (which) {
		case 0 :
			iFastTime = special ? value : value/2; // special == 0 means half values. Derp.
			if (g_verbose) cout << "Fast " << iFastTime << '\n';
			break;
		case 1 :
			iSlowTime = special ? value : value/2;
			if (g_verbose) cout << "Slow " << iSlowTime << '\n';
			break;
		case 2 :
			iPGASamples = special ? value : value/2;
			if (g_verbose) cout << "PGA " << iPGASamples << '\n';
			break;
		default : break;
	}
}

void CCM::root_init(TTree* tree_in) {
	if (!CCM::sbInitialized) {
		CCM::tree = unique_ptr<TTree>(tree_in);
		
		CCM::tree->Branch("FullWaveform",	CCM::sbFullWave,		"full_waveform[8]/O");
		CCM::tree->Branch("Saturated",		CCM::sbSaturated,		"saturated[8]/O");
		CCM::tree->Branch("Truncated",		CCM::sbTruncated,		"trunc[8]/O");
		
		CCM::tree->Branch("Risetime",		CCM::ssRise,			"risetime[8]/S");
		CCM::tree->Branch("Decaytime",		CCM::ssDecay,			"decay_time[8]/S");
		CCM::tree->Branch("Faststop",		CCM::ssFastStop,		"fstop[8]/S");
		CCM::tree->Branch("Slowstop",		CCM::ssSlowStop,		"sstop[8]/S");
		CCM::tree->Branch("Peakx",			CCM::ssPeakX,			"peakx[8]/S");
		
		CCM::tree->Branch("Peakheight0",	CCM::sdPeak0,			"pk_height0[8]/D");
		CCM::tree->Branch("Peakheight1",	CCM::sdPeak1,			"pk_height1[8]/D");
		CCM::tree->Branch("Peakheight2",	CCM::sdPeak2,			"pk_height1[8]/D");
		CCM::tree->Branch("PeakPos",		CCM::sdPeakP,			"pk_pos[8]/D");
		CCM::tree->Branch("Baseline",		CCM::sdBaseline,		"baseline[8]/D");
		CCM::tree->Branch("BaseSigma",		CCM::sdBaseSigma,		"baseSigma[8]/D");
		CCM::tree->Branch("BasePost",		CCM::sdBasePost,		"basePost[8]/D");
		CCM::tree->Branch("BasePostSigma",	CCM::sdBasePostSigma,	"basePostSigma[8]/D");
		CCM::tree->Branch("BasePeakPos",	CCM::sdBasePeakP,		"basePeakPos[8]/D");
		CCM::tree->Branch("BasePeakNeg",	CCM::sdBasePeakN,		"basepeakNeg[8]/D");
		CCM::tree->Branch("Integral",		CCM::sdFullInt,			"integrals[8]/D");
		CCM::tree->Branch("SlowInt",		CCM::sdSlowInt,			"slowintegral[8]/D");
		CCM::tree->Branch("FastInt",		CCM::sdFastInt,			"fastintegral[8]/D");
		
		CCM::tree->Branch("Gradient",		CCM::sdGradient,		"gradient[8]/D");
		
		CCM::sbInitialized = true;
	}
}

void CCM::Analyze() {
	auto iStart(0), iStop(iEventlength-1), i(0), iFast(0), iSlow(0);
	auto lFastint(0l), lSlowint(0l), lFullint(0l);
	auto dThreshold(event->Baseline() - 3*event->BaseSigma()), dTemp(0.);
	
	// normalizing baseline values
	CCM::sdBaseline[id] = (event->Baseline() - event->Zero())*dScaleV;
	CCM::sdBaseSigma[id] = event->BaseSigma()*dScaleV;
	CCM::sdBasePeakP[id] = (event->B_pk_p() - event->Baseline())*dScaleV;
	CCM::sdBasePeakN[id] = (event->Baseline() - event->B_pk_n())*dScaleV;
	CCM::sdPeakP[id] = (event->Peak_pos() - event->Baseline())*dScaleV;
	CCM::sdBasePost[id] = (event->BasePost() - event->Zero())*dScaleV;
	CCM::sdBasePostSigma[id] = event->BasePostSigma()*dScaleV;

	for (i = 0; i < iEventlength; i++) { // determine integration bounds
		if ((iStop == (iEventlength-1)) && ((event->Peak_x() + i) < iEventlength) && (event->Trace(event->Peak_x() + i) > dThreshold)) iStop = (event->Peak_x() + i);
		if ((iStart == 0) && ((event->Peak_x() - i) > -1) && (event->Trace(event->Peak_x() - i) > dThreshold)) iStart = (event->Peak_x() - i);
		// first checks to see if start/stop has been found, then if the array index is valid, then checks the value
		if ((iStart != 0) && (iStop != (iEventlength-1))) break;
	}
	
	// boolean results
	CCM::sbFullWave[id] = (iStop != (iEventlength-1));
	CCM::sbSaturated[id] = (event->Peak_y() == 0);
	CCM::sbTruncated[id] = ((iStart + iSlowTime) >= iEventlength);
	
	iFast = min(iFastTime, iEventlength -1 - iStart);
	iSlow = min(iSlowTime, iEventlength -1 - iStart); // local integration limits for fast and slow

	CCM::sdPeak0[id] = (event->Baseline() - event->Peak_y()) * dScaleV;

	if (((event->Peak_x() + 2) < iEventlength) && (event->Peak_x() > 1) && !(CCM::sbSaturated[id])) { // peak averaging
		for (i = -1; i < 2; i++) dTemp += event->Trace(event->Peak_x() + i);
		CCM::sdPeak1[id] = (event->Baseline() - (0.333*dTemp))*dScaleV; // averaged with adjacent samples
		dTemp += (event->Trace(event->Peak_x() - 2) + event->Trace(event->Peak_x() + 2));
		CCM::sdPeak2[id] = (event->Baseline() - (0.2*dTemp))*dScaleV; // averaged with two adjacent samples
	} else CCM::sdPeak2[id] = CCM::sdPeak1[id] = CCM::sdPeak0[id];
	
	for (i = iStart; i < iEventlength; i++) { // integrator
		if (i <= iStop) lFullint += event->Trace(i);
		if (i <= (iStart + iFast)) lFastint += event->Trace(i);
		if (i <= (iStart + iSlow)) lSlowint += event->Trace(i);
		if ((i > iStop) && (i > (iStart + iSlow))) break;
	}
	lFullint <<= 1;
	lFastint <<= 1; // trapezoid rule: integral = f(0) + 2*f(1) + ... + 2*f(n-1) + f(n)
	lSlowint <<= 1; // faster to do sum f(i), double, and subtract endpoints
	lFullint -= (event->Trace(iStart) + event->Trace(iStop));
	lFastint -= (event->Trace(iStart) + event->Trace(iStart + iFast));
	lSlowint -= (event->Trace(iStart) + event->Trace(iStart + iSlow));
	
	CCM::ssRise[id] = (event->Peak_x() - iStart)*dScaleT;
	CCM::ssDecay[id] = (iStop - event->Peak_x())*dScaleT;

	CCM::ssPeakX[id] = event->Peak_x()*dScaleT;
	CCM::ssFastStop[id] = (iStart + iFast)*dScaleT;
	CCM::ssSlowStop[id] = (iStart + iSlow)*dScaleT;
	
	CCM::sdFullInt[id] = ((event->Baseline() * (iStop - iStart)) - (0.5*lFullint)) * dScaleV * dScaleT;
	CCM::sdSlowInt[id] = ((event->Baseline() * (iSlow)) - (0.5*lSlowint)) * dScaleV * dScaleT; // baseline subtraction
	CCM::sdFastInt[id] = ((event->Baseline() * (iFast)) - (0.5*lFastint)) * dScaleV * dScaleT;
	
	if ((event->Peak_x() + iPGASamples + 1) < iEventlength) { // PGA
		dTemp = 0;
		for (i = -1; i < 2; i++) dTemp += event->Trace(event->Peak_x() + iPGASamples + i); // average with adjacent points to reduce statistical fluctuations
		dTemp *= 0.333;
		CCM::sdGradient[id] = (iPGASamples * (event->Baseline() - event->Peak_y()) == 0) ? -1 : (dTemp - event->Peak_y())/(double)(iPGASamples * (event->Baseline() - event->Peak_y()));
	} else CCM::sdGradient[id] = -1;
	
	return;
}
