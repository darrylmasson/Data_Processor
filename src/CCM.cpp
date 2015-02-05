#include "CCM.h"
#include <cstdlib>
#include <cmath>
#include <iostream>

float CCM::version = 2.8;
bool CCM::initialized = false;
shared_ptr<TTree> CCM::tree = nullptr;
int CCM::howmany = 0;

bool CCM::fullwave[8]			= {0,0,0,0,0,0,0,0};
bool CCM::saturated[8]			= {0,0,0,0,0,0,0,0};
bool CCM::truncated[8]			= {0,0,0,0,0,0,0,0};

short CCM::decay[8]				= {0,0,0,0,0,0,0,0};
short CCM::fstop[8]				= {0,0,0,0,0,0,0,0};
short CCM::rise[8]				= {0,0,0,0,0,0,0,0};
short CCM::peakx[8]				= {0,0,0,0,0,0,0,0};
short CCM::sstop[8]				= {0,0,0,0,0,0,0,0};

double CCM::baseline[8]			= {0,0,0,0,0,0,0,0};
double CCM::baseSigma[8]		= {0,0,0,0,0,0,0,0};
double CCM::basePost[8]			= {0,0,0,0,0,0,0,0};
double CCM::basePostSigma[8]	= {0,0,0,0,0,0,0,0};
double CCM::basePeakp[8]		= {0,0,0,0,0,0,0,0};
double CCM::basePeakn[8]		= {0,0,0,0,0,0,0,0};
double CCM::fullint[8]			= {0,0,0,0,0,0,0,0};
double CCM::fastint[8]			= {0,0,0,0,0,0,0,0};
double CCM::slowint[8]			= {0,0,0,0,0,0,0,0};
double CCM::peak0[8]			= {0,0,0,0,0,0,0,0};
double CCM::peak1[8]			= {0,0,0,0,0,0,0,0};
double CCM::peak2[8]			= {0,0,0,0,0,0,0,0};
double CCM::peakp[8]			= {0,0,0,0,0,0,0,0};

double CCM::gradient[8]			= {0,0,0,0,0,0,0,0};

CCM::CCM(const int ch, const int fast, const int slow, const int samples, const shared_ptr<Digitizer> digitizer) {
	failed = 0;
	eventlength = Event::Length();
	bool special = (digitizer->Special() == 0);
	fastTime = special ? fast/2 : fast;
	slowTime = special ? slow/2 : slow;
	gradSamples = special ? samples/2 : samples;
	id = ch;
	if ((id >= MAX_CH) || (id < 0)) failed |= method_error;
	CCM::howmany++;
	scaleV = digitizer->ScaleV();
	scaleT = digitizer->ScaleT();
}

CCM::~CCM() {
	std::cout << " CCM " << id << " d'tor ";
	CCM::howmany--;
//	CCM::tree = nullptr;
}

void CCM::root_init(shared_ptr<TTree> tree_in) {
	if (!CCM::initialized) {
		CCM::tree = tree_in;
		
		CCM::tree->Branch("FullWaveform",	CCM::fullwave,		"full_waveform[8]/O");
		CCM::tree->Branch("Saturated",		CCM::saturated,		"saturated[8]/O");
		CCM::tree->Branch("Truncated",		CCM::truncated,		"trunc[8]/O");
		
		CCM::tree->Branch("Risetime",		CCM::rise,			"risetime[8]/S");
		CCM::tree->Branch("Decaytime",		CCM::decay,			"decay_time[8]/S");
		CCM::tree->Branch("Faststop",		CCM::fstop,			"fstop[8]/S");
		CCM::tree->Branch("Slowstop",		CCM::sstop,			"sstop[8]/S");
		CCM::tree->Branch("Peakx",			CCM::peakx,			"peakx[8]/S");
		
		CCM::tree->Branch("Peakheight0",	CCM::peak0,			"pk_height0[8]/D");
		CCM::tree->Branch("Peakheight1",	CCM::peak1,			"pk_height1[8]/D");
		CCM::tree->Branch("Peakheight2",	CCM::peak2,			"pk_height1[8]/D");
		CCM::tree->Branch("PeakPos",		CCM::peakp,			"pk_pos[8]/D");
		CCM::tree->Branch("Baseline",		CCM::baseline,		"baseline[8]/D");
		CCM::tree->Branch("BaseSigma",		CCM::baseSigma,		"baseSigma[8]/D");
		CCM::tree->Branch("BasePost",		CCM::basePost,		"basePost[8]/D");
		CCM::tree->Branch("BasePostSigma",	CCM::basePostSigma,	"basePostSigma[8]/D");
		CCM::tree->Branch("BasePeakPos",	CCM::basePeakp,		"basePeakPos[8]/D");
		CCM::tree->Branch("BasePeakNeg",	CCM::basePeakn,		"basepeakNeg[8]/D");
		CCM::tree->Branch("Integral",		CCM::fullint,		"integrals[8]/D");
		CCM::tree->Branch("SlowInt",		CCM::slowint,		"slowintegral[8]/D");
		CCM::tree->Branch("FastInt",		CCM::fastint,		"fastintegral[8]/D");
		
		CCM::tree->Branch("Gradient",		CCM::gradient,		"gradient[8]/D");
		
		CCM::initialized = true;
	}
}

void CCM::evaluate(const shared_ptr<Event> event) {
	int start(0), stop(eventlength-1), i(0);
	int temp(0), fast(0), slow(0);
	long Tfast(0), Tslow(0), Tfull(0);
	double threshold(event->baseline - 3*event->baseSigma);
	
	// normalizing baseline values
	CCM::baseline[id] = (event->baseline - event->zero)*scaleV;
	CCM::baseSigma[id] = event->baseSigma*scaleV;
	CCM::basePeakp[id] = (event->b_pk_p - event->baseline)*scaleV;
	CCM::basePeakn[id] = (event->baseline - event->b_pk_n)*scaleV;
	CCM::peakp[id] = (event->peak_pos - event->baseline)*scaleV;
	CCM::basePost[id] = (event->basePost - event->zero)*scaleV;
	CCM::basePostSigma[id] = event->basePostSigma*scaleV;

	for (i = 0; i < eventlength; i++) { // determine integration bounds
		if ((stop == eventlength-1) && (event->peak_x + i < eventlength) && (event->trace[event->peak_x + i] > threshold)) stop = event->peak_x + i;
		if ((start == 0) && (event->peak_x - i > -1) && (event->trace[event->peak_x - i] > threshold)) start = event->peak_x - i;
		if ((start != 0) && (stop != eventlength-1)) break;
	}
	
	// boolean results
	CCM::fullwave[id] = (stop != eventlength-1);
	CCM::saturated[id] = (event->peak_y == 0);
	CCM::truncated[id] = (start + slowTime >= eventlength);
	
	fast = min(fastTime, eventlength - start);
	slow = min(slowTime, eventlength - start);

	CCM::peak0[id] = (event->baseline - event->peak_y) * scaleV;

	if ((event->peak_x + 2 < eventlength) && (event->peak_x > 1) && !(CCM::saturated[id])) { // peak averaging
		for (i = -1; i < 2; i++) temp += event->trace[event->peak_x + i];
		CCM::peak1[id] = (event->baseline - 0.333*temp)*scaleV;
		temp += (event->trace[event->peak_x - 2] + event->trace[event->peak_x + 2]);
		CCM::peak2[id] = (event->baseline - 0.2*temp)*scaleV;
	} else CCM::peak2[id] = CCM::peak1[id] = CCM::peak0[id];
	
	for (i = start; i < eventlength; i++) { // integrator
		if (i <= stop) Tfull += event->trace[i];
		if (i <= start + fast) Tfast += event->trace[i];
		if (i <= start + slow) Tslow += event->trace[i];
		if ((i > stop) && (i > start + slow)) break;
	}
	Tfull <<= 1;
	Tfast <<= 1;
	Tslow <<= 1;
	Tfull -= (event->trace[start] + event->trace[stop]);
	Tfast -= (event->trace[start] + event->trace[start + fast]);
	Tslow -= (event->trace[start] + event->trace[start + slow]);
	
	CCM::rise[id] = (event->peak_x - start)*scaleT;
	CCM::decay[id] = (stop - event->peak_x)*scaleT;

	CCM::peakx[id] = event->peak_x*scaleT;
	CCM::fstop[id] = (start + fast)*scaleT;
	CCM::sstop[id] = (start + slow)*scaleT;
	
	CCM::fullint[id] = (event->baseline * (stop - start) - Tfull/2.) * scaleV * scaleT;
	CCM::slowint[id] = (event->baseline * (slow) - Tslow/2.) * scaleV * scaleT;
	CCM::fastint[id] = (event->baseline * (fast) - Tfast/2.) * scaleV * scaleT;
	
	if (event->peak_x + gradSamples + 1 < eventlength) {
		temp = 0;
		for (i = -1; i < 2; i++) temp += event->trace[event->peak_x + gradSamples + i]; // average with adjacent points to reduce statistical fluctuations
		temp /= 3;
		CCM::gradient[id] = (gradSamples * (event->baseline - event->peak_y) == 0) ? -1 : (temp - event->peak_y)/(double)(gradSamples * (event->baseline - event->peak_y));
	} else CCM::gradient[id] = -1;
	
	return;
}
