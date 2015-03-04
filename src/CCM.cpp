#include "CCM.h"
#include <cstdlib>
#include <cmath>
#include <iostream>

float CCM::sf_version = 2.81;
bool CCM::sb_initialized = false;
unique_ptr<TTree> CCM::tree = nullptr;
int CCM::si_howmany = 0;

bool CCM::sb_fullwave[8]		= {0,0,0,0,0,0,0,0};
bool CCM::sb_saturated[8]		= {0,0,0,0,0,0,0,0};
bool CCM::sb_truncated[8]		= {0,0,0,0,0,0,0,0};

short CCM::ss_decay[8]			= {0,0,0,0,0,0,0,0};
short CCM::ss_fstop[8]			= {0,0,0,0,0,0,0,0};
short CCM::ss_rise[8]			= {0,0,0,0,0,0,0,0};
short CCM::ss_peakx[8]			= {0,0,0,0,0,0,0,0};
short CCM::ss_sstop[8]			= {0,0,0,0,0,0,0,0};

double CCM::sd_baseline[8]		= {0,0,0,0,0,0,0,0};
double CCM::sd_baseSigma[8]		= {0,0,0,0,0,0,0,0};
double CCM::sd_basePost[8]		= {0,0,0,0,0,0,0,0};
double CCM::sd_basePostSigma[8]	= {0,0,0,0,0,0,0,0};
double CCM::sd_basePeakp[8]		= {0,0,0,0,0,0,0,0};
double CCM::sd_basePeakn[8]		= {0,0,0,0,0,0,0,0};
double CCM::sd_fullint[8]		= {0,0,0,0,0,0,0,0};
double CCM::sd_fastint[8]		= {0,0,0,0,0,0,0,0};
double CCM::sd_slowint[8]		= {0,0,0,0,0,0,0,0};
double CCM::sd_peak0[8]			= {0,0,0,0,0,0,0,0};
double CCM::sd_peak1[8]			= {0,0,0,0,0,0,0,0};
double CCM::sd_peak2[8]			= {0,0,0,0,0,0,0,0};
double CCM::sd_peakp[8]			= {0,0,0,0,0,0,0,0};

double CCM::sd_gradient[8]		= {0,0,0,0,0,0,0,0};

CCM::CCM(const int ch, const int fast, const int slow, const int samples, const shared_ptr<Digitizer> digitizer) {
	failed = 0;
	eventlength = Event::Length();
	bool special = (digitizer->Special() == 0);
	i_fastTime = (special ? fast/2 : fast);
	i_slowTime = (special ? slow/2 : slow);
	i_gradSamples = (special ? samples/2 : samples);
	id = ch;
	if ((id >= MAX_CH) || (id < 0)) failed |= method_error;
	CCM::si_howmany++;
	d_scaleV = digitizer->ScaleV();
	d_scaleT = digitizer->ScaleT();
}

CCM::~CCM() {
	if (g_verbose) cout << " CCM " << id << " d'tor ";
	CCM::si_howmany--;
}

void CCM::root_init(TTree* tree_in) {
	if (!CCM::sb_initialized) {
		CCM::tree = unique_ptr<TTree>(tree_in);
		
		CCM::tree->Branch("FullWaveform",	CCM::sb_fullwave,		"full_waveform[8]/O");
		CCM::tree->Branch("Saturated",		CCM::sb_saturated,		"saturated[8]/O");
		CCM::tree->Branch("Truncated",		CCM::sb_truncated,		"trunc[8]/O");
		
		CCM::tree->Branch("Risetime",		CCM::ss_rise,			"risetime[8]/S");
		CCM::tree->Branch("Decaytime",		CCM::ss_decay,			"decay_time[8]/S");
		CCM::tree->Branch("Faststop",		CCM::ss_fstop,			"fstop[8]/S");
		CCM::tree->Branch("Slowstop",		CCM::ss_sstop,			"sstop[8]/S");
		CCM::tree->Branch("Peakx",			CCM::ss_peakx,			"peakx[8]/S");
		
		CCM::tree->Branch("Peakheight0",	CCM::sd_peak0,			"pk_height0[8]/D");
		CCM::tree->Branch("Peakheight1",	CCM::sd_peak1,			"pk_height1[8]/D");
		CCM::tree->Branch("Peakheight2",	CCM::sd_peak2,			"pk_height1[8]/D");
		CCM::tree->Branch("PeakPos",		CCM::sd_peakp,			"pk_pos[8]/D");
		CCM::tree->Branch("Baseline",		CCM::sd_baseline,		"baseline[8]/D");
		CCM::tree->Branch("BaseSigma",		CCM::sd_baseSigma,		"baseSigma[8]/D");
		CCM::tree->Branch("BasePost",		CCM::sd_basePost,		"basePost[8]/D");
		CCM::tree->Branch("BasePostSigma",	CCM::sd_basePostSigma,	"basePostSigma[8]/D");
		CCM::tree->Branch("BasePeakPos",	CCM::sd_basePeakp,		"basePeakPos[8]/D");
		CCM::tree->Branch("BasePeakNeg",	CCM::sd_basePeakn,		"basepeakNeg[8]/D");
		CCM::tree->Branch("Integral",		CCM::sd_fullint,		"integrals[8]/D");
		CCM::tree->Branch("SlowInt",		CCM::sd_slowint,		"slowintegral[8]/D");
		CCM::tree->Branch("FastInt",		CCM::sd_fastint,		"fastintegral[8]/D");
		
		CCM::tree->Branch("Gradient",		CCM::sd_gradient,		"gradient[8]/D");
		
		CCM::sb_initialized = true;
	}
}

void CCM::evaluate(const shared_ptr<Event> event) {
	int i_start(0), i_stop(eventlength-1), i(0);
	int i_temp(0), i_fast(0), i_slow(0);
	long l_fastint(0), l_slowint(0), l_fullint(0);
	double d_threshold(event->baseline - 3*event->baseSigma);
	
	// normalizing baseline values
	CCM::sd_baseline[id] = (event->baseline - event->zero)*d_scaleV;
	CCM::sd_baseSigma[id] = event->baseSigma*d_scaleV;
	CCM::sd_basePeakp[id] = (event->b_pk_p - event->baseline)*d_scaleV;
	CCM::sd_basePeakn[id] = (event->baseline - event->b_pk_n)*d_scaleV;
	CCM::sd_peakp[id] = (event->peak_pos - event->baseline)*d_scaleV;
	CCM::sd_basePost[id] = (event->basePost - event->zero)*d_scaleV;
	CCM::sd_basePostSigma[id] = event->basePostSigma*d_scaleV;

	for (i = 0; i < eventlength; i++) { // determine integration bounds
		if ((i_stop == (eventlength-1)) && ((event->peak_x + i) < eventlength) && (event->trace[event->peak_x + i] > d_threshold)) i_stop = (event->peak_x + i);
		if ((i_start == 0) && ((event->peak_x - i) > -1) && (event->trace[event->peak_x - i] > d_threshold)) i_start = (event->peak_x - i);
		if ((i_start != 0) && (i_stop != (eventlength-1))) break;
	}
	
	// boolean results
	CCM::sb_fullwave[id] = (i_stop != (eventlength-1));
	CCM::sb_saturated[id] = (event->peak_y == 0);
	CCM::sb_truncated[id] = ((i_start + i_slowTime) >= eventlength);
	
	i_fast = min(i_fastTime, eventlength -1 - i_start);
	i_slow = min(i_slowTime, eventlength -1 - i_start);

	CCM::sd_peak0[id] = (event->baseline - event->peak_y) * d_scaleV;

	if (((event->peak_x + 2) < eventlength) && (event->peak_x > 1) && !(CCM::sb_saturated[id])) { // peak averaging
		for (i = -1; i < 2; i++) i_temp += event->trace[event->peak_x + i];
		CCM::sd_peak1[id] = (event->baseline - (0.333*i_temp))*d_scaleV;
		i_temp += (event->trace[event->peak_x - 2] + event->trace[event->peak_x + 2]);
		CCM::sd_peak2[id] = (event->baseline - (0.2*i_temp))*d_scaleV;
	} else CCM::sd_peak2[id] = CCM::sd_peak1[id] = CCM::sd_peak0[id];
	
	for (i = i_start; i < eventlength; i++) { // integrator
		if (i <= i_stop) l_fullint += event->trace[i];
		if (i <= (i_start + i_fast)) l_fastint += event->trace[i];
		if (i <= (i_start + i_slow)) l_slowint += event->trace[i];
		if ((i > i_stop) && (i > (i_start + i_slow))) break;
	}
	l_fullint <<= 1;
	l_fastint <<= 1;
	l_slowint <<= 1;
	l_fullint -= (event->trace[i_start] + event->trace[i_stop]);
	l_fastint -= (event->trace[i_start] + event->trace[i_start + i_fast]);
	l_slowint -= (event->trace[i_start] + event->trace[i_start + i_slow]);
	
	CCM::ss_rise[id] = (event->peak_x - i_start)*d_scaleT;
	CCM::ss_decay[id] = (i_stop - event->peak_x)*d_scaleT;

	CCM::ss_peakx[id] = event->peak_x*d_scaleT;
	CCM::ss_fstop[id] = (i_start + i_fast)*d_scaleT;
	CCM::ss_sstop[id] = (i_start + i_slow)*d_scaleT;
	
	CCM::sd_fullint[id] = ((event->baseline * (i_stop - i_start)) - (0.5*l_fullint)) * d_scaleV * d_scaleT;
	CCM::sd_slowint[id] = ((event->baseline * (i_slow)) - (0.5*l_slowint)) * d_scaleV * d_scaleT;
	CCM::sd_fastint[id] = ((event->baseline * (i_fast)) - (0.5*l_fastint)) * d_scaleV * d_scaleT;
	
	if ((event->peak_x + i_gradSamples + 1) < eventlength) {
		i_temp = 0;
		for (i = -1; i < 2; i++) i_temp += event->trace[event->peak_x + i_gradSamples + i]; // average with adjacent points to reduce statistical fluctuations
		i_temp /= 3;
		CCM::sd_gradient[id] = (i_gradSamples * (event->baseline - event->peak_y) == 0) ? -1 : (i_temp - event->peak_y)/(double)(i_gradSamples * (event->baseline - event->peak_y));
	} else CCM::sd_gradient[id] = -1;
	
	return;
}
