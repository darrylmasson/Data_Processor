#include "Event_ave.h"
#include <cmath>
#include <algorithm>
#include <iostream>

Event_ave::Event_ave(int len, std::shared_ptr<Digitizer> dig, int dc_offset, int threshold_in, int average_in) {
	if (Event::howmany == 1) Event::length = len;
	digitizer = dig;
	special = digitizer->Special();
	baselength = digitizer->Baselength();
	trace = nullptr;
	threshold = threshold_in;
	average = average_in;
	zero = digitizer->Resolution()*(1. - (double)dc_offset/65535.); // conversion from wavedump documentation
	failed = 0;
	if ((special == 0) && (Event::howmany == 1)) Event::length >>= 1;
	if ((average > 0) && (Event::howmany == 1)) Event::length -= average;
	eventlength = Event::length;
	try { trace.reset(new double[eventlength]);}
	catch (bad_alloc& ba) {failed |= alloc_error; return;}
}

Event_ave::~Event_ave() {
	if (g_verbose) std::cout << " event_ave " << Event::howmany << " d'tor ";
	trace.reset();
	digitizer.reset();
}

void Event_ave::Set(unsigned short* in) {
	int i(0), j(0);
	if (special > 0) for (i = 0; i < eventlength; i++) trace[i] = in[i] >> special; // special resolution
	if (special == 0) for (i = 0; i < eventlength; i++) trace[i] = (in[2*i] + in[2*i+1]) >> 1; // special samplerate
	if (average > 0) {
		for (i = 0; i < eventlength; i++) {
			trace[i] = 0;
			for (j = 0; j < average; j++) trace[i] += in[i+j];
			trace[i] *= 1./average;
	}	}
	baseline = 0;
	baseSigma = 0;
	peak_y = 15000; // some arbitratily high number as -1 doesn't work for floats
	peak_x = 0;
	b_pk_p = 0;
	b_pk_n = 15000;
	peak_pos = 0;
	trigger = 0;
	basePost = 0;
	basePostSigma = 0;
	double temp(0);
	for (i = 0; i < eventlength; i++) {
		peak_pos = std::max(peak_pos, trace[i]);
		if (i < baselength) {
			baseline += trace[i];
			b_pk_p = std::max(trace[i],b_pk_p);
			b_pk_n = std::min(trace[i],b_pk_n);
			basePost += trace[eventlength-baselength+i];
		}
		if (peak_y > trace[i]) {
			peak_y = trace[i];
			peak_x = i;
		}
		if ((trigger == 0) && (trace[i] < threshold)) trigger = i;
	}
	baseline /= baselength;
	basePost /= baselength;
	for (i = 0; i < baselength; i++) {
		temp = trace[i] - baseline;
		baseSigma += temp*temp;
		temp = trace[eventlength-baselength+i] - basePost;
		basePostSigma += temp*temp;
	}
	baseSigma = sqrt(baseSigma/baselength);
	basePostSigma = sqrt(basePostSigma/baselength);
}