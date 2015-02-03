#include "Event.h"
#include <cmath>
#include <algorithm>

int Event::length = 0;

Event::Event(int len, std::shared_ptr<Digitizer> dig, int dc_offset, int threshold_in) {
	if (Event::length == 0) Event::length = len;
	digitizer = dig;
	special = digitizer->Special();
	baselength = digitizer->Baselength();
	trace = nullptr;
	threshold = threshold_in;
	zero = digitizer->Resolution()*(1. - (double)dc_offset/65535.);
	failed = 0;
	if ((special == 0) && (Event::length == len)) Event::length >>= 1;
	eventlength = Event::length;
}

Event::~Event() {
	trace = nullptr;
	digitizer.reset();
}

void Event::Set(unsigned short* in) {
	trace = in;
	int i(0);
	if (special > 0) for (i = 0; i < eventlength; i++) trace[i] >>= special;
	if (special == 0) for (i = 0; i < eventlength; i++) trace[i] = (trace[2*i] + trace[2*i+1]) >> 1;
	baseline = 0;
	baseSigma = 0;
	peak_y = -1;
	peak_x = 0;
	b_pk_p = 0;
	b_pk_n = -1;
	peak_pos = 0;
	trigger = 0;
	basePost = 0;
	basePostSigma = 0;
	sat_end = 0;
	double temp(0);
	for (i = 0; i < eventlength; i++) {
		peak_pos = std::max(peak_pos, trace[i]);
		if (i < baselength) {
			baseline += trace[i];
			b_pk_p = std::max(trace[i],b_pk_p);
			b_pk_n = std::min(trace[i],b_pk_n);
			basePost += trace[eventlength-1-i];
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
		temp = trace[eventlength-1-i] - basePost;
		basePostSigma += temp*temp;
	}
	baseSigma = sqrt(baseSigma/baselength);
	basePostSigma = sqrt(basePostSigma/baselength);
	if (peak_y == 0) {
		sat_end = peak_x;
		while (((sat_end+1) < eventlength) && (trace[sat_end+1] == 0)) sat_end++;
	}
}