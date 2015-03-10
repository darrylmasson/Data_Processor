#include "Event.h"
#include <cmath>
#include <algorithm>
#include <iostream>

int Event::length = 0;
int Event::howmany = 0;

Event::Event(int len, std::shared_ptr<Digitizer> dig, int dc_offset, int threshold_in) {
	Event::howmany++;
	if (Event::howmany == 1) Event::length = len;
	digitizer = dig;
	special = digitizer->Special();
	baselength = digitizer->Baselength();
	us_trace = nullptr;
	threshold = threshold_in;
	zero = digitizer->Resolution()*(1. - (double)dc_offset/65535.); // conversion from wavedump documentation
	failed = 0;
	if ((special == 0) && (Event::howmany == 1)) Event::length >>= 1;
	eventlength = Event::length;
}

Event::~Event() {
	if (g_verbose) std::cout << " event " << --Event::howmany << " d'tor ";
	us_trace = nullptr;
	digitizer.reset();
}

void Event::Set(unsigned short* in) {
	us_trace = in;
	int i(0);
	if (special > 0) for (i = 0; i < eventlength; i++) us_trace[i] >>= special; // special resolution
	if (special == 0) for (i = 0; i < eventlength; i++) us_trace[i] = (us_trace[2*i] + us_trace[2*i+1]) >> 1; // special samplerate
	d_baseline = 0;
	d_baseSigma = 0;
	us_peak_y = -1;
	us_peak_x = 0;
	d_b_pk_p = 0;
	d_b_pk_n = -1;
	us_peak_pos = 0;
	us_trigger = 0;
	d_basePost = 0;
	d_basePostSigma = 0;
	double d_temp(0);
	for (i = 0; i < eventlength; i++) {
		us_peak_pos = std::max(us_peak_pos, us_trace[i]);
		if (i < baselength) {
			d_baseline += us_trace[i];
			us_b_pk_p = std::max(us_trace[i],us_b_pk_p);
			us_b_pk_n = std::min(us_trace[i],us_b_pk_n);
			d_basePost += us_trace[eventlength-baselength+i];
		}
		if (us_peak_y > us_trace[i]) {
			us_peak_y = us_trace[i];
			us_peak_x = i;
		}
		if ((us_trigger == 0) && (us_trace[i] < threshold)) us_trigger = i;
	}
	d_baseline /= baselength;
	d_basePost /= baselength;
	for (i = 0; i < baselength; i++) {
		d_temp = us_trace[i] - d_baseline;
		d_baseSigma += d_temp*d_temp;
		d_temp = us_trace[eventlength-baselength+i] - d_basePost;
		d_basePostSigma += d_temp*d_temp;
	}
	d_baseSigma = sqrt(d_baseSigma/baselength);
	d_basePostSigma = sqrt(d_basePostSigma/baselength);
}