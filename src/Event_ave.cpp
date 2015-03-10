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
	try { trace = std::unique_ptr<double[]>(new double[eventlength]);}
	catch (std::bad_alloc& ba) {failed |= alloc_error; return;}
}

Event_ave::~Event_ave() {
	if (g_verbose) std::cout << " event_ave " << Event::howmany << " d'tor ";
	trace.reset();
	digitizer.reset();
}

void Event_ave::Set(unsigned short* in) {
	int i(0), j(0);
	if (special > 0) for (i = 0; i < eventlength; i++) d_trace[i] = in[i] >> special; // special resolution
	if (special == 0) for (i = 0; i < eventlength; i++) d_trace[i] = (in[2*i] + in[2*i+1]) >> 1; // special samplerate
	if (average > 0) {
		for (i = 0; i < eventlength; i++) {
			d_trace[i] = 0;
			for (j = 0; j < average; j++) d_trace[i] += in[i+j];
			d_trace[i] *= 1./average;
	}	}
	d_baseline = 0;
	d_baseSigma = 0;
	d_peak_y = 15000; // some arbitratily high number as -1 doesn't work for floats
	us_peak_x = 0;
	d_b_pk_p = 0;
	d_b_pk_n = 15000;
	d_peak_pos = 0;
	us_trigger = 0;
	d_basePost = 0;
	d_basePostSigma = 0;
	double d_temp(0);
	for (i = 0; i < eventlength; i++) {
		d_peak_pos = std::max(d_peak_pos, d_trace[i]);
		if (i < baselength) {
			d_baseline += d_trace[i];
			d_b_pk_p = std::max(d_trace[i],d_b_pk_p);
			d_b_pk_n = std::min(d_trace[i],d_b_pk_n);
			d_basePost += d_trace[eventlength-baselength+i];
		}
		if (d_peak_y > d_trace[i]) {
			d_peak_y = d_trace[i];
			us_peak_x = i;
		}
		if ((us_trigger == 0) && (d_trace[i] < threshold)) us_trigger = i;
	}
	d_baseline /= baselength;
	d_basePost /= baselength;
	for (i = 0; i < baselength; i++) {
		d_temp = d_trace[i] - d_baseline;
		d_baseSigma += d_temp*d_temp;
		d_temp = d_trace[eventlength-baselength+i] - d_basePost;
		d_basePostSigma += d_temp*d_temp;
	}
	d_baseSigma = sqrt(d_baseSigma/baselength);
	d_basePostSigma = sqrt(d_basePostSigma/baselength);
}