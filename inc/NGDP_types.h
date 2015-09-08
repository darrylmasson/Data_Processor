#ifndef NGDP_TYPES_H
#define NGDP_TYPES_H

#include <memory>
#include <string>
#include <cstring>
#include <cstdint>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>

// note on variable names:
// most variables have their type as a one or two character prefix
// using normal names, ie u = unsigned, s = static/short, d = double, etc

using namespace std;

const int MAX_CH = 8;
const string sWorkingDir(WORKING_DIR, string(WORKING_DIR).find_last_of('/'));
const int sizeof_f_header = 86; // sizeof(f_header_t) returns the wrong value due to padding
const int sizeof_ev_header = 16;

extern bool g_verbose; // for debugging purposes

class Digitizer;
class Event; // forward declarations
class Event_ave;
class Method;
class CCM;
class DFT;
class XSQ_TF1;
class XSQ_NEW;
class LAP;

enum ret_codes_t {
	no_error = 0,
	method_error,
	dig_error,
	root_error,
	file_error,
	alloc_error,
	thread_error,
	config_file_error,

	err_dummy_last
};

enum dig_id_t { // easy way for internal identification of digitizer
	dt5751 = 0,
	dt5751des,
	dt5730,
	v1724,
	invalid_dig
};

enum methods_t {
	CCM_t = 0,
	DFT_t,
	XSQ_TF1_t,
	LAP_t,
	XSQ_NEW_t,

	NUM_METHODS // add others before this
};

enum discrims_t {
	d_CCM_t = 0,
	d_DFT_t,
	d_NGM_t,
	d_LAP_t,
	d_PGA_t,
	d_WBS_t,

	NUM_DISCRIMS
}; // discrete discrimination methods

#endif // NGDP_TYPES_H