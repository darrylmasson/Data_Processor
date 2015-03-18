#ifndef NGDP_TYPES_H
#define NGDP_TYPES_H

#include <memory>

// note on variable names:
// most variables have their type as a one or two character prefix
// using normal names, ie u = unsigned, s = static, d = double, etc

using namespace std;

const int MAX_CH = 8;
const char path[] = {"/data/NeutronGenerator"};
const int sizeof_f_header = 86; // sizeof(f_header_t) returns the wrong value
const int sizeof_ev_header = 16;

extern bool g_verbose; // for debugging purposes

class Digitizer;
class Event; // forward declarations
class Event_ave;
class Method;
class CCM;
class DFT;
class XSQ;
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

enum particle {n = 0, y, P}; // used exclusively by XSQ

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
	XSQ_t,
	LAP_t,
	
	NUM_METHODS // add others before this
};

const char method_names[NUM_METHODS][12] = {
	"CCM_PGA",
	"FOURIER",
	"CHISQUARED",
	"LAPLACE"
};

const char error_message[err_dummy_last][64] = {
	"No error",
	"Error in method constructor",
	"Invalid digitizer",
	"Root error",
	"File not found",
	"Allocation error",
	"Thread error",
	"Bad value in file"
};

struct config_t {
	bool			method_active[NUM_METHODS]; // active for this run
	bool			already_done; // for timestamps
	bool			method_done[NUM_METHODS]; // previously processed
	unsigned short	mask;
	int				nchan;
	int				numEvents;
	int				eventsize; // bytes, incl event header
	int				eventlength; // samples
	int				trig_post;
	int				average;
	int				chan[MAX_CH]; // only the first nchan entries used
	int				pga_samples[MAX_CH];
	int				fastTime[MAX_CH];
	int				slowTime[MAX_CH];
	unsigned int	dc_offset[MAX_CH];
	unsigned int	threshold[MAX_CH];
	float			gain[MAX_CH][P]; // for fitter
};

struct thread_data_t {
	unsigned short*		data;
	shared_ptr<Event>	event;
	shared_ptr<Method>	methods[NUM_METHODS];
	const bool*			activated;
};

struct f_header_t {
	char dig_name[12];
	unsigned short mask;
	unsigned int ev_len;
	int trig_post;
	unsigned int dc_off[8];
	unsigned int threshold[8];
};

#endif // NGDP_TYPES_H