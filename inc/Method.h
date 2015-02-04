#ifndef METHOD_H // ABC for processing codes
#define METHOD_H // also some fundamental types

#ifndef EVENT_H
#include "Event.h"
#endif
#include <string>
#include "TFile.h"
#include "TTree.h"

using namespace std;

const int MAX_CH = 8;
const char path[] = {"/data/NeutronGenerator/"};
const int sizeof_f_header = 86; // sizeof(f_header_t) returns the wrong value
const int sizeof_ev_header = 16;

class Method {
	protected:
		int failed;
		int id;
		int eventlength;
	
	public:
		Method() {};
		virtual ~Method() {};
		virtual void evaluate(const shared_ptr<Event>) =0;
		int Failed() {return failed;}
		int GetID() {return id;}
};

enum methods_t {
	CCM_t = 0,
	DFT_t,
	XSQ_t,
	
	NUM_METHODS // add others before this
};

const char method_names[NUM_METHODS][12] = {
	"CCM_PGA",
	"FOURIER",
	"CHISQUARED"
};

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

enum particle {n = 0, y, P};

struct config_t {
	unsigned short	mask;
	int				nchan;
	int				chan[MAX_CH]; // only the first nchan entries used
	int				numEvents;
	int				eventsize; // bytes, incl event header
	int				eventlength; // sampes
	unsigned int	dc_offset[MAX_CH];
	unsigned int	threshold[MAX_CH];
	float			gain[MAX_CH][P]; // for fitter
	int				trig_post;
	int				pga_samples[MAX_CH];
	int				fastTime[MAX_CH];
	int				slowTime[MAX_CH];
	bool			method_active[NUM_METHODS]; // active for this run
	bool			already_done; // for timestamps
	bool			method_done[NUM_METHODS]; // previously processed
};

struct thread_data_t {
	unsigned short*		data;
	shared_ptr<Event>	event;
	unique_ptr<Method>	methods[NUM_METHODS];
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

#endif // METHOD_H
