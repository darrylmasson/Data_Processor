#include "LAPL.h"
#include <cmath>

std::string LAPL::version = "1.0";
bool LAPL::initialized = false;

LAPL::LAPL(int eventlength) {
	results = new lapl_results_t;
	failed = 0;
	int s(0), t(0);
	if (!initialized) {
		initialized = true;
		for (s = 0; s < length; s++) {
			ex[s] = new (std::nothrow) double[eventlength];
			if (ex[s] == NULL) failed |= 1;
			for (t = 0; t < eventlength; t++) ex[s][t] = exp(-1.*s*t);
	}	}
	for (s = 0; s < length; s++) e[s] = ex[s];
}

void LAPL::evaluate(const Event* event, const Digitizer* digitizer) {
	static int s(0), t(0);
	for (s = 0; s < length; s++) {
		results->lapl[s] = 0;
		for (t = 0; t < event->length; t++) results->lapl[s] += e[s][t]*event->trace[t];
	}
}

LAPL::~LAPL() {
	if (initialized) {
		delete[] *ex;
		initialized = false;
	}
	delete results;
}
