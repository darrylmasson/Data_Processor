#ifndef LAPL_H // laplace transform
#define LAPL_H

#include "Process.h"

struct lapl_results_t {
	double lapl[256];
};

class LAPL : public Process {
	private:
		lapl_results_t* results;
		const int length = 256;
		const double* e[256];
	
	public:
		LAPL(int);
		~LAPL();
		virtual void* GetResults() {return (void*)results;}
		virtual codes_t GetType() {return LAPL_t;}
		virtual void evaluate(const Event*, const Digitizer*);

		static std::string version;
		static bool initialized;
		static double* ex[256];
};


#endif // LAPL_H
