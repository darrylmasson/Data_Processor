#ifndef METHOD_H // ABC for processing codes
#define METHOD_H // also some fundamental types

#ifndef EVENT_H
#include "Event.h"
#endif

#ifndef EVENT_AVE_H
#include "Event_ave.h"
#endif

#include <string>
#include "TFile.h"
#include "TTree.h"

class Method {
	protected:
		int failed;
		int id; // which entry in the array to use
		int eventlength;
		double d_scaleT; // scale factors from Digitizer class
		double d_scaleV;
	
	public:
		Method() {};
		virtual ~Method() {};
		virtual void evaluate(const shared_ptr<Event>) = 0;
		int Failed() {return failed;}
		int GetID() {return id;}
};

#endif // METHOD_H
