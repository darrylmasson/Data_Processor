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
		int iFailed;
		int id; // which entry in the array to use
		int iEventlength;
		double dScaleT; // scale factors from Digitizer class
		double dScaleV;
	
	public:
		Method() {};
		Method(int ch, int length, shared_ptr<Digitizer> digitizer) : id(ch), iEventlength(length), dScaleT(digitizer->ScaleT()), dScaleV(digitizer->ScaleV()), iFailed(0) {}
		virtual ~Method() {};
		virtual void Analyze(const shared_ptr<Event>) = 0;
		virtual void SetParameters(void* val, int which, shared_ptr<Digitizer> digitizer) = 0;
		int Failed() {return iFailed;}
		int GetID() {return id;}
};

#endif // METHOD_H
