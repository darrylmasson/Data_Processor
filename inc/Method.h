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
		shared_ptr<Event> event;
	
	public:
		Method(int ch, int length, shared_ptr<Digitizer> digitizer) : iFailed(0), id(ch), iEventlength(length), dScaleT(digitizer->ScaleT()), dScaleV(digitizer->ScaleV()) {}
		virtual ~Method() {event.reset();}
		virtual void Analyze() = 0;
		virtual void SetEvent(shared_ptr<Event> ev) {event = ev;}
		virtual void SetParameters(void* val, int which, shared_ptr<Digitizer> digitizer) = 0;
		int Failed() {return iFailed;}
		int GetID() {return id;}
};

#endif // METHOD_H
