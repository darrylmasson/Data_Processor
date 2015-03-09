#ifndef EVENT_H
#define EVENT_H

#include <memory>
#ifndef DIGITIZER_H
#include "Digitizer.h"
#endif

class Event {
	protected:
		int failed;
		int threshold;
		int eventlength;
		int baselength;
		int special;
		static int length;
		static int howmany;

	public:
		Event() {++Event::howmany;} // for Event_ave
		Event(int len, std::shared_ptr<Digitizer> dig, int dc_offset, int threshold_in);
		~Event();
		virtual void Set(unsigned short* in);
		int Failed() {return failed;}
		static const int& Length() {return Event::length;}
		std::shared_ptr<Digitizer> digitizer;
		unsigned short trigger;
		virtual unsigned short* trace;
		unsigned short peak_x;
		virtual unsigned short peak_y; // everything measured in ADC counts
		virtual unsigned short b_pk_p;
		virtual unsigned short b_pk_n;
		virtual unsigned short peak_pos;
		double zero;
		double baseline;
		double baseSigma;
		double basePost;
		double basePostSigma;
};

#endif // EVENT_H
