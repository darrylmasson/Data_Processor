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
		Event(int len, std::shared_ptr<Digitizer> dig, int dc_offset, int threshold_in);
		~Event();
		void Set(unsigned short* in);
		int Failed() {return failed;}
		static const int& Length() {return Event::length;}
		std::shared_ptr<Digitizer> digitizer;
		unsigned short trigger;
		unsigned short* trace;
		unsigned short peak_x;
		unsigned short peak_y; // everything measured in ADC counts
		unsigned short b_pk_p;
		unsigned short b_pk_n;
		unsigned short peak_pos;
		double zero;
		double baseline;
		double baseSigma;
		double basePost;
		double basePostSigma;
};

#endif // EVENT_H
