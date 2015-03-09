#ifndef EVENT_AVE_H
#define EVENT_AVE_H

#ifndef EVENT_H
#include "Event.h"
#endif

class Event_ave : public Event {
	private:
		int average;

	public:
		Event_ave(int len, std::shared_ptr<Digitizer> dig, int dc_offset, int threshold_in, int average_in);
		~Event_ave();
		virtual void Set(unsigned short* in);
		unique_ptr<double[]> trace;
		double peak_y; // these need to be floating point, not integers
		double b_pk_p;
		double b_pk_n; // everything else inherited
		double peak_pos;
};

#endif // EVENT_AVE_H
