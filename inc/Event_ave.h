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
};

#endif // EVENT_AVE_H
