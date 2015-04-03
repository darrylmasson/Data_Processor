#ifndef EVENT_AVE_H
#define EVENT_AVE_H

#ifndef EVENT_H
#include "Event.h"
#endif

class Event_ave : public Event {
	private:
		int iAverage;
		double dPeakY; // these need to be floating point, not integers
		double dBasePkP; // but still technically ADC counts
		double dBasePkN; // everything else inherited
		double dPeakPos;
		double dScale; // 1/average
		std::unique_ptr<double[]> dTrace;

	public:
		Event_ave(int len, shared_ptr<Digitizer>);
		~Event_ave();
		virtual void SetAverage(int average);
		virtual void Analyze();
		virtual double Trace(int i)		{return ((i < 0) ? dBaseline : (i < ciSamples ? dTrace[i] : dBasePost));}
		virtual double Peak_y()			{return dPeakY;}
		virtual double B_pk_p()			{return dBasePkP;}
		virtual double B_pk_n()			{return dBasePkN;}
		virtual double Peak_pos()		{return dPeakPos;}
};

#endif // EVENT_AVE_H
