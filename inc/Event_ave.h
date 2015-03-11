#ifndef EVENT_AVE_H
#define EVENT_AVE_H

#ifndef EVENT_H
#include "Event.h"
#endif

class Event_ave : public Event {
	private:
		int average;
		double d_peak_y; // these need to be floating point, not integers
		double d_b_pk_p; // but still technically ADC counts
		double d_b_pk_n; // everything else inherited
		double d_peak_pos;
		std::unique_ptr<double[]> d_trace;

	public:
		Event_ave(int len, std::shared_ptr<Digitizer> dig, int dc_offset, int threshold_in, int average_in);
		~Event_ave();
		virtual void Set(unsigned short* in);
	//	virtual unsigned short& Trigger()	{return us_trigger;}
		virtual double Trace(int i)			{return d_trace[i];}
	//	virtual unsigned short& Peak_x()	{return us_peak_x;}
		virtual double Peak_y()				{return d_peak_y;}
		virtual double B_pk_p()				{return d_b_pk_p;}
		virtual double B_pk_n()				{return d_b_pk_n;}
		virtual double Peak_pos()			{return d_peak_pos;}
	//	virtual auto Zero()					{return d_zero;}
	//	virtual auto Baseline()				{return d_baseline;}
	//	virtual auto BaseSigma()			{return d_baseSigma;}
	//	virtual auto BasePost()				{return d_basePost;}
	//	virtual auto BasePostSigma()		{return d_basePostSigma;}
};

#endif // EVENT_AVE_H
