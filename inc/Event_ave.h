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
	//	virtual auto Trigger()		->	decltype(us_trigger)		{return us_trigger;}
		virtual auto Trace(int i)	->	decltype(d_trace[0])		{return d_trace[i];}
	//	virtual auto Peak_x()		->	decltype(us_peak_x)			{return us_peak_x;}
		virtual auto Peak_y()		->	decltype(d_peak_y)			{return d_peak_y;}
		virtual auto B_pk_p()		->	decltype(d_b_pk_p			{return d_b_pk_p;}
		virtual auto B_pk_n()		->	decltype(d_b_pk_n			{return d_b_pk_n;}
		virtual auto Peak_pos()		->	decltype(d_peak_pos)		{return d_peak_pos;}
	//	virtual auto Zero()			->	decltype(d_zero)			{return d_zero;}
	//	virtual auto Baseline()		->	decltype(d_baseline)		{return d_baseline;}
	//	virtual auto BaseSigma()	->	decltype(d_baseSigma)		{return d_baseSigma;}
	//	virtual auto BasePost()		->	decltype(d_basePost)		{return d_basePost;}
	//	virtual auto BasePostSigma()->	decltype(d_basePostSigma)	{return d_basePostSigma;}
};

#endif // EVENT_AVE_H
