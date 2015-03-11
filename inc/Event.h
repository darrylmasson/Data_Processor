#ifndef EVENT_H
#define EVENT_H

#include <memory>
#ifndef DIGITIZER_H
#include "Digitizer.h"
#endif

class Event {
	protected: // these get inherited
		int failed;
		int threshold;
		int eventlength;
		int baselength;
		int special;
		static int length;
		static int howmany;
		
		unsigned short us_peak_y; // everything measured in ADC counts
		unsigned short us_b_pk_p;
		unsigned short us_b_pk_n;
		unsigned short us_peak_pos;
		unsigned short* us_trace;
		unsigned short us_trigger;
		unsigned short us_peak_x;
		double d_zero;
		double d_baseline;
		double d_baseSigma;
		double d_basePost;
		double d_basePostSigma;

	public:
		Event() {++Event::howmany;} // for Event_ave
		Event(int len, std::shared_ptr<Digitizer> dig, int dc_offset, int threshold_in);
		~Event();
		virtual void Set(unsigned short* in);
		int Failed() {return failed;}
		static const int& Length() {return Event::length;}
		std::shared_ptr<Digitizer> digitizer;
		virtual auto& Trigger()		->	decltype(us_trigger)		{return us_trigger;}
		virtual auto& Trace(int i)	->	decltype(us_trace[0])		{return us_trace[i];}
		virtual auto& Peak_x()		->	decltype(us_peak_x)			{return us_peak_x;}
		virtual auto& Peak_y()		->	decltype(us_peak_y)			{return us_peak_y;}
		virtual auto& B_pk_p()		->	decltype(us_b_pk_p			{return us_b_pk_p;}
		virtual auto& B_pk_n()		->	decltype(us_b_pk_n			{return us_b_pk_n;}
		virtual auto& Peak_pos()	->	decltype(us_peak_pos)		{return us_peak_pos;}
		virtual auto& Zero()		->	decltype(d_zero)			{return d_zero;}
		virtual auto& Baseline()	->	decltype(d_baseline)		{return d_baseline;}
		virtual auto& BaseSigma()	->	decltype(d_baseSigma)		{return d_baseSigma;}
		virtual auto& BasePost()	->	decltype(d_basePost)		{return d_basePost;}
		virtual auto& BasePostSigma()->	decltype(d_basePostSigma)	{return d_basePostSigma;}
};

#endif // EVENT_H
