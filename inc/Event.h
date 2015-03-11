#ifndef EVENT_H
#define EVENT_H

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
		virtual unsigned short& Trigger()	{return us_trigger;}
		virtual double& Trace(int i)		{return us_trace[i];}
		virtual unsigned short& Peak_x()	{return us_peak_x;}
		virtual double& Peak_y()			{return us_peak_y;}
		virtual double& B_pk_p()			{return us_b_pk_p;}
		virtual double& B_pk_n()			{return us_b_pk_n;}
		virtual double& Peak_pos()			{return us_peak_pos;}
		virtual double& Zero()				{return d_zero;}
		virtual double& Baseline()			{return d_baseline;}
		virtual double& BaseSigma()			{return d_baseSigma;}
		virtual double& BasePost()			{return d_basePost;}
		virtual double& BasePostSigma()		{return d_basePostSigma;}
};

#endif // EVENT_H
