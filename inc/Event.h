#ifndef EVENT_H
#define EVENT_H

#ifndef DIGITIZER_H
#include "Digitizer.h"
#endif

class Event {
	protected: // these are inheritable
		int iFailed;
		int iThreshold; // trigger threshold
		int iEventlength;
		int iBaselength;
		int iSpecial;
		const int ciSamples; // number of samples of pulse in memory (useful for Special or Average instances)
		static int siLength; // number of samples in the waveform (generally not strictly ciSamples)
		static int siHowMany;
		
		unsigned short usPeakY; // primary pulse peak, y coord.
		unsigned short usBasePkP; // positive peak in baseline samples
		unsigned short usBasePkN; // negative peak
		unsigned short usPeakPos; // positive peak
		unsigned short* uspTrace;
		unsigned short usTrigger; // triggering sample
		unsigned short usPeakX; // primary pulse peak, x coord
		double dZero; // ADC bin for ground
		double dBaseline;
		double dBaseSigma;
		double dBasePost;
		double dBasePostSigma;
		double dBaseScale; // 1/baselength

	public:
		Event();
		Event(int len, shared_ptr<Digitizer> digitizer);
		~Event();
		virtual void SetAverage(int average)										{}
		virtual void SetDCOffset(shared_ptr<Digitizer> digitizer, int dc_offset)	{dZero = digitizer->Resolution()*(1.-(double)dc_offset/65535.);}
		virtual void SetThreshold(int threshold)									{iThreshold = threshold;}
		virtual void SetTrace(unsigned short* trace)								{uspTrace = trace;}
		virtual void Analyze();
		int Failed() {return iFailed;}
		static const int& Length()			{return Event::siLength;}
		virtual unsigned short& Trigger()	{return usTrigger;}
		virtual double Trace(int i)			{return ((i < 0) ? dBaseline : (i < ciSamples ? uspTrace[i] : dBasePost));} // hopefully the most efficient way to do the checks
		virtual unsigned short& Peak_x()	{return usPeakX;}
		virtual double Peak_y()				{return usPeakY;}
		virtual double B_pk_p()				{return usBasePkP;}
		virtual double B_pk_n()				{return usBasePkN;}
		virtual double Peak_pos()			{return usPeakPos;}
		virtual double& Zero()				{return dZero;}
		virtual double& Baseline()			{return dBaseline;}
		virtual double& BaseSigma()			{return dBaseSigma;}
		virtual double& BasePost()			{return dBasePost;}
		virtual double& BasePostSigma()		{return dBasePostSigma;}
};

#endif // EVENT_H
