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
		static int siLength; // number of samples in the waveform
		static int siHowMany;
		
		unsigned short usPeakY; // everything measured in ADC counts
		unsigned short usBasePkP; // peak in baseline samples
		unsigned short usBasePkN;
		unsigned short usPeakPos; // positive peak
		unsigned short* usTrace;
		unsigned short usTrigger; // triggering sample
		unsigned short usPeakX;
		double dZero; // ADC bin for ground
		double dBaseline;
		double dBaseSigma;
		double dBasePost;
		double dBasePostSigma;

	public:
		Event() {++Event::siHowMany;} // for Event_ave
		Event(int len, std::shared_ptr<Digitizer> dig, int dc_offset, int threshold_in);
		~Event();
		virtual void Set(unsigned short* in);
		int Failed() {return iFailed;}
		static const int& Length() {return Event::siLength;}
		std::shared_ptr<Digitizer> digitizer;
		virtual unsigned short& Trigger()	{return usTrigger;}
		virtual double Trace(int i)			{return usTrace[i];}
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
