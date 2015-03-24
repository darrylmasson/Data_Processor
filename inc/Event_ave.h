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
		Event_ave(int len, std::shared_ptr<Digitizer> dig, int dc_offset, int threshold_in, int average_in);
		~Event_ave();
		virtual void Set(unsigned short* in);
	//	virtual unsigned short& Trigger()	{return usTrigger;}
		virtual double Trace(int i)			{return dTrace[i];}
	//	virtual unsigned short& Peak_x()	{return usPeakX;}
		virtual double Peak_y()				{return dPeakY;}
		virtual double B_pk_p()				{return dBasePkP;}
		virtual double B_pk_n()				{return dBasePkN;}
		virtual double Peak_pos()			{return dPeakPos;}
	//	virtual auto Zero()					{return dZero;}
	//	virtual auto Baseline()				{return dBaseline;}
	//	virtual auto BaseSigma()			{return dBaseSigma;}
	//	virtual auto BasePost()				{return dBasePost;}
	//	virtual auto BasePostSigma()		{return dBasePostSigma;}
};

#endif // EVENT_AVE_H
