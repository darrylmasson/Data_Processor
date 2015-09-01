#ifndef EVENT_H
#define EVENT_H

#ifndef NGDP_TYPES_H
#include "NGDP_types.h"
#endif

class Event {
	protected: // these are inheritable
		int iFailed;
		int iThreshold; // trigger threshold
		int iEventlength;
		int iBaselength;
		int iSpecial;
		int iLength; // number of samples in the waveform
		int iAverage; // 2*x+1 total samples averaged
		vector<unsigned short> uspTrace;

	public:
		Event();
		Event(int eventlength, int baselength, int average, unsigned short* start, unsigned short* end, int threshold, int chan);
		~Event();
		virtual void SetDCOffset(short sResolution, int dc_offset)	{dZero = sResolution*(1.-(double)dc_offset/65535.);}
		virtual void Analyze();
		virtual inline void PreAnalyze();
		virtual int& GetAverage()			{return iAverage;}
		const int& Length()					{return iLength;}
		int Failed()						{return iFailed;}

		vector<double> vTrace;
		vector<double>::iterator itPeakY; // primary pulse peak
		vector<double>::iterator itBasePkP; // positive peak in baseline samples
		vector<double>::iterator itBasePkN; // negative peak
		vector<double>::iterator itTrigger; // triggering sample
		vector<double>::iterator itBegin; // front of waveform
		vector<double>::iterator itEnd; // end of waveform
		vector<double>::iterator itSatEnd; // end of saturation, = dPeakY if not saturated
		vector<double>::iterator itPulseStart;
		vector<double>::iterator itPulseEnd;
		double dZero; // ADC bin for ground
		double dBaseline;
		double dBaseSigma;
		double dBasePost;
		double dBasePostSigma;
		double dIntegral;

		bool bPileUp;
		bool bFullWaveform;
		bool bSaturated;

		const int ciChan;
};

#endif // EVENT_H
