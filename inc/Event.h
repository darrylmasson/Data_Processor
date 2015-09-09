#ifndef EVENT_H
#define EVENT_H

#ifndef NGDP_TYPES_H
#include "NGDP_types.h"
#endif

class Event {
	protected:
		inline void Average();
		int Peakfinder();
		int iFailed;
		int iThreshold; // trigger threshold
		int iEventlength;
		int iBaselength;
		int iSpecial;
		int iLength; // number of samples in the waveform
		int iAverage; // 2*x+1 total samples averaged
		unsigned short* const usTrace;
		double dScaleV;
		double dScaleT;
		struct Peak_t {
			double* itPeak;
			double* itStart;
			double* itEnd;
		};

	public:
		Event();
		Event(int eventlength, int baselength, int average, int threshold, int chan, unsigned short* usStart, double* dStart);
		~Event();
		void Analyze();
		inline void PreAnalyze();
		int& GetAverage()					{return iAverage;}
		const int& Length()					{return iLength;}
		int Failed()						{return iFailed;}
		void SetAddresses(vector<void*> add);
		void SetScales(double dV, double dT)		{dScaleV = dV; dScaleT = dT;}

		Peak_t Peak; // primary pulse peak
		double* itBasePkP; // positive peak in baseline samples
		double* itBasePkN; // negative peak
		double* const itBegin; // front of waveform
		double* const itEnd; // end of waveform
		double* itSatEnd; // end of saturation, = dPeakY if not saturated

		double* dBaseline;
		double* dBaseSigma;
		double* dBasePost;
		double* dBasePostSigma;
		double* dIntegral;
		double* dPeak0;
		double* dBasePeakP;
		double* dBasePeakN;

		bool* bPileUp;
		bool* bFullWaveform;
		bool* bSaturated;

		short* sDecay;
		short* sRise;
		short* sPeakX;
		short* sTrigger;

		const int ciChan;

		static float sfVersion;
};

#endif // EVENT_H
