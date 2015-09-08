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
		double dScaleInt; // V*ns/(bin*clock cycle)
		double dScaleV;
		double dScaleT;
		struct Peak_t {
			vector<double>::iterator itPeak;
			vector<double>::iterator itStart;
			vector<double>::iterator itEnd;
		};

	public:
		Event();
		Event(int eventlength, int baselength, int average, unsigned short* start, int threshold, int chan);
		~Event();
		void Analyze();
		inline void PreAnalyze();
		int& GetAverage()					{return iAverage;}
		const int& Length()					{return iLength;}
		int Failed()						{return iFailed;}
		void SetAddresses(vector<void*> add);
		void SetScales(double dInt, double dV, double dT)		{dScaleInt = dInt; dScaleV = dV; dScaleT = dT;}

		vector<double> vTrace;
		Peak_t Peak; // primary pulse peak
		vector<double>::iterator itBasePkP; // positive peak in baseline samples
		vector<double>::iterator itBasePkN; // negative peak
		vector<double>::iterator itBegin; // front of waveform
		vector<double>::iterator itEnd; // end of waveform
		vector<double>::iterator itSatEnd; // end of saturation, = dPeakY if not saturated

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
