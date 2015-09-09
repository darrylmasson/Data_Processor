#ifndef EVENT_H
#define EVENT_H

#ifndef NGDP_TYPES_H
#include "NGDP_types.h"
#endif

class Event {
	private:
		inline void Average();
		int Peakfinder(bool keep);
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
			Peak_t() : itPeak(nullptr), itStart(nullptr), itEnd(nullptr) {};
			Peak_t(double* ptr) {itPeak = ptr; itStart = ptr; itEnd = ptr;}
			Peak_t(const Peak_t& rhs) : itPeak(rhs.itPeak), itStart(rhs.itStart), itEnd(rhs.itEnd) {};
			~Peak_t() {itPeak = itStart = itEnd = nullptr;}
			bool operator== (const Peak_t& rhs) {return (itPeak == rhs.itPeak);} // x-position
			bool operator!= (const Peak_t& rhs) {return (itPeak != rhs.itPeak);}
			bool operator> (const Peak_t& rhs) {return ((*itStart-*itPeak) > (*rhs.itStart-*rhs.itPeak));} // amplitude
			bool operator< (const Peak_t& rhs) {return ((*itStart-*itPeak) < (*rhs.itStart-*rhs.itPeak));}
			bool operator>= (const Peak_t& rhs) {return (rhs.itPeak-itPeak > 0);} // x-position
			bool operator<= (const Peak_t& rhs) {return (rhs.itPeak-itPeak < 0);}
			Peak_t& operator= (const Peak_t& rhs) {itPeak = rhs.itPeak; itStart = rhs.itStart; itEnd = rhs.itEnd; return *this;}
			Peak_t& operator= (double* ptr) {itPeak = itStart = itEnd = ptr; return *this;}
//			int FWHM();
			double* itPeak;
			double* itStart;
			double* itEnd;
		};

	public:
		Event();
		Event(int eventlength, int baselength, int average, int threshold, int chan, unsigned short* usStart, double* dStart);
		~Event();
		void Analyze();
		int& GetAverage()					{return iAverage;}
		const int& Length()					{return iLength;}
		int Failed()						{return iFailed;}
		void SetAddresses(vector<void*> add);
		void SetScales(double dV, double dT)		{dScaleV = dV; dScaleT = dT;}

		Peak_t Peak; // primary pulse peak
		Peak_t PeakS; // secondary pulse peak for pileup
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
		double* dPeak0s;
		double* dBasePeakP;
		double* dBasePeakN;

		bool* bPileUp;
		bool* bFullWaveform;
		bool* bSaturated;

		short* sDecay;
		short* sRise;
		short* sPeakX;
		short* sPeakXs;
		short* sTrigger;

		const int ciChan;

		static float sfVersion;
};

#endif // EVENT_H
