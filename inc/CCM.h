#ifndef CCM_H // charge comparison method
#define CCH_H

#ifndef METHOD_H
#include "Method.h"
#endif

class CCM : public Method { // also includes PGA
	private:
		int iFastTime;
		int iSlowTime; // fast and slow integral window lengths
		int iPGASamples; // PGA method
		
		static unique_ptr<TTree> tree;
		static int siHowMany;
		static bool sbInitialized;
		
		static bool sbFullWave[8]; // waveform decays before end of event
		static bool sbSaturated[8]; // voltage saturation on digitizer
		static bool sbTruncated[8]; // slow end > eventlength
		
		static short ssDecay[8]; // decay time
		static short ssFastStop[8]; // stop point for fast integral
		static short ssRise[8]; // risetime
		static short ssPeakX[8]; // x location of peak
		static short ssSlowStop[8]; // stop point for slow integral
		
		static double sdBaseline[8];
		static double sdBaseSigma[8];
		static double sdBasePost[8];
		static double sdBasePostSigma[8]; // baseline stuff from Event
		static double sdBasePeakP[8];
		static double sdBasePeakN[8];
		static double sdFullInt[8];
		static double sdFastInt[8]; // integral values
		static double sdSlowInt[8];
		static double sdPeak0[8];
		static double sdPeak1[8]; // peakheights
		static double sdPeak2[8];
		static double sdPeakP[8];
		
		static double sdGradient[8]; // PGA value
	
	public:
		CCM();
		CCM(int ch, int length, shared_ptr<Digitizer> digitizer);
		virtual ~CCM();
		virtual void Analyze();
		virtual void SetParameters(void* val, int which, shared_ptr<Digitizer> digitizer);
		static void root_fill()		{CCM::tree->Fill();}
		static void root_init(TTree* tree_in);
		static void root_write()	{CCM::tree->Write();}
		static void root_deinit()	{CCM::tree.reset();} // friending is handled after the fact, writing by the TFile
		static int HowMany() {return CCM::siHowMany;}
		static float sfVersion;
};

#endif // CCM_H
