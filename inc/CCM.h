#ifndef CCM_H // charge comparison method
#define CCH_H

#ifndef METHOD_H
#include "Method.h"
#endif

class CCM : public Method { // also includes PGA
	private:
		int i_fastTime;
		int i_slowTime; // fast and slow integral window lengths
		int i_gradSamples; // PGA method
		
		static unique_ptr<TTree> tree;
		static int si_howmany;
		static bool sb_initialized;
		
		static bool sb_fullwave[8]; // waveform decays before end of event
		static bool sb_saturated[8]; // voltage saturation on digitizer
		static bool sb_truncated[8]; // slow end > eventlength
		
		static short ss_decay[8]; // decay time
		static short ss_fstop[8]; // stop point for fast integral
		static short ss_rise[8]; // risetime
		static short ss_peakx[8]; // x location of peak
		static short ss_sstop[8]; // stop point for slow integral
		
		static double sd_baseline[8];
		static double sd_baseSigma[8];
		static double sd_basePost[8];
		static double sd_basePostSigma[8]; // baseline stuff from Event
		static double sd_basePeakp[8];
		static double sd_basePeakn[8];
		static double sd_fullint[8];
		static double sd_fastint[8]; // integral values
		static double sd_slowint[8];
		static double sd_peak0[8];
		static double sd_peak1[8]; // peakheights
		static double sd_peak2[8];
		static double sd_peakp[8];
		
		static double sd_gradient[8]; // PGA value
	
	public:
		CCM(const int ch, const int fast, const int slow, const int samples, const shared_ptr<Digitizer> digitizer);
		virtual ~CCM();
		virtual void evaluate(const shared_ptr<Event> event);
		static void root_fill() {CCM::tree->Fill();}
		static void root_init(TTree* tree_in);
		static TTree* root_deinit() {return CCM::tree.release();} // returns the TTree pointer to Processor() for finishing and stuff
		static int HowMany() {return CCM::si_howmany;}
		static float sf_version;
};

#endif // CCM_H
