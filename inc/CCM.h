#ifndef CCM_H // charge comparison method
#define CCH_H

#ifndef METHOD_H
#include "Method.h"
#endif

class CCM : public Method { // also includes PGA
	private:
		int fastTime;
		int slowTime;
		int gradSamples;
		double scaleT;
		double scaleV;
		
		static weak_ptr<TTree> tree;
		static int howmany;
		static bool initialized;
		
		static bool fullwave[8];
		static bool saturated[8]; // up to 8 channels
		static bool truncated[8];
		
		static short decay[8];
		static short fstop[8];
		static short rise[8];
		static short peakx[8];
		static short sstop[8];
		
		static double baseline[8];
		static double baseSigma[8];
		static double basePost[8];
		static double basePostSigma[8];
		static double basePeakp[8];
		static double basePeakn[8];
		static double fullint[8];
		static double fastint[8];
		static double slowint[8];
		static double peak0[8];
		static double peak1[8];
		static double peak2[8];
		static double peakp[8];
		
		static double gradient[8];
	
	public:
		CCM(const int ch, const int fast, const int slow, const int samples, const shared_ptr<Digitizer> digitizer);
		virtual ~CCM();
		virtual void evaluate(const shared_ptr<Event> event);
		static void root_init(shared_ptr<TTree> tree_in);
		static int HowMany() {return CCM::howmany;}
		static float version;
};

#endif // CCM_H
