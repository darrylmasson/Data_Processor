#ifndef XSQ_NEW_H // chi-squared test
#define XSQ_NEW_H

#ifndef METHOD_H
#include "Method.h"
#endif

class XSQ_NEW : public Method {
	private:
		enum { ciNPar = 3}; // peakheight, baseline, trigger offset (and particle)
		enum { n = 0, y, P};
		float fGain[P];
		int iPeakCut; // 4mV in bins
		int iResolutionScale;
		int iStdLength; // 450ns
		int iStdTrig; // 64ns
		int iIterations; // 5, usually
		unique_ptr<double[]> dStdWave[P]; // both
		double dConvergence; // Newton's
		double dGradient[ciNPar]; // Newton's
		double dHessianInv[ciNPar][ciNPar]; // Newton's
		double dStdNorm[P];
		double dStdPeak[P];
		double dStep[ciNPar]; // Newton's
		
		static unique_ptr<TTree> tree;
		static int siHowMany;
		static bool sbInitialized;
		// only three eljen detectors but may use CH3 in DT5751DES
		static double sdXsq[2][4]; // chisquared
		static double sdPeakheight[2][4]; // peak scale factor
		static double sdBaseline[2][4]; // baseline
		static double sdOffset[2][4]; // trigger shift
		
		static double sdPeak_err[2][4]; // errors
		static double sdBase_err[2][4];
		static double sdOff_err[2][4];
		
		static short ssIterations[2][4]; // for debugging, mainly
		static double sdConvergence[2][4];

	public:
		XSQ_NEW();
		XSQ_NEW(int ch, int length, shared_ptr<Digitizer> digitizer);
		virtual ~XSQ_NEW();
		virtual void Analyze();
		virtual void SetParameters(void* val, int which, shared_ptr<Digitizer> digitizer);
		static void root_init(TTree* tree_in);
		static void root_fill()		{XSQ_NEW::tree->Fill();}
		static void root_write()	{XSQ_NEW::tree->Write();}
		static void root_deinit()	{XSQ_NEW::tree.reset();} // friending is handled after the fact, writing by the TFile or root_write
		static int HowMany()		{return XSQ_NEW::siHowMany;}
		double FindChisquare(int p, double dPeak, double dBase, int iOff);
		void SetDefaultParameters();
		static float sfVersion;
		
};

#endif // XSQ_NEW_H