#ifndef XSQ_H // chi-squared test
#define XSQ_H

#ifndef METHOD_H
#include "Method.h"
#endif
#include "TGraph.h"
#include "TF1.h"
#include <vector>

struct fit_results_t { // for later, keeping track of used parameters to eliminate unnecessary evaluations and speed convergence
	double peak;
	double base;
	int offset;
	double xsq;
};

class XSQ : public Method {
	private:
		enum { ciNPar = 3, TF1_pars}; // peakheight, baseline, trigger offset (and particle)
		enum { n = 0, y, P};
		float fGain[P];
		int iPeakCut; // 4mV
		int iResolutionScale;
		int iStdLength; // 450ns
		int iStdTrig; // 64ns
		int iIterations; // 5, usually
		const int ciVariant;
		unique_ptr<double[]> dInputWave; // only used for TF1 variant
		unique_ptr<double[]> dStdWave[P]; // both
		unique_ptr<double[]> dX; // only TF1
		double dConvergence; // Newton's
		double dGradient[ciNPar]; // Newton's
		double dHessianInv[ciNPar][ciNPar]; // Newton's
		double dStdNorm[P];
		double dStdPeak[P];
		double dStep[ciNPar]; // Newton's
		unique_ptr<TF1> fit;
		unique_ptr<TGraph> graph; // only TF1
//		vector<fit_results_t> fit_results_v;
		
		static unique_ptr<TTree> tree;
		static int siHowMany;
		static bool sbInitialized;
		// only three eljen detectors but may use CH3 in DT5751DES
		static double sdXsq_n[4]; // chisquared for neutron
		static double sdPeakheight_n[4]; // peak scale factor
		static double sdBaseline_n[4]; // baseline
		static double sdOffset_n[4]; // trigger shift
		
		static double sdXsq_y[4]; // same, but for gamma
		static double sdPeakheight_y[4];
		static double sdBaseline_y[4];
		static double sdOffset_y[4];
		
		static double sdPeak_err_n[4]; // errors
		static double sdBase_err_n[4];
		static double sdOff_err_n[4];
		static double sdPeak_err_y[4];
		static double sdBase_err_y[4];
		static double sdOff_err_y[4];

	public:
		XSQ();
		XSQ(int ch, int length, shared_ptr<Digitizer> digitizer, int variant);
		virtual ~XSQ();
		virtual void Analyze();
		virtual void SetParameters(void* val, int which, shared_ptr<Digitizer> digitizer);
		static void root_init(TTree* tree_in);
		static void root_fill()		{XSQ::tree->Fill();}
		static void root_write()	{XSQ::tree->Write();}
		static void root_deinit()	{XSQ::tree.reset();} // friending is handled after the fact, writing by the TFile or root_write
		static int HowMany() {return XSQ::siHowMany;}
		double FindChisquare(int p, double dPeak, double dBase, int iOff);
		void SetDefaultParameters();
		double TF1_fit_func(double* x, double* par);
		static float sfVersion;
		enum { VAR_TF1 = 0, VAR_NEW };
		
};

#endif // XSQ_H
