#ifndef XSQ_TF1_H // chi-squared test
#define XSQ_TF1_H

#ifndef METHOD_H
#include "Method.h"
#endif
#include "TGraph.h"
#include "TF1.h"

class XSQ_TF1 : public Method {
	private:
		enum { ciNPar = 4}; // peakheight, baseline, trigger offset, and particle
		enum { n = 0, y, P};
		float fGain[P];
		int iPeakCut; // 4mV in bins
		int iResolutionScale;
		int iStdLength; // 450ns
		int iStdTrig; // 64ns
		int iIterations; // 5, usually
		unique_ptr<double[]> dInputWave; // only used for TF1 variant
		unique_ptr<double[]> dStdWave[P]; // both
		unique_ptr<double[]> dX; // only TF1
		double dStdNorm[P];
		double dStdPeak[P];
		unique_ptr<TF1> fit;
		unique_ptr<TGraph> graph; // only TF1
		
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
		XSQ_TF1();
		XSQ_TF1(int ch, int length, shared_ptr<Digitizer> digitizer);
		virtual ~XSQ_TF1();
		virtual void Analyze();
		virtual void SetParameters(void* val, int which, shared_ptr<Digitizer> digitizer);
		static void root_init(TTree* tree_in);
		static void root_fill()		{XSQ_TF1::tree->Fill();}
		static void root_write()	{XSQ_TF1::tree->Write();}
		static void root_deinit()	{XSQ_TF1::tree.reset();} // friending is handled after the fact, writing by the TFile or root_write
		static int HowMany()		{return XSQ_TF1::siHowMany;}
		void SetDefaultParameters();
		double TF1_fit_func(double* x, double* par);
		static float sfVersion;
};

#endif // XSQ_TF1_H
