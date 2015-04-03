#ifndef XSQ_H // chi-squared test
#define XSQ_H

#ifndef METHOD_H
#include "Method.h"
#endif
#include "TGraph.h"
#include "TF1.h"

class XSQ : public Method {
	private:
		unique_ptr<TGraph> graph;
		const int ciNPar; // 4: peakheight, baseline offset, trigger offset, particle (fixed)
		float fGain[P];
		unique_ptr<double[]> dPars;
		unique_ptr<double[]> dInputWave; // doubles in case we're using averaged events
		unique_ptr<double[]> dX;
		int iStdLength; // 450ns
		int iStdTrig; // 64ns
		int iStdPeakX; // 71 ns;
		unique_ptr<double[]> dStdWave[P];
		unique_ptr<TF1> fit;
		double dStdPeak[P];
		double dStdBase[P];
		double dStdNorm[P];
		
		static unique_ptr<TTree> tree;
		static int siHowMany;
		static bool sbInitialized;
		// only three eljen detectors but may use CH3 in DT5751DES
		static double sdXsq_n[4]; // chisquared for neutron
		static double sdPeakheight_n[4]; // peak scale factor
		static double sdBaseline_n[4]; // baseline shift
		static double sdOffset_n[4]; // trigger shift
		static double sdPeakErr_n[4]; // errors in fit params
		static double sdBaseErr_n[4];
		static double sdOffsetErr_n[4];
		static double sdProb_n[4]; // fit probability
		
		static double sdXsq_y[4]; // same, but for gamma
		static double sdPeakheight_y[4];
		static double sdBaseline_y[4];
		static double sdOffset_y[4];
		static double sdPeakErr_y[4];
		static double sdBaseErr_y[4];
		static double sdOffsetErr_y[4];
		static double sdProb_y[4];
		
		static int siFitStatus_n[4]; // fit status. Usually 4
		static int siFitStatus_y[4];

	public:
		XSQ();
		XSQ(int ch, int length, shared_ptr<Digitizer> digitizer);
		virtual ~XSQ();
		virtual void Analyze();
		virtual void SetParameters(void* val, int which, shared_ptr<Digitizer> digitizer);
		static void root_fill() {XSQ::tree->Fill();}
		static void root_init(TTree* tree_in);
		static TTree* root_deinit() {return XSQ::tree.release();}
		static int Std_Wave_init(const shared_ptr<Digitizer> digitizer);
		static int HowMany() {return XSQ::siHowMany;}
		double fitter(double* x, double* par);
		static float sfVersion;
};

#endif // XSQ_H
