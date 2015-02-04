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
		const int nPar; // 4: peakheight, baseline offset, trigger offset, particle (fixed)
		float gain[P];
		unique_ptr<double[]> pars;
		unique_ptr<int[]> input_wave;
		unique_ptr<int[]> x;
		int std_length; // 450ns
		int std_trig; // 64ns
		int std_peak_x; // 71 ns;
		unique_ptr<double[]> std_wave[P];
		unique_ptr<TF1> fit;
		double std_peak[P];
		double std_base[P];
		double std_norm[P];
		
		static weak_ptr<TTree> tree;
		static int howmany;
		static bool initialized;
		
		static double xsq_n[4]; // only three eljen detectors
		static double peakheight_n[4]; // but may use CH3
		static double baseline_n[4]; // in DT5751DES
		static double offset_n[4];
		static double peak_err_n[4];
		static double base_err_n[4];
		static double offset_err_n[4];
		
		static double xsq_y[4];
		static double peakheight_y[4];
		static double baseline_y[4];
		static double offset_y[4];
		static double peak_err_y[4];
		static double base_err_y[4];
		static double offset_err_y[4];
		
		static int fit_status_n[4];
		static int fit_status_y[4];

	public:
		XSQ(const int ch, const int len, const float gain_in[], const shared_ptr<Digitizer> digitizer);
		virtual ~XSQ();
		virtual void evaluate(const shared_ptr<Event> event);
		static void root_init(shared_ptr<TTree> tree_in);
		static int Std_Wave_init(const shared_ptr<Digitizer> digitizer);
		static int HowMany() {return XSQ::howmany;}
		double fitter(double* x, double* par);
		static float version;
};

#endif // XSQ_H
