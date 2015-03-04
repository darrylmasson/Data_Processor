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
		const int ci_nPar; // 4: peakheight, baseline offset, trigger offset, particle (fixed)
		float f_gain[P];
		unique_ptr<double[]> d_pars;
		unique_ptr<int[]> i_input_wave;
		unique_ptr<int[]> i_x;
		int i_std_length; // 450ns
		int i_std_trig; // 64ns
		int i_std_peak_x; // 71 ns;
		unique_ptr<double[]> d_std_wave[P];
		unique_ptr<TF1> fit;
		double d_std_peak[P];
		double d_std_base[P];
		double d_std_norm[P];
		
		static unique_ptr<TTree> tree;
		static int si_howmany;
		static bool sb_initialized;
		
		static double sd_xsq_n[4]; // only three eljen detectors
		static double sd_peakheight_n[4]; // but may use CH3
		static double sd_baseline_n[4]; // in DT5751DES
		static double sd_offset_n[4];
		static double sd_peak_err_n[4];
		static double sd_base_err_n[4];
		static double sd_offset_err_n[4];
		static double sd_prob_n[4];
		
		static double sd_xsq_y[4];
		static double sd_peakheight_y[4];
		static double sd_baseline_y[4];
		static double sd_offset_y[4];
		static double sd_peak_err_y[4];
		static double sd_base_err_y[4];
		static double sd_offset_err_y[4];
		static double sd_prob_y[4];
		
		static int si_fit_status_n[4];
		static int si_fit_status_y[4];

	public:
		XSQ(const int ch, const int len, const float gain_in[], const shared_ptr<Digitizer> digitizer);
		virtual ~XSQ();
		virtual void evaluate(const shared_ptr<Event> event);
		static void root_fill() {XSQ::tree->Fill();}
		static void root_init(TTree* tree_in);
		static TTree* root_deinit() {return XSQ::tree.release();}
		static int Std_Wave_init(const shared_ptr<Digitizer> digitizer);
		static int HowMany() {return XSQ::si_howmany;}
		double fitter(double* x, double* par);
		static float sf_version;
};

#endif // XSQ_H
