#ifndef DFT_H
#define DFT_H

#ifndef METHOD_H
#include "Method.h"
#endif

class DFT : public Method {
	private:
		const int ci_order; // 3 non-zero
		double d_scalefactor;
		unique_ptr<double[]> d_COS;
		unique_ptr<double[]> d_SIN;

		static unique_ptr<TTree> tree;
		static int si_howmany;
		static bool sb_initialized;
		
		static double sd_magnitude[8][4]; // order = 3
		static double sd_phase[8][4]; // but we want 0th as well

		
	public:
		DFT(const int ch, const int len, const shared_ptr<Digitizer> digitizer);
		virtual ~DFT();
		virtual void evaluate(const shared_ptr<Event> event);
		static void root_fill() {DFT::tree->Fill();}
		static void root_init(TTree* tree_in);
		static TTree* root_deinit() {return DFT::tree.release();}
		static int HowMany() {return DFT::si_howmany;}
		static float sf_version;
};

#endif // DFT_H
