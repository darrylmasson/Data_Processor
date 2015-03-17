#ifndef LAP_H
#define LAP_H

#ifndef METHOD_H
#include "Method.h"
#endif

class LAP : public Method {
	private:
		static unique_ptr<TTree> tree;
		static int si_howmany;
		static bool sb_initialized;
		static double sd_laplace[8];
		static double sd_longint[8]; // integral over the full event window
		
		unique_ptr<double[]> d_EXP;
		
	public:
		LAP(int ch, int len, shared_ptr<Digitizer> digitizer);
		virtual ~LAP();
		virtual void evaluate(const shared_ptr<Event> event);
		static void root_fill() {LAP::tree->Fill();}
		static void root_init(TTree* tree_in);
		static TTree* root_deinit() {return LAP::tree.release();}
		static int HowMany() {return LAP::si_howmany;}
		static float sf_version;
		static float sf_s;
};

#endif // LAP_H
