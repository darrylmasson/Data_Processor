#ifndef DFT_H
#define DFT_H

#ifndef METHOD_H
#include "Method.h"
#endif

class DFT : public Method {
	private:
		const int ciOrder; // 3 non-zero
		double dScalefactor;
		unique_ptr<double[]> dCos;
		unique_ptr<double[]> dSin;

		static unique_ptr<TTree> tree;
		static int siHowMany;
		static bool sbInitialized;
		
		static double sdMagnitude[8][4]; // order = 3
		static double sdPhase[8][4]; // but we want 0th as well

		
	public:
		DFT(const int ch, const int len, const shared_ptr<Digitizer> digitizer);
		virtual ~DFT();
		virtual void evaluate(const shared_ptr<Event> event);
		static void root_fill() {DFT::tree->Fill();}
		static void root_init(TTree* tree_in);
		static TTree* root_deinit() {return DFT::tree.release();}
		static int HowMany() {return DFT::siHowMany;}
		static float sfVersion;
};

#endif // DFT_H
