#ifndef DFT_H
#define DFT_H

#ifndef METHOD_H
#include "Method.h"
#endif

class DFT : public Method {
	private:
		enum { ciOrder = 17}; // not including zero
		double dScalefactor;
		vector<double> dCos[ciOrder];
		vector<double> dSin[ciOrder];

		static unique_ptr<TTree> tree;
		static int siHowMany;
		static bool sbInitialized;
		
		double dMagnitude[ciOrder];
		double dPhase[ciOrder];
		static double sdOdd[8]; // discrimination parameters
		static double sdEven[8];

		
	public:
		DFT();
		DFT(int ch, int len, shared_ptr<Digitizer> digitizer);
		virtual ~DFT();
		virtual void Analyze();
		virtual void SetParameters(void* val, int which, shared_ptr<Digitizer> digitizer) {}
		static void root_fill()		{DFT::tree->Fill();}
		static void root_write()	{DFT::tree->Write();}
		static void root_init(TTree* tree_in);
		static void root_deinit()	{DFT::tree.reset();} // friending is handled after the fact, writing by the TFile
		static int HowMany() {return DFT::siHowMany;}
		static float sfVersion;
};

#endif // DFT_H
