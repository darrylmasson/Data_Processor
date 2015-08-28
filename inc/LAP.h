#ifndef LAP_H
#define LAP_H

#ifndef METHOD_H
#include "Method.h"
#endif

class LAP : public Method {
	private:
		static unique_ptr<TTree> tree;
		static int siHowMany;
		static bool sbInitialized;
		static double sdLaplaceLow[8];
		static double sdLaplaceHigh[8];
		const static double sdSlow;
		const static double sdShigh;
		const static int ciNpts;
		int iAve;
		double dAveScale;

		vector<vector<double>> dExp;
		vector<double> dTrace;
		vector<double> dS;
		vector<double> dXform;

	public:
		LAP();
		LAP(int ch, int length, shared_ptr<Digitizer> digitizer);
		virtual ~LAP();
		virtual void Analyze();
		virtual void SetParameters(void* val, int which, shared_ptr<Digitizer> digitizer) {}
		virtual void SetEvent(shared_ptr<Event> ev);
		static void root_fill()		{LAP::tree->Fill();}
		static void root_write()	{LAP::tree->Write();}
		static void root_init(TTree* tree_in);
		static void root_deinit()	{LAP::tree.reset();} // friending is handled after the fact, writing by the TFile
		static int HowMany() {return LAP::siHowMany;}
		static float sfVersion;
};

#endif // LAP_H
