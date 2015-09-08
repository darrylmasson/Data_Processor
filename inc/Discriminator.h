#ifndef DISCRIMINATOR_H
#define DISCRIMINATOR_H

#ifndef METHOD_H
#include "Method.h"
#endif

#include "TTree.h"

class Discriminator {
	private:
		enum { mean = 0, sig1, sig3, NUM_BANDS};

		const double gain[3];

		char			cDiscrimNames[NUM_DISCRIMS][4];
		char			cBandNames[NUM_BANDS][5];

		vector<double> 		vDiscrimBand[NUM_DISCRIMS][4][NUM_BANDS];

		static unique_ptr<TTree> tree;

		static bool sbCCMCutPass[4][NUM_BANDS];
		static bool sbPGACutPass[4][NUM_BANDS];
		static bool sbNGMCutPass[4][NUM_BANDS];
		static bool sbWBSCutPass[4][NUM_BANDS];
		static bool sbLAPCutPass[4][NUM_BANDS];
		static bool sbDFTCutPass[4][NUM_BANDS];

		static bool sbCutPass[NUM_DISCRIMS][4][NUM_BANDS];

		int			iChan[4];
		int			iDiscrimBins;
		int			iFailed;
		const int	ch;
		const int	ciChan;

		const double* fastint;
		const double* slowint;
		const double* sample;
		const double* peakheight2;
		const double* xsq_n;
		const double* xsq_y;
		const double* peakscale_n;
		const double* peakscale_y;
		const double* baseshift_n;
		const double* baseshift_y;
		const double* dft_even;
		const double* dft_odd;
		const double* lap_high;
		const double* lap_low;
		const double* integral;

		static double dDiscrim[NUM_DISCRIMS][4];

public:
		Discriminator(int channel);
		~Discriminator();
		void SetDiscriminationValue();
		void Discriminate();
		int Failed() {return iFailed;}
		void SetAddresses(vector<void*> add);
		static void FriendshipIsMagic();
		static void CutsTree_init(TTree* tree_cuts);
		static void Cuts_fill();
		static void Cuts_write()		{Discriminator::tree->Write();}
		static void Cuts_deinit()		{Discriminator::tree.reset();}

		static float sfVersion;
};

#endif // DISCRIMINATOR_H