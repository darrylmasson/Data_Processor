#ifndef DISCRIMINATOR_H
#define DISCRIMINATOR_H

#ifndef NGDP_TYPES_H
#include "NGDP_types.h"
#endif

#include "TTree.h"
#include "TFile.h"

class Discriminator {
	private:
		enum { mean = 0, sig1, sig3, NUM_BANDS, NUM_CHANS = 4, iDiscrimBins = 500};

		const double gain[3];

		char			cDiscrimNames[NUM_DISCRIMS][4];
		char			cBandNames[NUM_BANDS][5];

		double 		dDiscrimBand[NUM_DISCRIMS][NUM_CHANS][NUM_BANDS][iDiscrimBins];

		unique_ptr<TTree> T0;
		unique_ptr<TTree> T2;
		unique_ptr<TFile> f;

		bool sbCutPass[NUM_DISCRIMS][NUM_CHANS][NUM_BANDS];
		bool bNewPeak2;

		int			iChan[NUM_CHANS];
		int			iNChan;
		int			iFailed;
		long		lNumEvents;

		double fastint[8];
		double slowint[8];
		double sample[8];
		vector<vector<double>> vPeak2;
		vector<vector<double>>* pPeak2;
		double dPeak2[8]; // for backwards compatibility
		double xsq_n[4];
		double xsq_y[4];
		double peakscale_n[4];
		double peakscale_y[4];
		double baseshift_n[4];
		double baseshift_y[4];
		double dft_even[8];
		double dft_odd[8];
		double lap_high[8];
		double lap_low[8];
		double integral[8];

		double dDiscrim[NUM_DISCRIMS][NUM_CHANS];

public:
		Discriminator();
		~Discriminator();
		void Discriminate();
		int Failed() {return iFailed;}
		void Setup(string filein);

		static float sfVersion;
};

#endif // DISCRIMINATOR_H