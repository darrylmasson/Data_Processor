#ifndef METHOD_H // class for processing codes
#define METHOD_H

#ifndef EVENT_H
#include "Event.h"
#endif
#include "TGraph.h"
#include "TF1.h"

class Method {
protected:
	int iFailed;
	int iEventlength;
	double dScaleT; // scale factors from Digitizer class
	double dScaleV;
	double dZero;

	enum {n = 0, y = 1, P = 2, ciDFTOrder = 17, ciLAPNpts = 50};
	const double dSlow;
	const double dShigh;

	int iFastTime;
	int iSlowTime; // fast and slow integral window lengths
	int iPGASamples; // PGA method
	int iPGAAverage;

	vector<double> dCos[ciDFTOrder];
	vector<double> dSin[ciDFTOrder];

	vector<double> dExp[ciLAPNpts];
	vector<double> dTrace; // averaged
	vector<double> dS;
	vector<double> dXform;
	int iLAPAverage;

	float fGain[P];
	int iPeakCut; // 4mV in bins
	int iResolutionScale;
	int iStdLength; // 450ns
	int iStdTrig; // 64ns
	vector<double> dStdWave[P];
	vector<double> dX;
	double dStdNorm[P];
	double dStdPeak[P];
	unique_ptr<TF1> fit;
	unique_ptr<TGraph> graph;

public:
	Method();
	Method(int length, int fast, int slow, int samples, float gain[2], double scaleT, double scaleV, shared_ptr<Event> ev);
	~Method();
	void Analyze();
	void SetDCOffset(short sResolution, int dc_offset)	{dZero = sResolution*(1.-(double)dc_offset/65535.);}
	void SetDefaultValues();
	void SetAddresses(vector<void*> add);
	int Failed() {return iFailed;}
	double TF1_fit_func(double* x, double* par);
	shared_ptr<Event> event;

	//CCM
	bool *bTruncated;
	double *dBaseline;
	double *dBaseSigma;
	double *dBasePost;
	double *dBasePostSigma;
	double *dBasePeakN;
	double *dBasePeakP;
	double *dPeak1;
	double *dPeak2;
	double *dSlowInt;
	double *dFastInt;

	//PGA
	double *dSample;

	//DFT
	double dMagnitude[ciDFTOrder];
	double *dOdd;
	double *dEven;

	//LAP
	double *dLaplaceHigh;
	double *dLaplaceLow;

	//NGM
	double *dXsq_n;
	double *dPeakheight_n;
	double *dBaseline_n;
	double *dOffset_n;
	double *dXsq_y;
	double *dPeakheight_y;
	double *dBaseline_y;
	double *dOffset_y;
	double *dPeak_err_n;
	double *dPeak_err_y;
	double *dBase_err_n;
	double *dBase_err_y;
	double *dOff_err_n;
	double *dOff_err_y;

	static float sfVersion;
};

#endif // METHOD_H
