#ifndef METHOD_H // class for processing codes
#define METHOD_H

#ifndef EVENT_H
#include "Event.h"
#endif
#include "TGraph.h"
#include "TF1.h"

using std::unique_ptr;
using std::shared_ptr;
using std::array;
using std::vector;

class Method {
protected:
	int iFailed;
	int iEventlength;
	double dNsPerSample; // scale factors from Digitizer class
	double dVoltsPerBin;
	double dZero;

	enum {n = 0, y = 1, P = 2, ciDFTOrder = 17, ciLAPNpts = 50};
	const double dSlow;
	const double dShigh;

	int iFastTime;
	int iSlowTime; // fast and slow integral window lengths
	int iPGASamples; // PGA method
	int iPGAAverage;

	array<vector<double>,ciDFTOrder> dCos;
	array<vector<double>,ciDFTOrder> dSin;

	array<vector<double>,ciLAPNpts> dExp;
	vector<double> dTrace; // averaged
	array<double,ciLAPNpts> dS;
	array<double,ciLAPNpts> dXform;
	int iLAPAverage;

	array<float,P> fGain;
	int iPeakCut; // 4mV in bins
	int iResolutionScale;
	int iStdLength; // 450ns
	int iStdTrig; // 64ns
	array<vector<double>,P> dStdWave;
	vector<double> dX;
	array<double,P> dStdNorm;
	array<double,P> dStdPeak;
	unique_ptr<TF1> fit_n;
	unique_ptr<TF1> fit_y;
	unique_ptr<TF1> fit_n_f;
	unique_ptr<TF1> fit_y_f;
	unique_ptr<TGraph> graph;
	std::string sConfigDir;

public:
	Method();
	Method(const int length, const int fast, const int slow, const int samples, const array<float,2>& gain, const double scaleT, const double scaleV, shared_ptr<Event> ev);
	~Method();
	void Analyze();
	void SetDCOffset(const short sResolution, const int dc_offset)	{dZero = sResolution*(1.-(double)dc_offset/65535.);} // conversion from wavedump code
	void SetDefaultValues();
	void SetAddresses(const vector<void*>& add);
	void SetConfigDir(const std::string& in) {sConfigDir = in;}
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
	double *dSlowInt;
	double *dFastInt;

	//PGA
	double *dSample;

	//DFT
	array<double,ciDFTOrder> dMagnitude;
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

	double *dXsq_n_f;
	double *dPeakheight_n_f;
	double *dOffset_n_f;
	double *dXsq_y_f;
	double *dPeakheight_y_f;
	double *dOffset_y_f;
	double *dPeak_err_n_f;
	double *dPeak_err_y_f;
	double *dOff_err_n_f;
	double *dOff_err_y_f;

	static float sfVersion;
};

#endif // METHOD_H
