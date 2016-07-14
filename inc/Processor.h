#ifndef PROCESSOR_H
#define PROCESSOR_H

#include <fstream>

#ifndef METHOD_H
#include "Method.h"
#endif

#include "TFile.h"
#include "TTree.h"

class ProcessorException : public exception { // simpler than checking Failed() after every step in setup
	public:
	const char* what() const throw () {
		return "Setup error: ";
	}
};

class Processor {
	private:
		struct {
			char			cName[12];
			double			dSamplerate;
			short			sResolution;
			double			dVpp;
			double			dVoltsPerBin;
			double			dNsPerSample;
			int				iBaselength;
			int				iSpecial;
			dig_id_t		id;
		} digitizer;
		unique_ptr<TFile> f;
		unique_ptr<TTree> T0;
		unique_ptr<TTree> T1;
		unique_ptr<TTree> TS;
		vector<char> buffer;
		ifstream fin;

		string sConfigFileName;
		string sRawDataFile;
		string sRootFile;
		string sName;
		string sPositions;
		string sIODir;
		string sConfigDir;

		vector<vector<double> > dTrace;
		vector<shared_ptr<Event> > event;
		vector<shared_ptr<Method> > method;

		bool			bPositionsSet;
		bool			bForceOldFormat;

		char			cSource[12];

		unsigned short	usMask; // mask of enabled channels

		int					iAverage; // moving average
		vector<int>			iChan; // size == iNchan
		int					iEventlength; // samples
		int					iEventsize; // bytes, incl event header
		int					iFailed;
		int					iFastTime[MAX_CH];
		int					iLevel; // level of processing to be done
		int					iNchan; // number of enabled channels
		int					iPGASamples[MAX_CH];
		int					iSlowTime[MAX_CH];
		int					iTrigPost; // percentage of event after trigger
		int					iXSQ_ndf;

		unsigned int		uiDCOffset[MAX_CH]; // DC offset for each channel
		unsigned int		uiThreshold[MAX_CH]; // trigger thresholds

		array<array<float,2>,MAX_CH>	fGain; // for fitter
		float							fDetectorZ[3];
		float							fDetectorR[3];
		float							fHV;
		float							fCurrent;

		// These in T0
		bool bFullWave[8]; // waveform decays before end of event
		bool bSaturated[8]; // voltage saturation on digitizer

		short sDecay[8]; // decay time
		vector<vector<double>> sRise; // risetime
		vector<vector<double>>* pRise;
		vector<vector<double>> sPeakX; // x location of peak
		vector<vector<double>>* pPeakX;
		short sTrigger[8];
		vector<vector<double>> sHWHM;
		vector<vector<double>>* pHWHM;
		short sSaturation[8]; // length of saturation in ns
		short sNumPeaks[8];

		double dBase[8]; // these values in samples
		double dSigma[8];
		double dBaseP[8];
		double dBasePS[8];
		vector<vector<double>> dPeak0;
		vector<vector<double>>* pPeak0;
		vector<vector<double>> dPeak2;
		vector<vector<double>>* pPeak2;
		double dFullInt[8]; // except this
		double dBasePeakP[8];
		double dBasePeakN[8];

		// These in T1
		bool bTruncated[8]; // slow end >= eventlength

		double dBaseline[8];
		double dBaseSigma[8];
		double dBasePost[8];
		double dBasePostSigma[8]; // baseline stuff in volts
		double dFastInt[8]; // integral values
		double dSlowInt[8];

		double dSample[8]; // sample used for PGA

		double dOdd[8]; // discrimination parameters
		double dEven[8];

		double dLaplaceLow[8];
		double dLaplaceHigh[8];

		double dXsq_n[4]; // chisquared
		double dXsq_y[4];
		double dXsq_f_n[4]; // with fixed baseline
		double dXsq_f_y[4];
		double dPeak_scale_n[4]; // peak scale factor
		double dPeak_scale_y[4];
		double dPeak_scale_f_n[4];
		double dPeak_scale_f_y[4];
		double dBase_shift_n[4]; // baseline
		double dBase_shift_y[4];
		double dOffset_n[4]; // trigger shift
		double dOffset_y[4];
		double dOffset_f_n[4];
		double dOffset_f_y[4];

		double dPeak_err_n[4]; // errors
		double dPeak_err_y[4];
		double dPeak_err_f_n[4];
		double dPeak_err_f_y[4];
		double dBase_err_n[4];
		double dBase_err_y[4];
		double dOff_err_n[4];
		double dOff_err_y[4];
		double dOff_err_f_n[4];
		double dOff_err_f_y[4];


	public:
		Processor();
		~Processor();
		void BusinessTime(); // it's business, it's business time!
		int Failed()								{return iFailed;}
		void ForceOld() {bForceOldFormat = true;}
		vector<void*> GetAddresses(const int ch, const int level);
		void SetConfigFile(const string& in)				{sConfigFileName = in;}
		void SetDetectorPositions(const string& in);
		void SetNGSetpoint(const float HV, const float Current) {fHV = HV; fCurrent = Current;}
		void SetSource(const string& in)					{strcpy(cSource,in.c_str());}
		void SetParams(const int average, const int level)	{iAverage = average; iLevel = level;}
		void SetIODir(const string& in)						{sIODir = in;}
		void SetConfigDir(const string& in)					{sConfigDir = in;}
		void Setup(const string& in); // opens files, parses header and config file, does trees and alloc'ing

		int iNumEvents; // public for average rate counter
};

#endif // PROCESSOR_H
