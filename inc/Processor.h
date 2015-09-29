#ifndef PROCESSOR_H
#define PROCESSOR_H

#include <fstream>
#include <unistd.h>
#ifndef METHOD_H
#include "Method.h"
#endif
#ifndef DISCRIMINATOR_H
#include "Discriminator.h"
#endif
#include "TFile.h"

class ProcessorException : public exception { // simpler than checking Failed() after every step in setup
	public:
	const char* what() const throw () {
		return "Setup error: ";
	}
};

class Processor {
	private:
		struct {
			char		cName[12];
			double		dSamplerate;
			short		sResolution;
			double		dVpp;
			double		dScaleT;
			double		dScaleV;
			int			iBaselength;
			int			iSpecial;
			dig_id_t	id;
		} digitizer;
		unique_ptr<TFile> f;
		unique_ptr<TTree> T0;
		unique_ptr<TTree> T1;
		unique_ptr<TTree> TS; // discriminator handles its own tree
		unique_ptr<char[]> buffer;
		ifstream fin;

		string sConfigFileName;
		string sRawDataFile;
		string sRootFile;

		unique_ptr<double[]> dTrace[MAX_CH];
		shared_ptr<Event> event[MAX_CH];
		shared_ptr<Method> method[MAX_CH];
		shared_ptr<Discriminator> discriminator[MAX_CH];

		bool			bPositionsSet;

		char			cBuildID[21];
		char			cSource[12];

		unsigned short	usMask; // mask of enabled channels

		int				iAverage; // moving average
		int				iChan[MAX_CH]; // only the first iNchan entries used
		int				iEventlength; // samples
		int				iEventsize; // bytes, incl event header
		int				iFailed;
		int				iFastTime[MAX_CH];
		int				iLevel; // level of processing to be done
		int				iNchan; // number of enabled channels
		int				iNumEvents;
		int				iPGASamples[MAX_CH];
		int				iSlowTime[MAX_CH];
		int				iSpecial; // special processing options
		int				iTrigPost; // percentage of event after trigger
		int				iXSQ_ndf;

		unsigned int	uiDCOffset[MAX_CH]; // DC offset for each channel
		unsigned int	uiThreshold[MAX_CH]; // trigger thresholds

		float			fGain[MAX_CH][2]; // for fitter
		float			fDetectorZ[3];
		float			fDetectorR[3];

		// These in T0
		bool bFullWave[8]; // waveform decays before end of event
		bool bSaturated[8]; // voltage saturation on digitizer

		short sDecay[8]; // decay time
		short sRise[8]; // risetime
		short sPeakX[8]; // x location of peak
		short sPeakXs[8]; // pile-up peak
		short sTrigger[8];
		short sPeaks[8]; // number of peaks found

		double dBase[8]; // these values in samples
		double dSigma[8];
		double dBaseP[8];
		double dBasePS[8];
		double dPeak0[8];
		double dPeak0s[8];
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

		double dPeak1[8]; // peakheights
		double dPeak2[8];

		double dSample[8]; // sample used for PGA

		double dOdd[8]; // discrimination parameters
		double dEven[8];

		double dLaplaceLow[8];
		double dLaplaceHigh[8];

		double dXsq[2][4]; // chisquared
		double dPeak_scale[2][4]; // peak scale factor
		double dBase_shift[2][4]; // baseline
		double dOffset[2][4]; // trigger shift

		double dPeak_err[2][4]; // errors
		double dBase_err[2][4];
		double dOff_err[2][4];


	public:
		Processor();
		~Processor();
		void BusinessTime(); // it's business, it's business time!
		int Failed()								{return iFailed;}
		vector<void*> SetAddresses(int ch, int level);
		void SetConfigFile(string in)				{sConfigFileName = in;}
		void SetDetectorPositions(string in);
		void SetSource(string in)					{strcpy(cSource,in.c_str());}
		void SetSpecials(int special, int average)	{iSpecial = special; iAverage = average;}
		void Setup(string in); // opens files, parses header and config file, does trees and alloc'ing
};

#endif // PROCESSOR_H
