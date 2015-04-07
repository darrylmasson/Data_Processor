#ifndef PROCESSOR_H
#define PROCESSOR_H

#include <fstream>
#include <unistd.h>
#ifndef CCM_H
#include "CCM.h"
#endif
#ifndef DFT_H
#include "DFT.h"
#endif
#ifndef XSQ_H
#include "XSQ.h"
#endif
#ifndef LAP_H
#include "LAP.h"
#endif

class ProcessorException : public exception { // simpler than checking Failed() after every step in setup
	public:
	const char* what() const throw () {
		return "Setup error: ";
	}
};

struct thread_data_t {
	shared_ptr<Event>	event;
	shared_ptr<Method>	methods[NUM_METHODS];
	const bool*			cbpActivated;
};

// Argument is thread_data_t* passed as void*
void* Process(void* arg);

class Processor {
	private:
		shared_ptr<Digitizer> digitizer;
		unique_ptr<TFile> f;
		unique_ptr<TTree> tree;
		unique_ptr<char[]> buffer;
		ifstream fin;
		
		string sConfigFileName;
		string sRawDataFile;
		string sRootFile;
		
		bool			bMethodActive[NUM_METHODS]; // active for this run
		bool			bMethodDone[NUM_METHODS]; // previously processed
		bool			bPositionsSet;
		bool			bRecordTimestamps;
		
		char			cDigName[12];
		char			cMethodNames[NUM_METHODS][12];
		char			cSource[12];
		char			cTreename[NUM_METHODS][4];
		
		unsigned short	usMask; // mask of enabled channels
		
		int				iAverage; // moving average
		int				iChan[MAX_CH]; // only the first iNchan entries used
		int				iEventlength; // samples
		int				iEventsize; // bytes, incl event header
		int				iFailed;
		int				iFastTime[MAX_CH];
		int				iNchan;
		int				iNumEvents;
		int				iPGASamples[MAX_CH];
		int				iSlowTime[MAX_CH];
		int				iSpecial; // special processing options
		int				iTrigPost; // percentage of event after trigger
		int				iXSQ_ndf;
		
		unsigned int	uiDCOffset[MAX_CH]; // DC offset for each channel
		unsigned int	uiThreshold[MAX_CH]; // trigger thresholds
		
		float			fGain[MAX_CH][P]; // for fitter
		float			fDetectorZ[3];
		float			fDetectorR[3];
		
		thread_data_t	td[MAX_CH];
		
	public:
		Processor(int special = -1, int average = 0);
		~Processor(); // handles cleanup
		void BusinessTime(); // it's business, it's business time!
		void ClassAlloc(); // all allocs
		void ConfigTrees(); // includes version checking
		int Failed() {return iFailed;}
		void ParseFileHeader();
		void ParseConfigFile();
		void SetConfigFile(string in) {sConfigFileName = in;}
		void SetDetectorPositions(string in);
		void SetFileSet(string in);
		void SetSource(string in) {strcpy(cSource,in.c_str());}
};

#endif // PROCESSOR_H
