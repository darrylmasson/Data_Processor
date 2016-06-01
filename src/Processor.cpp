#include "Processor.h"
#include <ctime>
#include <pthread.h>
#include <ratio>
#include <chrono>
#include <csignal>

using namespace std::chrono;
using std::cout;

static int s_interrupted = 0;
static void s_signal_handler(int signal_value) {
	s_interrupted = 1;
}

static void s_catch_signals() {
	struct sigaction action;
	action.sa_handler = s_signal_handler;
	action.sa_flags = 0;
	sigemptyset(&action.sa_mask);
	sigaction(SIGINT, &action, nullptr);
	sigaction(SIGTERM, &action, nullptr);
}

long UnixConverter(const string& in) {
	//Read in file time
	tm timeinfo;
	try {
		timeinfo.tm_year = stoi(in.substr(0,2)) + 100; // counts from 1900
		timeinfo.tm_mon = stoi(in.substr(2,2)) - 1; // jan == 0
		timeinfo.tm_mday = stoi(in.substr(4,2));
		timeinfo.tm_hour = stoi(in.substr(7,2));
		timeinfo.tm_min = stoi(in.substr(9.2));
		timeinfo.tm_sec = 0;
	} catch (std::exception& e) {
		cout << e.what() << '\n';
		return -1;
	}

	return mktime(&timeinfo);

}

Processor::Processor() {

	iAverage = 0;
	iLevel = 0;

	sConfigFileName = "\0";
	sRawDataFile = "\0";
	sRootFile = "\0";
	sName = "\0";
	cSource[0] = '\0';

	bForceOldFormat = false;
	s_catch_signals();
	memset(&digitizer, 0, sizeof(digitizer));

	usMask = 0;

	iEventlength = 0;
	iEventsize = 0;
	iFailed = 0;
	memset(iFastTime, 0, sizeof(iFastTime));
	iLevel = 0;
	iNchan = 0;
	iNumEvents = 0;
	memset(iPGASamples, 0, sizeof(iPGASamples));
	memset(iSlowTime, 0, sizeof(iSlowTime));
	iTrigPost = 0;
	iXSQ_ndf = 0;

	fHV = 0;
	fCurrent = 0;

	memset(uiDCOffset, 0, sizeof(uiDCOffset));
	memset(uiThreshold, 0, sizeof(uiThreshold));

	for (auto& a : fGain) a.fill(0);
	memset(fDetectorZ, 0, sizeof(fDetectorZ));
	memset(fDetectorR, 0, sizeof(fDetectorR));

	// arrays for trees
	memset(bFullWave,		0,	sizeof(bFullWave));
	memset(bSaturated,		0,	sizeof(bSaturated));
	memset(sDecay,			0,	sizeof(sDecay));
	memset(sTrigger,		0,	sizeof(sTrigger));
	memset(sSaturation,		0,	sizeof(sSaturation));
	memset(sNumPeaks,		0,	sizeof(sNumPeaks));
	memset(dBase,			0,	sizeof(dBase));
	memset(dSigma,			0,	sizeof(dSigma));
	memset(dBaseP,			0,	sizeof(dBaseP));
	memset(dBasePS,			0,	sizeof(dBasePS));
	memset(dFullInt,		0,	sizeof(dFullInt));
	memset(dBasePeakP,		0,	sizeof(dBasePeakP));
	memset(dBasePeakN,		0,	sizeof(dBasePeakN));

	memset(bTruncated,		0,	sizeof(bTruncated));
	memset(dBaseline,		0,	sizeof(dBaseline));
	memset(dBaseSigma,		0,	sizeof(dBaseSigma));
	memset(dBasePost,		0,	sizeof(dBasePost));
	memset(dBasePostSigma,	0,	sizeof(dBasePostSigma));
	memset(dFastInt,		0,	sizeof(dFastInt));
	memset(dSlowInt,		0,	sizeof(dSlowInt));
	memset(dSample,			0,	sizeof(dSample));
	memset(dOdd,			0,	sizeof(dOdd));
	memset(dEven,			0,	sizeof(dEven));
	memset(dLaplaceLow,		0,	sizeof(dLaplaceLow));
	memset(dLaplaceHigh,	0,	sizeof(dLaplaceHigh));
	memset(dXsq_n,			0,	sizeof(dXsq_n));
	memset(dXsq_y,			0,	sizeof(dXsq_y));
	memset(dXsq_f_n,		0,	sizeof(dXsq_f_n));
	memset(dXsq_f_y,		0,	sizeof(dXsq_f_y));
	memset(dPeak_scale_n,	0,	sizeof(dPeak_scale_n));
	memset(dPeak_scale_y,	0,	sizeof(dPeak_scale_y));
	memset(dPeak_scale_f_n,	0,	sizeof(dPeak_scale_f_n));
	memset(dPeak_scale_f_y,	0,	sizeof(dPeak_scale_f_y));
	memset(dBase_shift_n,	0,	sizeof(dBase_shift_n));
	memset(dBase_shift_y,	0,	sizeof(dBase_shift_y));
	memset(dOffset_n,		0,	sizeof(dOffset_n));
	memset(dOffset_y,		0,	sizeof(dOffset_y));
	memset(dOffset_f_n,		0,	sizeof(dOffset_f_n));
	memset(dOffset_f_y,		0,	sizeof(dOffset_f_y));

	memset(dPeak_err_n,		0,	sizeof(dPeak_err_n));
	memset(dPeak_err_y,		0,	sizeof(dPeak_err_y));
	memset(dPeak_err_f_n,	0,	sizeof(dPeak_err_f_n));
	memset(dPeak_err_f_y,	0,	sizeof(dPeak_err_f_y));
	memset(dBase_err_n,		0,	sizeof(dBase_err_n));
	memset(dBase_err_y,		0,	sizeof(dBase_err_y));
	memset(dOff_err_n,		0,	sizeof(dOff_err_n));
	memset(dOff_err_y,		0,	sizeof(dOff_err_y));
	memset(dOff_err_f_n,	0,	sizeof(dOff_err_f_n));
	memset(dOff_err_f_y,	0,	sizeof(dOff_err_f_y));
	sRise.assign(MAX_CH,	vector<double>());
	sPeakX.assign(MAX_CH,	vector<double>());
	dPeak0.assign(MAX_CH,	vector<double>());
	dPeak2.assign(MAX_CH,	vector<double>());
	sHWHM.assign(MAX_CH,	vector<double>());
	for (auto& it : sRise)	{try {it.reserve(16);} catch (std::exception& e) {throw ProcessorException();}}
	for (auto& it : sPeakX)	{try {it.reserve(16);} catch (std::exception& e) {throw ProcessorException();}}
	for (auto& it : dPeak0)	{try {it.reserve(16);} catch (std::exception& e) {throw ProcessorException();}}
	for (auto& it : dPeak2)	{try {it.reserve(16);} catch (std::exception& e) {throw ProcessorException();}}
	for (auto& it : sHWHM)	{try {it.reserve(16);} catch (std::exception& e) {throw ProcessorException();}}
	pRise = &sRise;
	pPeakX = &sPeakX;
	pPeak0 = &dPeak0;
	pPeak2 = &dPeak2;
	pHWHM = &sHWHM;
}

Processor::~Processor() {
	for (auto& e : event) e.reset();
	for (auto& m : method) m.reset();
	if (g_verbose > 1) cout << "Processor d'tor\n";
	TS = nullptr;
	T0 = nullptr;
	T1 = nullptr;
	if (f) f->Close();
	f = nullptr;
	if (fin.is_open()) fin.close();
	pRise = pPeakX = pHWHM = nullptr;
	pPeak0 = pPeak2 = nullptr;
}

void Processor::BusinessTime() {
	int ev(0), iProgCheck(0), iRate(0), iTimeleft(0);

	unsigned long* ulpTimestamp = (unsigned long*)(buffer.data()+sizeof(long));
	unsigned long ulTSPrev(0);

	TS->Branch("Timestamp", ulpTimestamp, "time_stamp/l");
	TS->Branch("Timestamp_prev", &ulTSPrev, "time_stamp_prev/l");

	steady_clock::time_point t_this, t_that;
	duration<double> t_elapsed;
	iProgCheck = std::max(iNumEvents/(iLevel ? 100 : 10) + 1, (iLevel ? 1000 : 10000)); // too many print statements slows the process

	t_that = steady_clock::now();
	if (g_verbose) {
		cout << "Processing:\n";
		cout << "Completed\tRate (ev/s)\tTime left\n";
	}
	f->cd();
	if (bForceOldFormat) fin.seekg(sizeof_f_header-sizeof(long));
	for (ev = 0; ev < iNumEvents; ev++) {
		fin.read(buffer.data(), iEventsize);

		for (const auto& e : event) e->Analyze();
		if (iLevel) for (const auto& m : method) m->Analyze();

		T0->Fill();
		TS->Fill();
		if (iLevel > 0) T1->Fill();
		ulTSPrev = *ulpTimestamp;
		if ((g_verbose) && (ev % iProgCheck == iProgCheck/2)) { // progress updates
			cout << ev*100l/iNumEvents << "%\t\t";
			t_this = steady_clock::now();
			t_elapsed = duration_cast<duration<double>>(t_this-t_that);
			t_that = steady_clock::now();
			iRate = t_elapsed.count() == 0 ? 9001 : iProgCheck/t_elapsed.count(); // it's OVER 9000!
			cout << iRate << "\t\t";
			iTimeleft = (iNumEvents - ev)/iRate;
			cout << iTimeleft/3600 << "h" << (iTimeleft%3600)/60 << "m" << iTimeleft%60 << "s          \r";
			cout.flush();
		}
		if (s_interrupted) break;
	}
	if (g_verbose) cout << "Processing completed\n";
	T0->AddFriend("TS");
	if (iLevel > 0) {
		T1->AddFriend("TS");
		T1->AddFriend("T0");
		T0->AddFriend("T1");
		T1->Write("",TObject::kOverwrite);
	}
	T0->Write("",TObject::kOverwrite);
	TS->Write("",TObject::kOverwrite);
	if (g_verbose) cout << "Beginning cleanup: ";
	fin.close();
	TS.reset();
	T0.reset();
	T1.reset();
	f->Close();
	f.reset();
	if (g_verbose) cout << " done\n";
	return;
}

vector<void*> Processor::GetAddresses(const int ch, const int level) {
	vector<void*> add(50);
	int i(0);
	if (level == 0) {
		add[i++] = (void*)&bFullWave[ch];
		add[i++] = (void*)&bSaturated[ch];

		add[i++] = (void*)&sDecay[ch];
		add[i++] = (void*)&sRise[ch];
		add[i++] = (void*)&sPeakX[ch];
		add[i++] = (void*)&sTrigger[ch];
		add[i++] = (void*)&sHWHM[ch];
		add[i++] = (void*)&sSaturation[ch];
		add[i++] = (void*)&sNumPeaks[ch];

		add[i++] = (void*)&dBase[ch];
		add[i++] = (void*)&dSigma[ch];
		add[i++] = (void*)&dBaseP[ch];
		add[i++] = (void*)&dBasePS[ch];
		add[i++] = (void*)&dPeak0[ch];
		add[i++] = (void*)&dPeak2[ch];
		add[i++] = (void*)&dFullInt[ch];
		add[i++] = (void*)&dBasePeakP[ch];
		add[i++] = (void*)&dBasePeakN[ch];
	} else if (level == 1) {
		add[i++] = (void*)&bTruncated[ch];

		add[i++] = (void*)&dBaseline[ch];
		add[i++] = (void*)&dBaseSigma[ch];
		add[i++] = (void*)&dBasePost[ch];
		add[i++] = (void*)&dBasePostSigma[ch];
		add[i++] = (void*)&dFastInt[ch];
		add[i++] = (void*)&dSlowInt[ch];

		add[i++] = (void*)&dSample[ch];

		add[i++] = (void*)&dOdd[ch];
		add[i++] = (void*)&dEven[ch];

		add[i++] = (void*)&dLaplaceLow[ch];
		add[i++] = (void*)&dLaplaceHigh[ch];

		add[i++] = (void*)&dXsq_n[ch];
		add[i++] = (void*)&dXsq_y[ch];
		add[i++] = (void*)&dXsq_f_n[ch];
		add[i++] = (void*)&dXsq_f_y[ch];
		add[i++] = (void*)&dPeak_scale_n[ch];
		add[i++] = (void*)&dPeak_scale_y[ch];
		add[i++] = (void*)&dPeak_scale_f_n[ch];
		add[i++] = (void*)&dPeak_scale_f_y[ch];
		add[i++] = (void*)&dBase_shift_n[ch];
		add[i++] = (void*)&dBase_shift_y[ch];
		add[i++] = (void*)&dOffset_n[ch];
		add[i++] = (void*)&dOffset_y[ch];
		add[i++] = (void*)&dOffset_f_n[ch];
		add[i++] = (void*)&dOffset_f_y[ch];

		add[i++] = (void*)&dPeak_err_n[ch];
		add[i++] = (void*)&dPeak_err_y[ch];
		add[i++] = (void*)&dPeak_err_f_n[ch];
		add[i++] = (void*)&dPeak_err_f_y[ch];
		add[i++] = (void*)&dBase_err_n[ch];
		add[i++] = (void*)&dBase_err_y[ch];
		add[i++] = (void*)&dOff_err_n[ch];
		add[i++] = (void*)&dOff_err_y[ch];
		add[i++] = (void*)&dOff_err_f_n[ch];
		add[i++] = (void*)&dOff_err_f_y[ch];
	} else if (level == 2) {
		add[i++] = (void*)&dFastInt[ch];
		add[i++] = (void*)&dSlowInt[ch];
		add[i++] = (void*)&dSample[ch];
		add[i++] = (void*)&dPeak2[ch];
		add[i++] = (void*)&dXsq_n[ch];
		add[i++] = (void*)&dXsq_y[ch];
		add[i++] = (void*)&dPeak_scale_n[ch];
		add[i++] = (void*)&dPeak_scale_y[ch];
		add[i++] = (void*)&dBase_shift_n[ch];
		add[i++] = (void*)&dBase_shift_y[ch];
		add[i++] = (void*)&dOdd[ch];
		add[i++] = (void*)&dEven[ch];
		add[i++] = (void*)&dLaplaceLow[ch];
		add[i++] = (void*)&dLaplaceHigh[ch];
		add[i++] = (void*)&dFullInt[ch];
	}
	return add;
}

void Processor::SetDetectorPositions(const string& in) { // "z0=#,z1=#,z2=#,r0=#,r1=#,r2=#" or some permutation. Accepts z*=# or r*=#
	if (in == "\0") {
		bPositionsSet = false;
		sPositions = "-";
		return;
	}
	sPositions = in;
	if (g_verbose > 1) cout << "Setting detector positions\n";
	int i(0);
	array<unsigned int,2> iCommas{{0,1}};
	vector<string> vGlobs; // for dealing with each each substring
	for (i = 0; i < 6; i++) {
		iCommas[1] = in.find(',',iCommas[0]);
		vGlobs.push_back(in.substr(iCommas[0],iCommas[1]-iCommas[0])); // grabs input between commas
		iCommas[0] = iCommas[1]+1; // shift forward to next substring
		if (iCommas[1] >= in.length()) break;
	}
	for (const auto& it : vGlobs) {
		if ((it[0] == 'z') || (it[0] == 'Z')) {
			if (it[1] == '*') for (i = 0; i < iNchan; i++) fDetectorZ[iChan[i]] = stof(it.substr(it.find('=')+1));
			else fDetectorZ[stoi(it.substr(1,1))] = stof(it.substr(it.find('=')+1));
		} else if ((it[0] == 'r') || (it[0] == 'R')) {
			if (it[1] == '*') for (i = 0; i < iNchan; i++) fDetectorR[iChan[i]] = stof(it.substr(it.find('=')+1));
			else fDetectorR[stoi(it.substr(1,1))] = stof(it.substr(it.find('=')+1));
		}
	}
	bPositionsSet = true;
}

void Processor::Setup(const string& in) { // also opens raw and processed files
	buffer.resize(sizeof_f_header);
	long filesize;
	int ch(-1);
	long lUnixTS(0);
	sName = in;
	if (g_verbose > 1) cout << "Opening files\n";
	sRawDataFile = in;
	sRootFile = in;
	sRootFile.replace(sRootFile.find(".dat"),4,".root");
	sRootFile.replace(sRootFile.find("rawdata"),3,"pro");
	fin.open(sRawDataFile, std::ios::in | std::ios::binary);
	if (!fin.is_open()) {
		cout << "Error: " << sRawDataFile << " not found\n";
		cout << error_message[file_error];
		throw ProcessorException();
	}
	f = unique_ptr<TFile>(new TFile(sRootFile.c_str(), "RECREATE"));
	if (!f->IsOpen()) {
		cout << "Error: could not open " << sRootFile << '\n';
		cout << error_message[file_error];
		throw ProcessorException();
	}

	if (g_verbose) cout << "Parsing file header\n";

	fin.seekg(0, fin.end);
	filesize = fin.tellg();
	fin.seekg(0, fin.beg);
	fin.read(buffer.data(), sizeof_f_header);
	fin.seekg((in[1] >= '6' ? sizeof_f_header : sizeof_f_header-sizeof(long)),fin.beg);
	memcpy(digitizer.cName, buffer.data(), 12); // digitizer name
	memcpy(&usMask, buffer.data() + 12, sizeof(usMask)); // channel mask
	memcpy(&iEventlength, buffer.data() + 14, sizeof(iEventlength)); // Eventlength
	memcpy(&iTrigPost, buffer.data() + 18, sizeof(iTrigPost)); // post-trigger
	memcpy(uiDCOffset, buffer.data() + 22, sizeof(uiDCOffset)); // dc offsets
	memcpy(uiThreshold, buffer.data() + 54, sizeof(uiThreshold)); // trigger thresholds
	if (in[1] >= '6' && !bForceOldFormat) memcpy(&lUnixTS, buffer.data() + 86, sizeof(lUnixTS)); // unix timestamp

	if (string(digitizer.cName) == "DT5730") {
		digitizer.dSamplerate = 5E8;
		digitizer.sResolution = (1 << 14);
		digitizer.dVpp = 2;
		digitizer.iBaselength = 20;
		iFailed = 0;
		digitizer.id = dt5730;
	} else if (string(digitizer.cName) == "DT5751") {
		digitizer.dSamplerate = 1E9;
		digitizer.sResolution = (1 << 10);
		digitizer.dVpp = 1;
		digitizer.iBaselength = 40;
		iFailed = 0;
		digitizer.id = dt5751;
	} else if (string(digitizer.cName) == "DT5751DES") {
		digitizer.dSamplerate = 2E9;
		digitizer.sResolution = (1 << 10);
		digitizer.dVpp = 1;
		digitizer.iBaselength = 80;
		iFailed = 0;
		digitizer.id = dt5751des;
	} else if (string(digitizer.cName) == "V1724") {
		digitizer.dSamplerate = 1E8;
		digitizer.sResolution = (1 << 14);
		digitizer.dVpp = 2.25;
		digitizer.iBaselength = 15;
		iFailed = 0;
		digitizer.id = v1724;
	} else {
		digitizer.dSamplerate = -1;
		digitizer.sResolution = -1;
		digitizer.dVpp = -1;
		digitizer.iBaselength = 1;
		cout << error_message[dig_error];
		digitizer.id = invalid_dig;
		throw ProcessorException();
	}

	digitizer.dVoltsPerBin = digitizer.dVpp/(double)digitizer.sResolution;
	digitizer.dNsPerSample = 1E9/digitizer.dSamplerate;

	iNchan = 0;
	for (auto i = 0; i < MAX_CH; i++) if (usMask & (1<<i)) iChan.push_back(i);
	iNchan = iChan.size();
	iEventsize = sizeof_ev_header + iNchan*iEventlength*sizeof(short); // each sample is size 2
	if (!bForceOldFormat) {
		iNumEvents = (filesize - (in[1] >= '6' ? sizeof_f_header : sizeof_f_header-sizeof(long)))/iEventsize;
	} else {
		iNumEvents = (filesize - sizeof_f_header - sizeof(long))/iEventsize;
	}

	if (g_verbose) cout << "Parsing config file\n";
	string sFilename = sConfigDir + "/config/" + sConfigFileName;
	std::ifstream fconf(sFilename,std::ios::in);
	if (!fconf.is_open()) {
		cout << "Config file " << sFilename << " not found\n";
		cout << error_message[file_error];
		throw ProcessorException();
	} else if (g_verbose > 1) cout << "Opened " << sFilename << '\n';

	buffer[0] = '\0';
	while (!fconf.eof()) {
		fconf.getline(buffer.data(), 64, '\n');
		if (buffer[0] == '#') continue;
		if (string(buffer.data()+10) ==  string(digitizer.cName)) {
			fconf.getline(buffer.data(), 64, '\n');
			while (strstr(buffer.data(), "END") == NULL) {
				if (strstr(buffer.data(), "CHANNEL") != NULL) { // loads processing parameters
					ch = atoi(&buffer[8]);
					if ((ch >= MAX_CH) || (ch < 0)) {cout << error_message[config_file_error]; throw ProcessorException();}
					fconf.getline(buffer.data(), 64, '\n');
					sscanf(buffer.data(), "SLOW %i FAST %i PGA %i GAIN_N %f GAIN_Y %f", &iSlowTime[ch], &iFastTime[ch], &iPGASamples[ch], &fGain[ch][0], &fGain[ch][1]);
				}
				fconf.getline(buffer.data(), 64, '\n');
			} // end of while
		} // end of if dig
	} // end of file
	fconf.close();

	int iDateNow(0), iTimeNow(0);
	string sTime;
	time_t t_rawtime;
	tm* t_today;
	time(&t_rawtime);
	t_today = localtime(&t_rawtime); // timestamp for processing
	iTimeNow = (t_today->tm_hour)*100 + t_today->tm_min; // hhmm
	iDateNow = (t_today->tm_year-100)*10000 + (t_today->tm_mon+1)*100 + t_today->tm_mday; // yymmdd
	//sprintf(cTime.data(), "%i_%i",iDateNow,iTimeNow);
	sTime = std::to_string(iDateNow) + "_" + std::to_string(iTimeNow);
	lUnixTS = in[1] >= '6' && !bForceOldFormat ? lUnixTS : UnixConverter(in.substr(in.find(".dat")-11,11));

	switch (digitizer.id) {
		case dt5751 :
			iXSQ_ndf = std::min(iEventlength, 450)-3;
			break;
		case dt5751des :
			iXSQ_ndf = std::min(iEventlength, 899)-3;
			break;
		case dt5730 :
			iXSQ_ndf = std::min(iEventlength, 225)-3;
			break;
		case v1724 :
		case invalid_dig :
		default :
			iXSQ_ndf = -1;
			break;
	}

	cout << "Creating config trees\n";
	try {T0.reset(new TTree("TI","Info"));}
	catch (std::bad_alloc& ba) {
		cout << "Could not create info tree\n";
		cout << error_message[alloc_error];
		throw ProcessorException();
	}
	if (string(cSource).find("252") != string::npos) { // Cf-252 (NG requirement relaxed)
		if (!bPositionsSet) {
			if (g_verbose == 0) { // assumed running in batch mode
				cout << "No detector positions specified\n";
				throw ProcessorException();
			}
			cout << "Enter detector positions:\n";
			for (auto i = 0; i < iNchan; i++) {
				cout << "Detector " << iChan[i] << " z: "; std::cin >> fDetectorZ[i];
				cout << "Detector " << iChan[i] << " r: "; std::cin >> fDetectorR[i];
			}	}
	}
	if (string(cSource).find("NG") != string::npos) {
		if ((fHV == 0) || (fCurrent == 0)) {
			if (g_verbose == 0) {
				cout << "No NG setpoint values\n";
				throw ProcessorException();
			}
			cout << "Enter NG setpoints\n";
			cout << "HV: "; std::cin >> fHV;
			cout << "Current: "; std::cin >> fCurrent;
		}
	}

	T0->Branch("Digitizer",			digitizer.cName,		"name[12]/B");
	T0->Branch("Source",			cSource ,			"source[12]/B");
	T0->Branch("ChannelMask",		&usMask,				"mask/s"); // general info on data run
	T0->Branch("TriggerThreshold",	uiThreshold ,		"threshold[8]/i");
	T0->Branch("DC_offset",			uiDCOffset ,		"dc_off[8]/i");
	T0->Branch("Posttrigger",		&iTrigPost,				"tri_post/I");
	T0->Branch("Eventlength",		&iEventlength,			"ev_len/I");
	T0->Branch("UnixTS",			&lUnixTS,				"unixts/L");
	T0->Branch("Level",				&iLevel,				"level/I");
	T0->Branch("Chisquared_NDF",	&iXSQ_ndf,				"ndf/I");
//	T0->Branch("Time",				sTime ,			"time[12]/B");
	T0->Branch("PGA_samples",		iPGASamples ,		"pga[8]/I");
	T0->Branch("Fast_window",		iFastTime ,		"fast[8]/I");
	T0->Branch("Slow_window",		iSlowTime ,		"slow[8]/I");
	T0->Branch("Detector_pos_z",	fDetectorZ ,		"z_pos[3]/F");
	T0->Branch("Detector_pos_r",	fDetectorR ,		"r_pos[3]/F");
	T0->Branch("NG_HV_setpoint",	&fHV,					"hv/F");
	T0->Branch("NG_Current_setpoint",&fCurrent,				"current/F");
	if (iAverage != 0) T0->Branch("Moving_average", &iAverage, "average/I");

	T0->Fill();
	T0->Write();
	T0.reset();

	if (g_verbose > 1) cout << "Alloc'ing trees\n";
	try {buffer.resize(iEventsize);}
	catch (std::exception& e) {
		cout << error_message[alloc_error] << "Buffer\n";
		cout << e.what() << '\n';
		throw ProcessorException();
	}
	memset(buffer.data(), 0, iEventsize);
	unsigned short* uspTrace = (unsigned short*)(buffer.data() + sizeof_ev_header);

	try {
		TS = unique_ptr<TTree>(new TTree("TS","Timestamps"));
		T0 = unique_ptr<TTree>(new TTree("T0","Event"));
		if (iLevel > 0) T1 = unique_ptr<TTree>(new TTree("T1","Method"));
	} catch (std::bad_alloc& ba) {
		cout << error_message[alloc_error] << "Trees\n";
		throw ProcessorException();
	}

	T0->Branch("FullWaveform",	bFullWave ,	"fullwave[8]/O");
	T0->Branch("Saturated",		bSaturated , 	"sat[8]/O");

	T0->Branch("Decaytime",		sDecay ,		"decay[8]/S");
	T0->Branch("Risetime",		"vector<vector<double>>", &pRise);
	T0->Branch("Peakx",			"vector<vector<double>>", &pPeakX);
	T0->Branch("Trigger",		sTrigger ,	"trig[8]/S");
	T0->Branch("HWHM",			"vector<vector<double>>", &pHWHM);
	T0->Branch("SatDur",		sSaturation ,"satdur[8]/S");

	T0->Branch("Base",			dBase ,		"base[8]/D");
	T0->Branch("Sigma",			dSigma ,		"sigma[8]/D");
	T0->Branch("BaseP",			dBaseP ,		"basep[8]/D");
	T0->Branch("BasePS",		dBasePS ,		"baseps[8]/D");
	T0->Branch("BasePkP",		dBasePeakP ,	"basepkp[8]/D");
	T0->Branch("BasePkN",		dBasePeakN ,	"basepkn[8]/D");
	T0->Branch("Peakheight0",	"vector<vector<double>>", &pPeak0);
	T0->Branch("Peakheight2",	"vector<vector<double>>", &pPeak2);
	T0->Branch("Integral",		dFullInt ,	"integral[8]/D");

	if (iLevel > 0) {
		T1->Branch("Truncated",		bTruncated , 		"trunc[8]/O");

		T1->Branch("Baseline",		dBaseline ,		"baseline[8]/D");
		T1->Branch("BaseSigma",		dBaseSigma ,		"basesigma[8]/D");
		T1->Branch("BasePost",		dBasePost ,		"basepost[8]/D");
		T1->Branch("BasePostSigma",	dBasePostSigma ,	"basepostsigma[8]/D");
		T1->Branch("FastInt",		dFastInt ,		"fastint[8]/D");
		T1->Branch("SlowInt",		dSlowInt ,		"slowint[8]/D");

		T1->Branch("Sample",		dSample ,			"sample[8]/D");

		T1->Branch("Odd",			dOdd ,			"odd[8]/D");
		T1->Branch("Even",			dEven ,			"even[8]/D");

		T1->Branch("LapLow",		dLaplaceLow ,		"laplow[8]/D");
		T1->Branch("LapHi",			dLaplaceHigh ,	"laphi[8]/D");

		T1->Branch("Xsq_n",			dXsq_n ,			"xsqn[4]/D");
		T1->Branch("Xsq_y",			dXsq_y ,			"xsqy[4]/D");
		T1->Branch("Xsq_n_f",		dXsq_f_n ,		"xsqnf[4]/D");
		T1->Branch("Xsq_y_f",		dXsq_f_y ,		"xsqyf[4]/D");
		T1->Branch("Peakscale_n",	dPeak_scale_n ,	"pkn[4]/D");
		T1->Branch("Peakscale_y",	dPeak_scale_y ,	"pky[4]/D");
		T1->Branch("Peakscale_n_f",	dPeak_scale_f_n ,	"pknf[4]/D");
		T1->Branch("Peakscale_y_f",	dPeak_scale_f_y ,	"pkyf[4]/D");
		T1->Branch("Base_shift_n",	dBase_shift_n ,	"bshn[4]/D");
		T1->Branch("Base_shift_y",	dBase_shift_y ,	"bshy[4]/D");
		T1->Branch("Offset_n",		dOffset_n ,		"offn[4]/D");
		T1->Branch("Offset_y",		dOffset_y ,		"offy[4]/D");
		T1->Branch("Offset_n_f",	dOffset_f_n ,		"offnf[4]/D");
		T1->Branch("Offset_y_f",	dOffset_f_y ,		"offyf[4]/D");

		T1->Branch("Peak_err_n",	dPeak_err_n ,		"pkerrn[4]/D");
		T1->Branch("Peak_err_y",	dPeak_err_y ,		"pkerry[4]/D");
		T1->Branch("Peak_err_n_f",	dPeak_err_f_n ,	"pkerrnf[4]/D");
		T1->Branch("Peak_err_y_f",	dPeak_err_f_y ,	"pkerryf[4]/D");
		T1->Branch("Base_err_n",	dBase_err_n ,		"baerrn[4]/D");
		T1->Branch("Base_err_y",	dBase_err_y ,		"baerry[4]/D");
		T1->Branch("Off_err_n",		dOff_err_n ,		"oferrn[4]/D");
		T1->Branch("Off_err_y",		dOff_err_y ,		"oferry[4]/D");
		T1->Branch("Off_err_n_f",	dOff_err_f_n ,	"oferrnf[4]/D");
		T1->Branch("Off_err_y_f",	dOff_err_f_y ,	"oferryf[4]/D");
	}
	cout << "Processing level " << iLevel << '\n';

	for (auto& ch : iChan) { // initializing all classes needed
		if (g_verbose > 1) cout << "CH" << ch << '\n';
		try {
			dTrace.push_back(vector<double>(iEventlength-iAverage));
			event.push_back(unique_ptr<Event>(new Event(iEventlength, digitizer.iBaselength, iAverage, uiThreshold[ch], ch, uspTrace + ch*iEventlength, dTrace.back().data())));
		} catch (std::bad_alloc& ba) {
			cout << error_message[alloc_error] << "Event\n";
			throw ProcessorException();
		}
		if (event.back()->Failed()) {
			cout << error_message[method_error] << "Event\n";
			throw ProcessorException();
		}
		event.back()->SetScales(digitizer.dVoltsPerBin, digitizer.dNsPerSample);
		event.back()->SetAddresses(GetAddresses(ch,0));
		if (iLevel > 0) {
			try {
				method.push_back(unique_ptr<Method>(new Method(event.back()->Length(), iFastTime[ch], iSlowTime[ch], iPGASamples[ch], fGain[ch], digitizer.dNsPerSample, digitizer.dVoltsPerBin, event.back(), sConfigDir)));
			} catch (std::bad_alloc& ba) {
				cout << error_message[alloc_error] << "Method\n";
				throw ProcessorException();
			}
			if (method.back()->Failed()) {
				cout << error_message[method_error] << "Method\n";
				throw ProcessorException();
			}
			method.back()->SetAddresses(GetAddresses(ch,1));
			method.back()->SetDCOffset(digitizer.sResolution, uiDCOffset[ch]);
		}
	} // ch loop
}
