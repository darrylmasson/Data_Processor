#include "Processor.h"
#include <ctime>
#include <pthread.h>
#include <ratio>
#include <chrono>

using namespace std::chrono;

const float cfMethodVersions[NUM_METHODS] = {
	CCM::sfVersion,
	DFT::sfVersion,
	XSQ_TF1::sfVersion,
	LAP::sfVersion,
	XSQ_NEW::sfVersion
};

auto (*root_init[])(TTree*) -> void = {
	CCM::root_init,
	DFT::root_init,
	XSQ_TF1::root_init,
	LAP::root_init,
	XSQ_NEW::root_init
};

auto (*root_fill[])() -> void = {
	CCM::root_fill,
	DFT::root_fill,
	XSQ_TF1::root_fill,
	LAP::root_fill,
	XSQ_NEW::root_fill
};

auto (*root_write[])() -> void = {
	CCM::root_write,
	DFT::root_write,
	XSQ_TF1::root_write,
	LAP::root_write,
	XSQ_NEW::root_write
};

auto (*root_deinit[])() -> void = {
	CCM::root_deinit,
	DFT::root_deinit,
	XSQ_TF1::root_deinit,
	LAP::root_deinit,
	XSQ_NEW::root_deinit
};

void* Process(void* arg) {
	thread_data_t* input = (thread_data_t*)arg;
	input->event->Analyze();
	for (int i = 0; i < NUM_METHODS; i++) if (input->cbpActivated[i]) input->methods[i]->Analyze();
	return nullptr;
}

Processor::Processor() {
	if (g_verbose) cout << "Processor c'tor\n";
	iSpecial = -1;
	iAverage = 0;
	strcpy(cMethodNames[0], "CCM_PGA");
	strcpy(cMethodNames[1], "FOURIER");
	strcpy(cMethodNames[2], "XSQ_TF1");
	strcpy(cMethodNames[3], "LAPLACE");
	strcpy(cMethodNames[4], "XSQ_NEW");

	sConfigFileName = "\0";
	sRawDataFile = "\0";
	sRootFile = "\0";

	memset(bMethodActive, 0, sizeof(bMethodActive));
	memset(bMethodDone, 0, sizeof(bMethodDone));
	bRecordTimestamps = false;

	cDigName[0] = '\0';
	cSource[0] = '\0';

	usMask = 0;

	memset(iChan, 0, sizeof(iChan));
	iEventlength = 0;
	iEventsize = 0;
	iFailed = 0;
	memset(iFastTime, 0, sizeof(iFastTime));
	iNchan = 0;
	iNumEvents = 0;
	memset(iPGASamples, 0, sizeof(iPGASamples));
	memset(iSlowTime, 0, sizeof(iSlowTime));
	iTrigPost = 0;
	iXSQ_ndf = 0;

	memset(uiDCOffset, 0, sizeof(uiDCOffset));
	memset(uiThreshold, 0, sizeof(uiThreshold));

	memset(fGain, 0, sizeof(fGain));
	memset(fDetectorZ, 0, sizeof(fDetectorZ));
	memset(fDetectorR, 0, sizeof(fDetectorR));

	memset(td, 0, sizeof(td));
}

Processor::~Processor() {
	if (g_verbose) cout << "Processor d'tor\n";
	for (auto ch = 0; ch < MAX_CH; ch++) {
		td[ch].event.reset();
		for (auto m = 0; m < NUM_METHODS; m++) td[ch].methods[m].reset();
		td[ch].cbpActivated = nullptr;
	}
	digitizer.reset();
	f.reset();
	tree.reset();
	buffer.reset();
	if (fin.is_open()) fin.close();
}

void Processor::BusinessTime() {
	int ch(0), ev(0), m(0), rc(0), iProgCheck(0), iRate(0), iTimeleft(0), iLivetime(0);
	pthread_t threads[MAX_CH];
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	void* status = nullptr;

	unsigned long* ulpTimestamp = (unsigned long*)(buffer.get() + sizeof(long));
	unsigned long ulTSFirst(0), ulTSLast(0), ulTSPrev(0);
	if (bRecordTimestamps) {
		tree->Branch("Timestamp", &ulpTimestamp[0], "time_stamp/l");
		tree->Branch("Timestamp_prev", &ulTSPrev, "time_stamp_prev/l");
	}

	steady_clock::time_point t_this, t_that;
	duration<double> t_elapsed;
	iProgCheck = max(iNumEvents/100 + 1, 10000); // too many print statements slows the process, every 1% or 10000 events

	if (memchr(bMethodActive, 1, NUM_METHODS) == nullptr) {
		cout << "No processing method activated\n";
		return;
	} else {
		cout << "Active methods: ";
		for (m = 0; m < NUM_METHODS; m++) if (bMethodActive[m]) cout << cMethodNames[m] << " ";
		cout << '\n';
	}

	t_that = steady_clock::now();
	cout << "Processing:\n";
	cout << "Completed\tRate (ev/s)\tTime left\n";
	fin.seekg(sizeof_f_header, fin.beg);
	f->cd();
	for (ev = 0; ev < iNumEvents; ev++) {
		fin.read(buffer.get(), iEventsize);
		if (bMethodActive[XSQ_TF1_t]) for (ch = 0; ch < iNchan; ch++) { // TF1 isn't thread-friendly
			td[ch].event->Analyze();
			for (m = 0; m < NUM_METHODS; m++) if (bMethodActive[m]) td[ch].methods[m]->Analyze();
		} else { // but Newton's method is
			for (ch = 0; ch < iNchan; ch++) if ( (rc = pthread_create(&threads[ch], &attr, Process, (void*)&td[ch])) ) {iFailed |= (1 << thread_error); return;}
			for (ch = 0; ch < iNchan; ch++) if ( (rc = pthread_join(threads[ch], &status)) ) {iFailed |= (1 << thread_error); return;}
		}
		for (m = 0; m < NUM_METHODS; m++) if (bMethodActive[m]) root_fill[m]();
		if (bRecordTimestamps) {
			tree->Fill();
			ulTSPrev = ulpTimestamp[0];
		}
		if (ev % iProgCheck == iProgCheck-1) { // progress updates
			cout << ev*100l/iNumEvents << "%\t\t";
			t_this = steady_clock::now();
			t_elapsed = duration_cast<duration<double>>(t_this-t_that);
			t_that = steady_clock::now();
			iRate = t_elapsed.count() == 0 ? 9001 : iProgCheck/t_elapsed.count(); // it's OVER 9000!
			cout << iRate << "\t\t";
			iTimeleft = (iNumEvents - ev)/iRate;
			if (iTimeleft > (1 << 12)) cout << iTimeleft/3600 << "h" << (iTimeleft%3600)/60 << "m\n";
			else if (iTimeleft > (1 << 7)) cout << iTimeleft/60 << "m" << iTimeleft%60 << "s\n";
			else cout << iTimeleft << "s\n";
		}
		if (ev == 0) ulTSFirst = ulpTimestamp[0];
	}
	tree->Write();
	cout << "Processing completed\n";
	ulTSLast = ulpTimestamp[0];
	iLivetime = (ulTSLast - ulTSFirst)/125e6;
	cout << "Acquisition livetime: " << iLivetime << "s\nBeginning cleanup: ";
	fin.close();
	for (m = 0; m < NUM_METHODS; m++) if (bMethodActive[m]) {
		root_write[m]();
		root_deinit[m]();
	}
	buffer.reset();
	tree.reset();
	auto i = system(("chmod g+w " + sRootFile).c_str());
	i++; // to keep g++ happy about unused variables
	cout << " done\n";
	return;
}

void Processor::ClassAlloc() {
	if (g_verbose) cout << "Alloc'ing\n";
	try {buffer = unique_ptr<char[]>(new char[iEventsize]);}
	catch (bad_alloc& ba) {
		iFailed |= (1 << alloc_error);
		throw ProcessorException();
	}
	unsigned short* uspTrace = (unsigned short*)(buffer.get() + sizeof_ev_header);

	f->cd();

	for (auto ch = 0; ch < iNchan; ch++) { // initializing all classes needed
		if (g_verbose) cout << "CH" << ch << '\n';
		try {
			if (iAverage == 0) td[ch].event = shared_ptr<Event>(new Event(iEventlength, digitizer));
			else td[ch].event = shared_ptr<Event_ave>(new Event_ave(iEventlength, digitizer));
		} catch (bad_alloc& ba) {
			iFailed |= (1 << alloc_error);
			throw ProcessorException();
		}
		td[ch].event->SetAverage(iAverage);
		td[ch].event->SetDCOffset(digitizer, uiDCOffset[iChan[ch]]);
		td[ch].event->SetThreshold(uiThreshold[iChan[ch]]);
		td[ch].event->SetTrace(uspTrace + ch*iEventlength);
		if (td[ch].event->Failed()) {
			iFailed |= (1 << method_error);
			throw ProcessorException();
		}
		if (bMethodActive[CCM_t]) {
			try {td[ch].methods[CCM_t] = shared_ptr<Method>(new CCM(iChan[ch],Event::Length(),digitizer));}
			catch (bad_alloc& ba) {
				iFailed |= (1 << alloc_error);
				bMethodActive[CCM_t] = false;
			}
			td[ch].methods[CCM_t]->SetParameters(&iFastTime[iChan[ch]], 0, digitizer);
			td[ch].methods[CCM_t]->SetParameters(&iSlowTime[iChan[ch]], 1, digitizer);
			td[ch].methods[CCM_t]->SetParameters(&iPGASamples[iChan[ch]], 2, digitizer);
			td[ch].methods[CCM_t]->SetEvent(td[ch].event);
			if (td[ch].methods[CCM_t]->Failed()) {
				iFailed |= (1 << method_error);
				bMethodActive[CCM_t] = false;
			}
		}
		if (bMethodActive[DFT_t]) {
			try {td[ch].methods[DFT_t] = shared_ptr<Method>(new DFT(iChan[ch], Event::Length(), digitizer));}
			catch (bad_alloc& ba) {
				iFailed |= (1 << alloc_error);
				bMethodActive[DFT_t] = false;
			}
			if (td[ch].methods[DFT_t]->Failed()) {
				iFailed |= (1 << method_error);
				bMethodActive[DFT_t] = false;
			}
			td[ch].methods[DFT_t]->SetEvent(td[ch].event);
		}
		if (bMethodActive[XSQ_TF1_t]) {
			try {td[ch].methods[XSQ_TF1_t] = shared_ptr<Method>(new XSQ_TF1(iChan[ch], Event::Length(), digitizer));}
			catch (bad_alloc& ba) {
				iFailed |= (1 << alloc_error);
				bMethodActive[XSQ_TF1_t] = false;
			}
			td[ch].methods[XSQ_TF1_t]->SetParameters(&fGain[iChan[ch]][0], 0, digitizer);
			td[ch].methods[XSQ_TF1_t]->SetParameters(&fGain[iChan[ch]][1], 1, digitizer);
			td[ch].methods[XSQ_TF1_t]->SetEvent(td[ch].event);
			if (td[ch].methods[XSQ_TF1_t]->Failed()) {
				iFailed |= (1 << method_error);
				bMethodActive[XSQ_TF1_t] = false;
			}
		}
		if (bMethodActive[LAP_t]) {
			try {td[ch].methods[LAP_t] = shared_ptr<Method>(new LAP(iChan[ch], Event::Length(), digitizer));}
			catch (bad_alloc& ba) {
				iFailed |= (1 << alloc_error);
				bMethodActive[LAP_t] = false;
			}
			td[ch].methods[LAP_t]->SetEvent(td[ch].event);
			if (td[ch].methods[LAP_t]->Failed()) {
				iFailed |= (1 << method_error);
				bMethodActive[LAP_t] = false;
			}
		}
		if (bMethodActive[XSQ_NEW_t]) {
			try {td[ch].methods[XSQ_NEW_t] = shared_ptr<Method>(new XSQ_NEW(iChan[ch], Event::Length(), digitizer));}
			catch (bad_alloc& ba) {
				iFailed |= (1 << alloc_error);
				bMethodActive[XSQ_NEW_t] = false;
			}
			td[ch].methods[XSQ_NEW_t]->SetParameters(&fGain[iChan[ch]][0], 0, digitizer);
			td[ch].methods[XSQ_NEW_t]->SetParameters(&fGain[iChan[ch]][1], 1, digitizer);
			td[ch].methods[XSQ_NEW_t]->SetEvent(td[ch].event);
			if (td[ch].methods[XSQ_NEW_t]->Failed()) {
				iFailed |= (1 << method_error);
				bMethodActive[XSQ_NEW_t] = false;
			}
		}
		td[ch].cbpActivated = bMethodActive;
		if (g_verbose) cout << '\n';
	}

	for (auto m = 0; m < NUM_METHODS; m++) { // setting up trees
		sprintf(cTreename[m], "T%i", m);
		if (bMethodActive[m]) {
			try {tree = unique_ptr<TTree>(new TTree(cTreename[m], cMethodNames[m]));}
			catch (bad_alloc& ba) {iFailed |= (1 << alloc_error); bMethodActive[m] = false;}
			if (tree->IsZombie()) {iFailed |= (1 << root_error); bMethodActive[m] = false;}
			root_init[m](tree.release());
	}	}

	if (bRecordTimestamps) {
		try {tree = unique_ptr<TTree>(new TTree("TS","Timestamps"));}
		catch (bad_alloc& ba) {
			iFailed |= (1 << alloc_error);
			return;
		} if (tree->IsZombie()) {
			iFailed |= (1 << root_error);
			return;
		}
	}

	if (iFailed) throw ProcessorException();
}

void Processor::ConfigTrees() {
	bool bUpdate(false), bChecked[NUM_METHODS];
	char overwrite(0), cMethodName[12];
	int iMethodID(0), iDateNow(0), iTimeNow(0), iPGACheck[MAX_CH], iSlowCheck[MAX_CH], iFastCheck[MAX_CH];
	float fVersion(0);
	unique_ptr<TTree> tc = nullptr;
	memset(bChecked, 0, sizeof(bChecked));
	time_t t_rawtime;
	tm* t_today;
	time(&t_rawtime);
	t_today = localtime(&t_rawtime); // timestamp for processing
	iTimeNow = (t_today->tm_hour)*100 + t_today->tm_min; // hhmm
	iDateNow = (t_today->tm_year-100)*10000 + (t_today->tm_mon+1)*100 + t_today->tm_mday; // yymmdd
	switch (digitizer->ID()) {
		case dt5751 :
			if (digitizer->Special() == 0) iXSQ_ndf = min(iEventlength/2, 225)-3; // length - number of free parameters
			else iXSQ_ndf = min(iEventlength, 450)-3;
			break;
		case dt5751des :
			iXSQ_ndf = min(iEventlength, 899)-3;
			break;
		case dt5730 :
			iXSQ_ndf = min(iEventlength, 225)-3;
			break;
		case v1724 :
		case invalid_dig :
		default :
			iXSQ_ndf = -1;
			break;
	}
	tree = unique_ptr<TTree>((TTree*)f->Get("Tx"));
	if (tree) { // check for old versions
		cout << "This file has been processed by an old version of NG_DP. Reprocessing will require removal of previous results. Continue <y|n>? ";
		cin >> overwrite;
		if (overwrite == 'y') f->Delete("*;*");
		else {
			iFailed |= (1 << root_error);
			throw ProcessorException();
		}
	}
	tree.reset((TTree*)f->Get("TV"));
	if (tree) { // already processed, checking versions
		tc.reset((TTree*)f->Get("TC"));
		bRecordTimestamps = false;
		cout << "Config trees already exist, checking versions\n";
		tree->SetBranchAddress("MethodID", &iMethodID);
		tree->SetBranchAddress("Version", &fVersion);
		for (auto i = tree->GetEntries()-1; i >= 0; i--) { // most recent entries will be last in the tree
			tree->GetEntry(i);
			if (bChecked[iMethodID] || !bMethodActive[iMethodID]) continue;
			if (g_verbose) cout << "Checking " << cMethodNames[iMethodID] << '\n';
			bChecked[iMethodID] = true;
			bUpdate = false;
			if (iMethodID == CCM_t) {
				tc->SetBranchAddress("PGA_samples", iPGACheck);
				tc->SetBranchAddress("Fast_window", iFastCheck);
				tc->SetBranchAddress("Slow_window", iSlowCheck);
				tc->GetEntry(tc->GetEntries()-1); // most recent entry is the one we want
				bUpdate |= (memcmp(iPGACheck, iPGASamples, sizeof(iPGACheck)) != 0);
				bUpdate |= (memcmp(iFastCheck, iFastTime, sizeof(iFastCheck)) != 0); // checks CCM parameters
				bUpdate |= (memcmp(iSlowCheck, iSlowTime, sizeof(iSlowCheck)) != 0);
				tc.reset();
			}
			if ((fVersion < cfMethodVersions[iMethodID]) || bUpdate) {
				cout << cMethodNames[iMethodID] << " will be reprocessed\n";
			} else {
				cout << cMethodNames[iMethodID] << " up to date, reprocess anyway <y|n>? ";
				cin >> overwrite;
				if (overwrite == 'y') {
					cout << "Reprocessing " << cMethodNames[iMethodID] << '\n';
				} else {
					bMethodActive[iMethodID] = false;
				}
			}
		} // end of checking
		for (auto i = 0; i < NUM_METHODS; i++) if (bMethodActive[i]) {
			iMethodID = i;
			fVersion = cfMethodVersions[i];
			tree->Fill();
			if (i == CCM_t) {
				memcpy(iPGACheck, iPGASamples, sizeof(iPGACheck));
				memcpy(iFastCheck, iFastTime, sizeof(iFastCheck));
				memcpy(iSlowCheck, iSlowTime, sizeof(iSlowCheck));
				tc->Fill();
			}
		}
		tree->Write("",TObject::kOverwrite);
		if (bMethodActive[CCM_t]) tc->Write("",TObject::kOverwrite);
	} else { // not already processed
		cout << "Creating config trees\n";
		bRecordTimestamps = true;
		try {tree.reset(new TTree("TI","Info"));}
		catch (bad_alloc& ba) {
			cout << "Could not create info tree\n";
			iFailed |= (1 << alloc_error);
			throw ProcessorException();
		}
		strcpy(cBuildID, sBuildID.c_str());
		tree->Branch("Digitizer", cDigName, "name[12]/B");
		tree->Branch("Source", cSource, "source[12]/B");
		tree->Branch("ChannelMask", &usMask, "mask/s"); // general info on data run
		tree->Branch("TriggerThreshold", uiThreshold, "threshold[8]/i");
		tree->Branch("DC_offset", uiDCOffset, "dc_off[8]/i"); // the numbers in this tree
		tree->Branch("Posttrigger", &iTrigPost, "tri_post/I"); // don't change, so it's only
		tree->Branch("Eventlength", &iEventlength, "ev_len/I"); // written out once
		tree->Branch("Chisquared_NDF", &iXSQ_ndf, "ndf/I");
		tree->Branch("Fitter", &iXSQ_fitter, "fitter/I");
		tree->Branch("Build_ID", cBuildID, "buildid[21]/B");
		if (iSpecial != -1) tree->Branch("Special", &iSpecial, "special/I");
		if (iAverage != 0) tree->Branch("Moving_average", &iAverage, "average/I");
		if ((strcmp(cSource, "NG") == 0) || (strstr(cSource, "252") != nullptr)) { // NG or Cf-252
			if (!bPositionsSet) {
				cout << "Enter detector positions:\n";
				for (auto i = 0; i < 3; i++) {
					cout << "Detector " << i << " z: "; cin >> fDetectorZ[i];
					cout << "Detector " << i << " r: "; cin >> fDetectorR[i];
			}	}
			tree->Branch("Detector_position_z", fDetectorZ, "z_pos[3]/F");
			tree->Branch("Detector_position_r", fDetectorR, "r_pos[3]/F");
		}
		tree->Fill();
		tree->Write();

		try {tree.reset(new TTree("TV","Versions"));}
		catch (bad_alloc& ba) {
			cout << "Could not create version tree\n";
			iFailed |= (1 << alloc_error);
			throw ProcessorException();
		}
		tree->Branch("MethodName", cMethodName, "codename[12]/B");
		tree->Branch("MethodID", &iMethodID, "codeid/I");
		tree->Branch("Date", &iDateNow, "date/I"); // this tree holds version info
		tree->Branch("Time", &iTimeNow, "time/I");
		tree->Branch("Version", &fVersion, "version/F");
		for (auto i = 0; i < NUM_METHODS; i++) if (bMethodActive[i]) {
			strcpy(cMethodName, cMethodNames[i]);
			fVersion = cfMethodVersions[i];
			iMethodID = i;
			tree->Fill();
		}
		tree->Write();

		try {tree.reset(new TTree("TC","CCM_info"));}
		catch (bad_alloc& ba) {
			cout << "Could not create CCM info tree\n";
			iFailed |= (1 << alloc_error);
			throw ProcessorException();
		}
		tree->Branch("PGA_samples", iPGASamples, "pga[8]/I");
		tree->Branch("Fast_window", iFastTime, "fast[8]/I"); // this tree holds configuration parameters for CCM
		tree->Branch("Slow_window", iSlowTime, "slow[8]/I");
		tree->Fill();
		tree->Write("",TObject::kOverwrite);

		tree.reset();
	}
	tree.reset();
	tc.reset();
}

void Processor::FriendshipIsMagic() {
	if (g_verbose) cout << "Making friends: ";
	bool bCuts(false), bTree[NUM_METHODS];
	int m(0), i(0);
	string sCutsFile = sRootFile;
	sCutsFile.insert(sRootFile.find('.'),"Cuts");
	if (!bRecordTimestamps) { // check for existence of cuts file - only if reprocessing
		fin.open(sCutsFile.c_str(), ios::in);
		if (fin.is_open()) {
			fin.close();
			bCuts = true;
		}
	}
	f->cd();
	for (m = 0; m < NUM_METHODS; m++) {
		tree = unique_ptr<TTree>((TTree*)f->Get(cTreename[m]));
		if (tree) bTree[m] = true;
		else bTree[m] = false;
	}
	for (m = 0; m < NUM_METHODS; m++) {
		if (bTree[m]) { // all existing trees
			if (g_verbose) cout << cTreename[m] << "\n";
			tree = unique_ptr<TTree>((TTree*)f->Get(cTreename[m]));
			tree->AddFriend("TS");
			for (i = 1; i < NUM_METHODS; i++) if (bTree[(m+i)%NUM_METHODS]) tree->AddFriend(cTreename[(m+i)%NUM_METHODS]);
			if (bCuts) tree->AddFriend("Tcuts",sCutsFile.c_str());
			tree->Write("",TObject::kOverwrite);
			tree.reset();
		}
	}
	tree.reset();
	return;
}

void Processor::ParseFileHeader() {
	if (g_verbose) cout << "Parsing file header\n";
	char cBuffer[sizeof_f_header];
	fin.seekg(0, fin.end);
	long filesize = fin.tellg();
	fin.seekg(0, fin.beg);
	fin.read(cBuffer, sizeof_f_header);

	strncpy(cDigName, cBuffer, sizeof(cDigName)); // digitizer name

	memcpy(&usMask, cBuffer + 12, sizeof(usMask)); // channel mask

	memcpy(&iEventlength, cBuffer + 14, sizeof(iEventlength)); // eventlength

	memcpy(&iTrigPost, cBuffer + 18, sizeof(iTrigPost)); // post-trigger

	memcpy(uiDCOffset, cBuffer + 22, sizeof(uiDCOffset)); // dc offsets

	memcpy(uiThreshold, cBuffer + 54, sizeof(uiThreshold)); // trigger thresholds

	iNchan = 0;
	for (auto i = 0; i < MAX_CH; i++) if (usMask & (1<<i)) iChan[iNchan++] = i;
	iEventsize = sizeof_ev_header + iNchan*iEventlength*sizeof(short); // each sample is size 2
	iNumEvents = (filesize - sizeof_f_header)/iEventsize;
}

void Processor::ParseConfigFile() {
	if (g_verbose) cout << "Parsing config file\n";
	string sFilename = sWorkingDir + "/Data_Processor/config/" + sConfigFileName;
	ifstream fconf(sFilename.c_str(),ios::in);
	if (!fconf.is_open()) {
		cout << "Config file " << sFilename << " not found\n";
		iFailed |= (1 << file_error);
		throw ProcessorException();
	} else if (g_verbose) cout << "Opened " << sFilename << '\n';

	char temp[32] = {'\0'}, cBuffer[64] = {'\0'};
	int ch(-1), code(0);
	while (!fconf.eof()) {
		fconf.getline(cBuffer, 64, '\n');
		if (cBuffer[0] == '#') continue;
		if (strcmp(cBuffer, "METHODS") == 0) { // checks for active processing methods
			fconf.getline(cBuffer, 64, '\n');
			while (strstr(cBuffer, "END") == NULL) {
				sscanf(cBuffer, "%s %i", temp, &code);
				for (int i = 0; i < NUM_METHODS; i++) if (strcmp(temp, cMethodNames[i]) == 0) bMethodActive[i] = code;
				fconf.getline(cBuffer, 64, '\n');
			} // end of while
		} // end of if method
		if (strcmp(cBuffer + 10, cDigName) == 0) {
			fconf.getline(cBuffer, 64, '\n');
			while (strstr(cBuffer, "END") == NULL) {
				if (strstr(cBuffer, "CHANNEL") != NULL) { // loads processing parameters
					ch = atoi(&cBuffer[8]);
					if ((ch >= MAX_CH) || (ch < 0)) {iFailed |= (1 << config_file_error); throw ProcessorException();}
					fconf.getline(cBuffer, 64, '\n');
					sscanf(cBuffer, "SLOW %i FAST %i PGA %i GAIN_N %f GAIN_Y %f", &iSlowTime[ch], &iFastTime[ch], &iPGASamples[ch], &fGain[ch][0], &fGain[ch][1]);
				}
				fconf.getline(cBuffer, 64, '\n');
			} // end of while
		} // end of if dig
	} // end of file
	fconf.close();
	try {digitizer = shared_ptr<Digitizer>(new Digitizer(cDigName, iSpecial));}
	catch (bad_alloc& ba) {iFailed |= (1 << alloc_error); throw ProcessorException();}
	if (digitizer->Failed()) iFailed |= (1 << method_error);
	if ((digitizer->ID() > 2) && (bMethodActive[XSQ_TF1_t] || bMethodActive[XSQ_NEW_t])) {
		cout << "Chisquare method not compatible with " << digitizer->Name() << '\n';
		bMethodActive[XSQ_TF1_t] = bMethodActive[XSQ_NEW_t] = false;
	}
	return;
}

void Processor::SetDetectorPositions(string in) { // "z0=#,z1=#,z2=#,r0=#,r1=#,r2=#" or some permutation
	if (in == "\0") {
		bPositionsSet = false;
		return;
	}
	if (g_verbose) cout << "Setting detector positions\n";
	unsigned int i(0), iCommas[2] = {0,1};
	string glob; // for dealing with each each substring
	for (i = 0; i < 6; i++) {
		iCommas[1] = in.find(',',iCommas[0]);
		glob = in.substr(iCommas[0],iCommas[1]-iCommas[0]); // grabs input between commas
		if ((glob[0] == 'z') || (glob[0] == 'Z')) fDetectorZ[atoi(glob.c_str()+1)] = atof(glob.c_str()+glob.find('=')+1);
		if ((glob[0] == 'r') || (glob[0] == 'R')) fDetectorR[atoi(glob.c_str()+1)] = atof(glob.c_str()+glob.find('=')+1);
		iCommas[0] = iCommas[1]+1; // shift forward to next substring
		if (iCommas[1] >= in.length()) break;
	}
	bPositionsSet = true;
}

void Processor::SetFileSet(string in) { // also opens raw and processed files
	if (g_verbose) cout << "Opening files\n";
	sRawDataFile = sWorkingDir + "/rawdata/" + in + ".dat";
	sRootFile = sWorkingDir + "/prodata/" + in;
	if ((iSpecial == -1) && (iAverage == 0)) sRootFile += ".root";
	else if ((iSpecial == -1) && (iAverage != 0)) sRootFile += "_a.root";
	else if ((iSpecial != -1) && (iAverage == 0)) sRootFile += "_s.root";
	else sRootFile += "_as.root";
	fin.open(sRawDataFile.c_str(), ios::in | ios::binary);
	if (!fin.is_open()) {
		cout << "Error: " << sRawDataFile << " not found\n";
		iFailed |= (1 << file_error);
		throw ProcessorException();
	}
	f = unique_ptr<TFile>(new TFile(sRootFile.c_str(), "RECREATE"));
	if (!f->IsOpen()) {
		cout << "Error: could not open " << sRootFile << '\n';
		iFailed |= (1 << file_error);
		throw ProcessorException();
	}
}

void Processor::SetSpecials(int special, int average) {
	iSpecial = special;
	iAverage = average;
}
