// Neutron generator (et al) data processor
// Written by Darryl, based on algorithms initially written by Jacques and some inspiration from Cassie

#include <cstdlib>
#include <iostream>
#include <unistd.h>
#include <ctime>

#ifndef PROCESSOR_H
#include "Processor.h"
#endif
#ifndef CONFIG_H
#include "Config.h"
#endif

const float method_versions[NUM_METHODS] = {
	CCM::version,
	DFT::version,
	XSQ::version
};

int main(int argc, char **argv) {
	cout << "Neutron generator data processor v3_5\n";
	int i(0), err_code(0), timenow(0), datenow(0), special(-1), method_id(0), pga_check[MAX_CH], fast_check[MAX_CH], slow_check[MAX_CH], XSQ_ndf(0);
	float version(0);
	string config_file = "\0", fileset = "\0";
	config_t config;
	f_header_t f_header;
	char filename[64], source[12], methodname[12];
	time_t rawtime;
	struct tm* today;
	bool update(0), checked[NUM_METHODS], verbose(0);
	unique_ptr<TFile> f = nullptr;
	unique_ptr<TTree> tx = nullptr, tc = nullptr, tv = nullptr;
	shared_ptr<Digitizer> digitizer;
	ifstream fin;
	memset(source, ' ', sizeof(source));
	memset(checked, 0, sizeof(checked));
	if (argc < 3) {
		cout << "Arguments: -f file [-s source -c config -x special]\n";
		return 0;
	}
	while ((i = getopt(argc, argv, "c:f:s:vx:")) != -1) {
		switch(i) {
			case 'c': config_file = optarg;		break;
			case 'f': fileset = optarg;			break;
			case 's': strcpy(source,optarg);	break;
			case 'v': verbose = true;			break;
			case 'x': special = atoi(optarg);	break;
			default: cout << "Arguments: -f file [-s source -c config -x special]\n";
				return -1;
	}	}
	if (fileset == "\0") {
		cout << "No file specified\n";
		return 0;
	}
	switch(special) {
		case -1: break; //cout << "No special options\n"; break; // default
		case 0: cout << "Special samplerate\n"; break; // 500 MSa/s samplerate
		case 1: cout << "Special resolution: 13-bit\n"; break; // 13-bit simulation
		case 2: cout << "Special resolution: 12-bit\n"; break; // 12-bit simulation
		case 3: cout << "Special resolution: 11-bit\n"; break; // 11-bit simulation
		case 4: cout << "Special resolution: 10-bit\n"; break; // 10-bit simulation
		default : cout << "Error: invalid special option specified\n"; return 0;
	}
	memset(&config, 0, sizeof(config_t));
	
	if (special == -1) sprintf(filename, "%sprodata/%s.root", path, fileset.c_str());
	else sprintf(filename, "%sprodata/%s_x.root", path, fileset.c_str());
	fin.open(filename, ios::in);
	config.already_done =  fin.is_open();
	if (fin.is_open()) fin.close();

	sprintf(filename, "%srawdata/%s.dat", path, fileset.c_str());
	fin.open(filename, ios::in | ios::binary);
	if (!fin.is_open()) {
		for (i = 0; i < MAX_CH; i++) {
			sprintf(filename, "%srawdata/%s_%i.dat", path, fileset.c_str(), i);
			fin.open(filename, ios::in);
			if (fin.is_open()) {
				cout << "Legacy files no longer supported, please run up-converter\n";
				fin.close();
				return 0;
		}	}
		cout << "Error: " << fileset << " not found\n";
		return 0;
	}
	GetFileHeader(&fin, &config, &f_header);
	if (config.numEvents == 0) {
		cout << "Error: file header not processed\n";
		return 0;
	}
	
	try {digitizer = shared_ptr<Digitizer>(new Digitizer(f_header.dig_name, special));}
	catch (bad_alloc& ba) {
		cout << "Could not instantiate Digitizer class\n";
		return 0;
	} if ( (err_code = digitizer->Failed()) ) {
		cout << error_message[err_code] << "\n";
		return 0;
	}
	
	if (config_file == "\0") config_file = "NG_dp_config.cfg";
	cout << "Parsing " << config_file << " with settings for " << digitizer->Name() << ": ";
	if ( (err_code = ParseConfigFile(config_file, &config, digitizer)) ) {
		cout << "failed: " << error_message[err_code] << '\n';
		return 0;
	} else cout << "done\n";
	
	cout << "Found " << config.numEvents << " events and " << config.nchan << " channels\n";
	
	time(&rawtime);
	today = localtime(&rawtime);
	timenow = (today->tm_hour)*100 + today->tm_min; // hhmm
	datenow = (today->tm_year-100)*10000 + (today->tm_mon+1)*100 + today->tm_mday; // yymmdd
	
	switch (digitizer->ID()) {
		case dt5751 : 
			if (digitizer->Special() == 0) XSQ_ndf = min(config.eventlength/2, 225)-3; // length - number of free parameters
			else XSQ_ndf = min(config.eventlength, 450)-3;
			break;
		case dt5751des :
			XSQ_ndf = min(config.eventlength, 899)-3;
			break;
		case dt5730 :
			XSQ_ndf = min(config.eventlength, 225)-3;
			break;
		case v1724 :
		default :
			XSQ_ndf = -1;
			break;
	}
	
	if (special == -1) sprintf(filename, "%sprodata/%s.root", path, fileset.c_str());
	else sprintf(filename, "%sprodata/%s_x.root", path, fileset.c_str());
	try {f = unique_ptr<TFile>(new TFile(filename, "UPDATE"));}
	catch (bad_alloc& ba) {
		cout << "Allocation error\n";
		return 0;
	} if (f->IsZombie()) {
		cout << "Error: could not open root file\n";
		return 0;
	}
	f->cd();
	if (config.already_done) {
		tx.reset((TTree*)f->Get("Tx"));
		if (tx) {
			cout << fileset << " already processed with an older version of NG_dp; v3.5 will require the removal of previous results.  Continue <y|n>? ";
			char overwrite(0);
			cin >> overwrite;
			if (overwrite == 'y') {
				f->Delete("*;*");
				config.already_done = false;
			} else {
				f->Close();
				f.reset();
				return 0;
			}
		}
	tx.reset();
	}
	
	if (config.already_done) {
		cout << "Processing already done, checking versions...\n";
		tv = unique_ptr<TTree>((TTree*)f->Get("TV"));
		tc = unique_ptr<TTree>((TTree*)f->Get("TC"));
		tv->SetBranchAddress("MethodID", &method_id);
		tv->SetBranchAddress("Version", &version);
		tc->SetBranchAddress("PGA_samples", pga_check);
		tc->SetBranchAddress("Fast_window", fast_check);
		tc->SetBranchAddress("Slow_window", slow_check);
		tc->GetEntry(tc->GetEntries()-1);
		for (i = tv->GetEntries() - 1; i >= 0; i--) { // most recent entry will be last
			tv->GetEntry(i); // so it's better to start from the end rather than the beginning
			update = false;
			config.method_done[method_id] = (checked[method_id] ? config.method_done[method_id] : true);
			if (config.method_active[method_id] && !checked[method_id]) {
				checked[method_id] = true;
				config.method_done[method_id] = true;
				update |= (version < method_versions[method_id]);
				update |= ((method_id == CCM_t) && (memcmp(pga_check, config.pga_samples, sizeof(pga_check)) != 0));
				update |= ((method_id == CCM_t) && (memcmp(fast_check, config.fastTime, sizeof(fast_check)) != 0));
				update |= ((method_id == CCM_t) && (memcmp(slow_check, config.slowTime, sizeof(slow_check)) != 0));
				if (update) {
					cout << "Method " <<  method_names[method_id] << " will be updated\n";
					sprintf(filename, "T%i;*", method_id);
					f->Delete(filename);
					config.method_done[method_id] = false;
				} else {
					cout << "Method " << method_names[method_id] << " does not need updating, reprocess anyway <y|n>? ";
					char overwrite(0);
					cin >> overwrite;
					if (overwrite == 'y') {
						sprintf(filename, "T%i;*", method_id);
						f->Delete(filename);
						config.method_done[method_id] = false;
					} else config.method_active[method_id] = false;
			}	}
		}
		if (memchr(config.method_active, 1, NUM_METHODS) == NULL) cout << "No processing methods activated\n";
		tv->SetBranchAddress("Date", &datenow);
		tv->SetBranchAddress("Time", &timenow);
		tv->SetBranchAddress("MethodName", methodname);
		tc->SetBranchAddress("PGA_samples", config.pga_samples);
		tc->SetBranchAddress("Fast_window", config.fastTime);
		tc->SetBranchAddress("Slow_window", config.slowTime);
		for (i = 0; i < NUM_METHODS; i++) if (config.method_active[i]) {
			method_id = i;
			version = method_versions[i];
			strcpy(methodname, method_names[i]);
			tv->Fill();
			if (i == CCM_t) tc->Fill();
		}
		tc->Write("",TObject::kOverwrite);
		tv->Write("",TObject::kOverwrite);
	} else {
		cout << "Creating config trees\n";
		try {
			tx = unique_ptr<TTree>(new TTree("TI", "Info"));
			tc = unique_ptr<TTree>(new TTree("TC","CCM,PGA info"));
			tv = unique_ptr<TTree>(new TTree("TV", "Versions"));
		} catch (bad_alloc& ba) {
			cout << "Could not create config trees\n";
			f->Close();
			return 0;
		} if ((tx->IsZombie()) || (tc->IsZombie()) || (tv->IsZombie())) {
			cout << "Could not create config trees\n";
			f->Close();
			return 0;
		}
		tx->Branch("Digitizer", f_header.dig_name, "name[12]/B");
		tx->Branch("Source", source, "source[12]/B");
		tx->Branch("ChannelMask", &config.mask, "mask/s");
		tx->Branch("TriggerThreshold", f_header.threshold, "threshold[8]/i");
		tx->Branch("DC_offset", f_header.dc_off, "dc_off[8]/i"); // the numbers in this tree
		tx->Branch("Posttrigger", &config.trig_post, "tri_post/I"); // don't change, so it's only
		tx->Branch("Eventlength", &config.eventlength, "ev_len/I"); // written out once
		tx->Branch("Chisquared_NDF", &XSQ_ndf, "ndf/I");
		if (special != -1) tx->Branch("Special", &special, "special/I");
		
		tv->Branch("MethodName", methodname, "codename[12]/B");
		tv->Branch("MethodID", &method_id, "codeid/I");
		tv->Branch("Date", &datenow, "date/I"); // this tree holds version info
		tv->Branch("Time", &timenow, "time/I");
		tv->Branch("Version", &version, "version/F");
		
		tc->Branch("PGA_samples", config.pga_samples, "pga[8]/I");
		tc->Branch("Fast_window", config.fastTime, "fast[8]/I");
		tc->Branch("Slow_window", config.slowTime, "slow[8]/I");
		for (i = 0; i < NUM_METHODS; i++) if (config.method_active[i]) {
			strcpy(methodname, method_names[i]);
			version = method_versions[i];
			method_id = i;
			tv->Fill();
		}
		tx->Fill();
		tc->Fill();
		tx->Write();
		tv->Write();
		tc->Write();
	}
	tx.reset();
	tc.reset();
	tv.reset();
	
	if (memchr(config.method_active, 1, NUM_METHODS) == NULL) {
		f->Close();
		f.reset();
		digitizer.reset();
		return 0;
	} else {
		cout << "Processing methods activated: ";
		for (i = 0; i < NUM_METHODS; i++) if (config.method_active[i]) cout << method_names[i] << " ";
		cout << '\n';
	}
	clock_t t = clock();
	if ( (err_code = Processor(&config, &fin, f.release(), digitizer, verbose)) ) cout << error_message[err_code] << '\n';

//	digitizer.reset();
	t = clock() - t;
	cout << "Total time elapsed: " << t/CLOCKS_PER_SEC << "sec\n";
	return 0;
}
