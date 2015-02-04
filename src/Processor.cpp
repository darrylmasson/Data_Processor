#include "Processor.h"
#include <iostream>
#include <ctime>
#include <pthread.h>

void* Process(void* arg) {
	thread_data_t* input = (thread_data_t*)arg;
	input->event->Set(input->data);
	for (int i = 0; i < NUM_METHODS; i++) if (input->activated[i]) input->methods[i]->evaluate(input->event);
	return nullptr;
}

int Processor(config_t* config, ifstream* fin, shared_ptr<TFile> f, shared_ptr<Digitizer> digitizer, bool verbose) {
	int ch(0), ev(0), m(0), rc(0), prog_check(0), rate(0), timeleft(0), ret(no_error), livetime(0);
	thread_data_t td[MAX_CH];
	pthread_t threads[MAX_CH];
	unique_ptr<char[]> buffer;
	char treename[NUM_METHODS][4];
	unique_ptr<TTree> TStree;
	shared_ptr<TTree> T_data[NUM_METHODS];
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	void* status = nullptr;
	try {buffer = unique_ptr<char[]>(new char[config->eventsize]);}
	catch (bad_alloc& ba) {
		ret |= alloc_error;
		return ret;
	}
	unsigned long* timestamp = (unsigned long*)(buffer.get() + sizeof(long));
	unsigned long ts_first(0), ts_last(0);
	unsigned short* trace = (unsigned short*)(buffer.get() + sizeof_ev_header);
	clock_t time_this = clock(), time_last = clock();
	prog_check = config->numEvents/100 + 1;
	
	f->cd();
	
	for (m = 0; m < NUM_METHODS; m++) {
		sprintf(treename[m], "T%i", m);
		if (config->method_active[m]) {
			try {T_data[m] = shared_ptr<TTree>(new TTree(treename[m], method_names[m]));}
			catch (bad_alloc& ba) {ret |= alloc_error; config->method_active[m] = false;}
			if (T_data[m]->IsZombie()) {ret |= root_error; config->method_active[m] = false;}
	}	}
	if (config->method_active[CCM_t]) CCM::root_init(T_data[CCM_t]);
	if (config->method_active[DFT_t]) DFT::root_init(T_data[DFT_t]);
	if (config->method_active[XSQ_t]) XSQ::root_init(T_data[XSQ_t]);
	
	for (ch = 0; ch < config->nchan; ch++) {
		if (verbose) cout << "CH" << ch << ": Event ";
		try {td[ch].event = shared_ptr<Event>(new Event(config->eventlength, digitizer, config->dc_offset[config->chan[ch]], config->threshold[config->chan[ch]]));}
		catch (bad_alloc& ba) {
			ret |= alloc_error;
			return ret;
		} if (td[ch].event->Failed()) {
			ret |= method_error;
			return ret;
		}
		if (config->method_active[CCM_t]) {
			if (verbose) cout << "CCM ";
			try {td[ch].methods[CCM_t] = unique_ptr<Method>(new CCM(config->chan[ch],config->fastTime[config->chan[ch]],config->slowTime[config->chan[ch]],config->pga_samples[config->chan[ch]],digitizer));}
			catch (bad_alloc& ba) {
				ret |= alloc_error;
				config->method_active[CCM_t] = false;
			} if (td[ch].methods[CCM_t]->Failed()) {
				ret |= method_error;
				config->method_active[CCM_t] = false;
			}
		}
		if (config->method_active[DFT_t]) {
			if (verbose) cout << "DFT ";
			try {td[ch].methods[DFT_t] = unique_ptr<Method>(new DFT(config->chan[ch], Event::Length()));}
			catch (bad_alloc& ba) {
				ret |= alloc_error;
				config->method_active[DFT_t] = false;
			} if (td[ch].methods[DFT_t]->Failed()) {
				ret |= method_error;
				config->method_active[DFT_t] = false;
			}
		}
		if (config->method_active[XSQ_t]) {
			if (verbose) cout << "XSQ ";
			try {td[ch].methods[XSQ_t] = unique_ptr<Method>(new XSQ(config->chan[ch], Event::Length(), config->gain[ch], digitizer));}
			catch (bad_alloc& ba) {
				ret |= alloc_error;
				config->method_active[XSQ_t] = false;
			} if (td[ch].methods[XSQ_t]->Failed()) {
				ret |= method_error;
				config->method_active[XSQ_t] = false;
			}
		}
		td[ch].activated = config->method_active;
		td[ch].data = trace + ch*config->eventlength;
		if (verbose) cout << '\n';
	}
	
	if (!config->already_done) {
		f->cd();
		try {TStree = unique_ptr<TTree>(new TTree("TS","Timestamps"));}
		catch (bad_alloc& ba) {
			ret |= alloc_error;
			return ret;
		} if (TStree->IsZombie()) {
			ret |= root_error;
			return ret;
		}
		TStree->Branch("Timestamp", &timestamp[0], "time_stamp/l");
	}
	time_last = clock();
	cout << "Processing:\n";
	cout << "Completed\tRate (ev/s)\tTime left (s)\n";
	for (ev = 0; ev < config->numEvents; ev++) {
		fin->read(buffer.get(), config->eventsize);
		if (config->method_active[XSQ_t]) for (ch = 0; ch < config->nchan; ch++) {
			td[ch].event->Set(td[ch].data);
			for (m = 0; m < NUM_METHODS; m++) if (config->method_active[m]) td[ch].methods[m]->evaluate(td[ch].event); // TF1 isn't thread-friendly
		} else {
			for (ch = 0; ch < config->nchan; ch++) if ( (rc = pthread_create(&threads[ch], &attr, Process, (void*)&td[ch])) ) {ret |= thread_error; return ret;}
			for (ch = 0; ch < config->nchan; ch++) if ( (rc = pthread_join(threads[ch], &status)) ) {ret |= thread_error; return ret;}
		}
		for (m = 0; m < NUM_METHODS; m++) if (config->method_active[m]) T_data[m]->Fill();
		if (!config->already_done) TStree->Fill();
		if (ev % prog_check == prog_check-1) {
			cout << ev*100l/config->numEvents << "%\t\t";
			time_this = clock();
			rate = prog_check*CLOCKS_PER_SEC/(time_this - time_last);
			time_last = time_this;
			cout << rate << "\t\t";
			timeleft = (config->numEvents - ev)/rate;
			cout << timeleft << "\n";
		}
		if (ev == 0) ts_first = timestamp[0];	
	}

	cout << "Processing completed\n";
	ts_last = timestamp[0];
	livetime = (ts_last - ts_first)/125e6;
	cout << "Acquisition livetime: " << livetime << "s\nBeginning cleanup: ";
	fin->close();
	if (!config->already_done) {
		TStree->Write("",TObject::kOverwrite);
		TStree.reset();
	}
	if (verbose) cout << "making friends: ";
	for (m = 0; m < NUM_METHODS; m++) {
		if (config->method_active[m]) {
			if (verbose) cout << treename[m] << " ";
			T_data[m]->AddFriend("TS");
			for (int i = 1; i < NUM_METHODS; i++) if ((config->method_done[(m+i)%NUM_METHODS]) || (config->method_active[(m+i)%NUM_METHODS])) T_data[m]->AddFriend(treename[(m+i)%NUM_METHODS]);
			T_data[m]->Write("",TObject::kOverwrite);
			T_data[m].reset();
		} else if ((config->method_done[m]) && !(config->method_active[m])) {
			if (verbose) cout << treename[m] << " ";
			T_data[m] = shared_ptr<TTree>((TTree*)f->Get(treename[m]));
			if (T_data[m].use_count() == 0) continue;
			for (int i = 1; i < NUM_METHODS; i++) if ((config->method_done[(m+i)%NUM_METHODS]) || (config->method_active[(m+i)%NUM_METHODS])) T_data[m]->AddFriend(treename[(m+i)%NUM_METHODS]);
			T_data[m]->Write("",TObject::kOverwrite);
			T_data[m].reset();
		}
	}
	f->Close();
	buffer.reset();
	if (verbose) cout << " d'toring classes ";
	for (ch = 0; ch < config->nchan; ch++) {
//		for (m = 0; m < NUM_METHODS; m++) td[ch].methods[m].reset();
//		td[ch].event.reset();
	}
//	digitizer.reset();
//	f.reset();
	cout << " done\n";
	return ret;
}
