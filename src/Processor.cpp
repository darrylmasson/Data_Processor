#include "Processor.h"
#include <iostream>
#include <ctime>
#include <pthread.h>
#include <ratio>
#include <chrono>

using namespace std::chrono;

auto (*root_init[])(TTree*) -> void = {
	CCM::root_init,
	DFT::root_init,
	XSQ::root_init,
	LAP::root_init
};

auto (*root_fill[])() -> void = {
	CCM::root_fill,
	DFT::root_fill,
	XSQ::root_fill,
	LAP::root_fill
};

auto (*root_deinit[])() -> TTree* = {
	CCM::root_deinit,
	DFT::root_deinit,
	XSQ::root_deinit,
	LAP::root_deinit
};

void* Process(void* arg) {
	thread_data_t* input = (thread_data_t*)arg;
	input->event->Set(input->data);
	for (int i = 0; i < NUM_METHODS; i++) if (input->activated[i]) input->methods[i]->evaluate(input->event);
	return nullptr;
}

int Processor(config_t* config, ifstream* fin, TFile* file, Digitizer* dig) {
	int ch(0), ev(0), m(0), rc(0), i_prog_check(0), i_rate(0), i_timeleft(0), ret(no_error), i_livetime(0);
	thread_data_t td[MAX_CH];
	pthread_t threads[MAX_CH];
	unique_ptr<char[]> buffer;
	char treename[NUM_METHODS][4];
	unique_ptr<TTree> TStree = nullptr;
	unique_ptr<TTree> T_data = nullptr; // only one as classes handle trees
	shared_ptr<Digitizer> digitizer(dig);
	unique_ptr<TFile> f(file);
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	void* status = nullptr;
	try {buffer = unique_ptr<char[]>(new char[config->eventsize]);}
	catch (bad_alloc& ba) {
		ret |= alloc_error;
		return ret;
	}
	unsigned long* ul_timestamp = (unsigned long*)(buffer.get() + sizeof(long));
	unsigned long ul_ts_first(0), ul_ts_last(0);
	unsigned short* us_trace = (unsigned short*)(buffer.get() + sizeof_ev_header);
	steady_clock::time_point t_this, t_that;
	duration<double> t_elapsed;
	i_prog_check = config->numEvents/100 + 1;
	f->cd();
	
	for (m = 0; m < NUM_METHODS; m++) {
		sprintf(treename[m], "T%i", m);
		if (config->method_active[m]) {
			try {T_data = unique_ptr<TTree>(new TTree(treename[m], method_names[m]));}
			catch (bad_alloc& ba) {ret |= alloc_error; config->method_active[m] = false;}
			if (T_data->IsZombie()) {ret |= root_error; config->method_active[m] = false;}
			root_init[m](T_data.release());
	}	}
		
	for (ch = 0; ch < config->nchan; ch++) {
		if (g_verbose) cout << "CH" << ch << ": Event ";
		try {
			if (config->average == 0) td[ch].event = shared_ptr<Event>(new Event(config->eventlength, digitizer, config->dc_offset[config->chan[ch]], config->threshold[config->chan[ch]]));
			else td[ch].event = shared_ptr<Event_ave>(new Event_ave(config->eventlength, digitizer, config->dc_offset[config->chan[ch]], config->threshold[config->chan[ch]], config->average));
		} catch (bad_alloc& ba) {
			ret |= alloc_error;
			return ret;
		} if (td[ch].event->Failed()) {
			ret |= method_error;
			return ret;
		}
		if (config->method_active[CCM_t]) {
			if (g_verbose) cout << "CCM ";
			try {td[ch].methods[CCM_t] = shared_ptr<Method>(new CCM(config->chan[ch],config->fastTime[config->chan[ch]],config->slowTime[config->chan[ch]],config->pga_samples[config->chan[ch]],digitizer));}
			catch (bad_alloc& ba) {
				ret |= alloc_error;
				config->method_active[CCM_t] = false;
			} if (td[ch].methods[CCM_t]->Failed()) {
				ret |= method_error;
				config->method_active[CCM_t] = false;
			}
		}
		if (config->method_active[DFT_t]) {
			if (g_verbose) cout << "DFT ";
			try {td[ch].methods[DFT_t] = shared_ptr<Method>(new DFT(config->chan[ch], Event::Length(), digitizer));}
			catch (bad_alloc& ba) {
				ret |= alloc_error;
				config->method_active[DFT_t] = false;
			} if (td[ch].methods[DFT_t]->Failed()) {
				ret |= method_error;
				config->method_active[DFT_t] = false;
			}
		}
		if (config->method_active[XSQ_t]) {
			if (g_verbose) cout << "XSQ ";
			try {td[ch].methods[XSQ_t] = shared_ptr<Method>(new XSQ(config->chan[ch], Event::Length(), config->gain[ch], digitizer));}
			catch (bad_alloc& ba) {
				ret |= alloc_error;
				config->method_active[XSQ_t] = false;
			} if (td[ch].methods[XSQ_t]->Failed()) {
				ret |= method_error;
				config->method_active[XSQ_t] = false;
			}
		}
		if (config->method_active[LAP_t]) {
			if (g_verbose) cout << "LAP ";
			try {td[ch].methods[LAP_t] = shared_ptr<Method>(new LAP(config->chan[ch], Event::Length(), digitizer));}
			catch (bad_alloc& ba) {
				ret |= alloc_error;
				config->method_active[LAP_t] = false;
			} if (td[ch].methods[LAP_t]->Failed()) {
				ret |= method_error;
				config->method_active[LAP_t] = false;
			}
		}
		td[ch].activated = config->method_active;
		td[ch].data = us_trace + ch*config->eventlength;
		if (g_verbose) cout << '\n';
	}
	
	if (!config->already_done) {
		try {TStree = unique_ptr<TTree>(new TTree("TS","Timestamps"));}
		catch (bad_alloc& ba) {
			ret |= alloc_error;
			return ret;
		} if (TStree->IsZombie()) {
			ret |= root_error;
			return ret;
		}
		TStree->Branch("Timestamp", &ul_timestamp[0], "time_stamp/l");
	}
	t_that = steady_clock::now();
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
		for (m = 0; m < NUM_METHODS; m++) if (config->method_active[m]) root_fill[m]();
		if (!config->already_done) TStree->Fill();
		if (ev % i_prog_check == i_prog_check-1) {
			cout << ev*100l/config->numEvents << "%\t\t";
			t_this = steady_clock::now();
			t_elapsed = duration_cast<duration<double>>(t_this-t_that);
			t_that = steady_clock::now();
			i_rate = t_elapsed.count() == 0 ? 9001 : i_prog_check/t_elapsed.count(); // it's OVER 9000!
			cout << i_rate << "\t\t";
			i_timeleft = (config->numEvents - ev)/i_rate;
			cout << i_timeleft << "\n";
		}
		if (ev == 0) ul_ts_first = ul_timestamp[0];	
	}

	cout << "Processing completed\n";
	ul_ts_last = ul_timestamp[0];
	i_livetime = (ul_ts_last - ul_ts_first)/125e6;
	cout << "Acquisition livetime: " << i_livetime << "s\nBeginning cleanup: ";
	fin->close();
	if (!config->already_done) {
		f->cd();
		TStree->Write("",TObject::kOverwrite);
		TStree.reset();
	}
	if (g_verbose) cout << "making friends: ";
	for (m = 0; m < NUM_METHODS; m++) {
		if (config->method_active[m]) { // processed this run
			if (g_verbose) cout << treename[m] << "a ";
			T_data = unique_ptr<TTree>(root_deinit[m]());
			T_data->AddFriend("TS");
			for (int i = 1; i < NUM_METHODS; i++) if ((config->method_done[(m+i)%NUM_METHODS]) || (config->method_active[(m+i)%NUM_METHODS])) T_data->AddFriend(treename[(m+i)%NUM_METHODS]);
			f->cd();
			T_data->Write("",TObject::kOverwrite);
			T_data.reset();
		} else if ((config->method_done[m]) && !(config->method_active[m])) { // processed sometime previous
			if (g_verbose) cout << treename[m] << "b ";
			T_data = unique_ptr<TTree>((TTree*)f->Get(treename[m]));
			for (int i = 1; i < NUM_METHODS; i++) if ((config->method_done[(m+i)%NUM_METHODS]) || (config->method_active[(m+i)%NUM_METHODS])) T_data->AddFriend(treename[(m+i)%NUM_METHODS]);
			f->cd();
			T_data->Write("",TObject::kOverwrite);
			T_data.reset();
		}
	}
	f->Close();
	buffer.reset();
	if (g_verbose) cout << " d'toring classes: ";
	for (ch = 0; ch < config->nchan; ch++) { // general d'tors
		if (g_verbose) cout << "CH" << ch << " ";
		for (m = 0; m < NUM_METHODS; m++) td[ch].methods[m] = nullptr;
		td[ch].event.reset();
	}
	digitizer.reset();
	f.reset();
	cout << " done\n";
	return ret;
}
