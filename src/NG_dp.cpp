// Neutron generator (et al) data processor
// Written by Darryl, based on algorithms initially written by Jacques and some inspiration from Cassie

#include <cstdlib>
#include <iostream>
#include <unistd.h>
#include <ctime>
#include <chrono>
#include <ratio>

#ifndef PROCESSOR_H
#include "Processor.h"
#endif

using namespace std::chrono;

char error_message[err_dummy_last][64] = {
	"No error",
	"Error in method constructor",
	"Invalid digitizer",
	"Root error",
	"File not found",
	"Allocation error",
	"Thread error",
	"Bad value in file"
};

bool g_verbose(false);

int main(int argc, char **argv) {
	cout << "Neutron generator data processor v3_5\n";
	int i(0), iSpecial(-1), iAverage(0), iErrCode(0);
	string sConfigFile = "\0", sFileset = "\0", sSource = "\0";
	steady_clock::time_point t_start, t_end;
	duration<double> t_elapsed;
	unique_ptr<Processor> processor = nullptr;
	if (argc < 3) {
		cout << "Arguments: -f file [-s source -c config -x special -a moving_average -v]\n";
		return 0;
	}
	while ((i = getopt(argc, argv, "a:c:f:s:vx:")) != -1) { // command line options
		switch(i) {
			case 'a': iAverage = atoi(optarg);	break;
			case 'c': sConfigFile = optarg;		break;
			case 'f': sFileset = optarg;		break;
			case 's': sSource = optarg;			break;
			case 'v': g_verbose = true;			break;
			case 'x': iSpecial = atoi(optarg);	break;
			default: cout << "Arguments: -f file [-s source -c config -x special -a moving_average -v]\n";
				return -1;
	}	}
	if (sFileset == "\0") {
		cout << "No file specified\n";
		return 0;
	}
	if (iAverage < 0) {
		cout << "Invalid moving average: negative\n";
		return 0;
	}
	switch(iSpecial) {
		case -1: break; //cout << "No special options\n"; break; // default
		case 0: cout << "Special samplerate\n"; break; // 500 MSa/s samplerate
		case 1: cout << "Special resolution: 13-bit\n"; break; // 13-bit simulation
		case 2: cout << "Special resolution: 12-bit\n"; break; // 12-bit simulation
		case 3: cout << "Special resolution: 11-bit\n"; break; // 11-bit simulation
		case 4: cout << "Special resolution: 10-bit\n"; break; // 10-bit simulation
		default : cout << "Error: invalid special option specified\n"; return 0;
	}
	try {processor = unique_ptr<Processor>(new Processor(iSpecial, iAverage));}
	catch (bad_alloc& ba) {
		cout << "Error allocating memory\n";
		return 0;
	}
	try { // general setup and preparatory steps
		processor->SetConfigFile(sConfigFile);
		processor->SetFileSet(sFileset);
		processor->SetSource(sSource);
		processor->ParseFileHeader();
		processor->ParseConfigFile();
		processor->ConfigTrees();
		processor->ClassAlloc();
	} catch (ProcessorException& e) {
		cout << e.what();
		cout << error_message[processor->Failed()] << '\n';
		processor.reset();
		return 0;
	}
	
	t_start = steady_clock::now();
	processor->BusinessTime();
	t_end = steady_clock::now();
	if (processor->Failed()) cout << error_message[processor->Failed()] << '\n';
	t_elapsed = duration_cast<duration<double>>(t_end-t_start);
	processor.reset();
	cout << "Total time elapsed: " << t_elapsed.count() << "sec\n";
	return 0;
}
