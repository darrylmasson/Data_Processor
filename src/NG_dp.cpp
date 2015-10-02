// Neutron generator (et al) data processor
// Written by Darryl, based on algorithms initially written by Jacques and some inspiration from Cassie

#include <ctime>
#include <chrono>
#include <ratio>

#ifndef PROCESSOR_H
#include "Processor.h"
#endif

using namespace std::chrono;

bool g_verbose(false);

void PrintVersions() {
	cout << "Versions installed:\n";
	cout << "Event: " << Event::sfVersion << '\n';
	cout << "Method: " << Method::sfVersion << '\n';
	cout << "Discriminator: " << Discriminator::sfVersion << '\n';
}

int main(int argc, char **argv) {
	cout << "Neutron generator raw data processor v4\n";
	int i(0), iElapsed(0), iAverage(0), iLevel(0);
	string sConfigFile = "NG_dp_config.cfg", sFileset = "\0", sSource = "\0", sDetectorPos = "\0";
	const string sArgs = "Arguments: -f file [-s source -l level -c config -a moving_average -p detector_positions -v -e]";
	steady_clock::time_point t_start, t_end;
	duration<double> t_elapsed;
	Processor processor;
	if (argc < 2) {
		cout << sArgs << '\n';
		return 0;
	}
	while ((i = getopt(argc, argv, "a:c:ef:s:p:vx:")) != -1) { // command line options
		switch(i) {
			case 'a': iAverage = atoi(optarg);	break;
			case 'c': sConfigFile = optarg;		break;
			case 'e': g_verbose = true;			break;
			case 'f': sFileset = optarg;		break;
			case 'l': iLevel = atoi(optarg);	break;
			case 'p': sDetectorPos = optarg;	break;
			case 's': sSource = optarg;			break;
			case 'v': PrintVersions();			return 0;
//			case 'x': iSpecial = atoi(optarg);	break;
			default: cout << sArgs << '\n';
				return -1;
	}	}
	if (sFileset == "\0") {
		cout << "No file specified\n";
		return 0;
	}
/*	switch(iSpecial) {
		case -1: break; //cout << "No special options\n"; break; // default
		case 0: cout << "Special samplerate\n"; break; // 500 MSa/s samplerate
		case 1: cout << "Special resolution: 13-bit\n"; break; // 13-bit simulation
		case 2: cout << "Special resolution: 12-bit\n"; break; // 12-bit simulation
		case 3: cout << "Special resolution: 11-bit\n"; break; // 11-bit simulation
		case 4: cout << "Special resolution: 10-bit\n"; break; // 10-bit simulation
		default : cout << "Error: invalid special option specified\n"; return 0;
	} */

	try { // general setup and preparatory steps
		processor.SetParams(iAverage, iLevel);
		processor.SetConfigFile(sConfigFile);
		processor.SetSource(sSource);
		processor.SetDetectorPositions(sDetectorPos);
		processor.Setup(sFileset);
	} catch (ProcessorException& e) {
		cout << e.what();
		cout << "Setup failed, exiting\n";
		return 0;
	}

	t_start = steady_clock::now();
	processor.BusinessTime();
	t_end = steady_clock::now();
	t_elapsed = duration_cast<duration<double>>(t_end-t_start);
	iElapsed = t_elapsed.count();
	cout << "Total time elapsed: " << iElapsed/3600 << 'h' << (iElapsed%3600)/60 << 'm' << iElapsed%60 << "s\n";
	return 0;
}