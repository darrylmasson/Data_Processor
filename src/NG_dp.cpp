// Neutron generator (et al) data processor
// Written by Darryl, based on algorithms initially written by Jacques and some inspiration from Cassie

#include <ctime>
#include <chrono>
#include <ratio>
#include <getopt.h> // getopt_long

#ifndef PROCESSOR_H
#include "Processor.h"
#endif
#ifndef DISCRIMINATOR_H
#include "Discriminator.h"
#endif

using namespace std::chrono;

int g_verbose(0);

void PrintVersions() {
	cout << "Versions installed:\n";
	cout << "Event: " << Event::sfVersion << '\n';
	cout << "Method: " << Method::sfVersion << '\n';
	cout << "Discriminator: " << Discriminator::sfVersion << '\n';
}

void Help() {
	cout << "Arguments:\n";
	cout << "--file\t\tSpecifies which file for processor to use. Required\n";
	cout << "--source\t\tSpecifies what source was used for file. Required\n";
	cout << "--moving_average\tSpecifies how many samples on each side of a given point to average. Optional\n";
	cout << "--config\t\tSpecifies a non-default configuration file to use. Must be in same directory as default. Optional\n";
	cout << "--verbose\t\tSets level of output during processing to 1 (default 0). Optional\n";
	cout << "--very_verbose\t\tSets level of output during processing to 2 (default 0). Optional\n";
	cout << "-h, --help\t\tPrints this message\n";
	cout << "--NG_current\tCurrent setpoint on neutron generator. Requried for NG runs, optional otherwise\n";
	cout << "--level\t\tSpecifies processing level to be done. Default 1 (Method and Event), 0 = Event, 2 = Event, Method, and Discriminator, 3 = Discriminator. Optional\n";
	cout << "--old\t\tForces old-format file header\n";
	cout << "--position\t\tSets positions of detectors. Required for Cf-252 runs, optional otherwise\n";
	cout << "--version\t\tPrints installed versions and exits\n";
	cout << "--NG_voltage\tVoltage setpoint on neutron generator. Requried for NG runs, optional otherwise\n";
	return;
}

int main(int argc, char **argv) {
	cout << "Neutron generator raw data processor v4\n";
	int i(0), iElapsed(0), iAverage(0), iLevel(1), option_index(0);
	string sConfigFile = "NG_dp_config.cfg", sFileset = "\0", sSource = "\0", sDetectorPos = "\0";
	steady_clock::time_point t_start, t_end;
	duration<double> t_elapsed;
	float fHV(0), fCurrent(0);
	Processor processor;
	Discriminator discriminator;
	option long_options[] = {
		{"very_verbose", no_argument, &g_verbose, 2},
		{"verbose", no_argument, &g_verbose, 1},
		{"help", no_argument, 0, 'h'},
		{"file", required_argument, 0, 'f'},
		{"source", required_argument, 0, 's'},
		{"moving_average", required_argument, 0, 'a'},
		{"config", required_argument, 0, 'c'},
		{"NG_current", required_argument, 0, 'I'},
		{"level", required_argument, 0, 'l'},
		{"position", required_argument, 0, 'p'},
		{"version", no_argument, 0, 'v'},
		{"NG_voltage", required_argument, 0, 'V'},
		{"old", no_argument, 0, 'o'},
		{0,0,0,0}
	};
	if (argc < 2) {
		Help();
		return 1;
	}
	while ((i = getopt_long(argc, argv, "a:c:e:f:hI:l:os:p:vV:", long_options, &option_index)) != -1) { // command line options
		switch(i) {
			case 0: break;
			case 'a': iAverage = atoi(optarg);	break;
			case 'c': sConfigFile = optarg;		break;
			case 'e': g_verbose = atoi(optarg);	break;
			case 'f': sFileset = optarg;		break;
			case 'h': Help();					return 1;
			case 'I': fCurrent = atof(optarg);	break;
			case 'l': iLevel = atoi(optarg);	break;
			case 'o': processor.ForceOld();		break;
			case 'p': sDetectorPos = optarg;	break;
			case 's': sSource = optarg;			break;
			case 'v': PrintVersions();			return 1;
			case 'V': fHV = atof(optarg);		break;
			default: Help();					return 1;
	}	}
	if (sFileset == "\0") {
		cout << "No file specified\n";
		return 1;
	} else if (sSource == "\0") {
		cout << "No source specified\n";
		return 1;
	}

	if (iLevel <= 2) {
		try { // general setup and preparatory steps
			processor.SetParams(iAverage, iLevel);
			processor.SetConfigFile(sConfigFile);
			processor.SetSource(sSource);
			processor.SetDetectorPositions(sDetectorPos);
			processor.SetNGSetpoint(fHV, fCurrent);
			processor.Setup(sFileset);
		} catch (ProcessorException& e) {
			cout << e.what();
			cout << "Setup failed, exiting\n";
			return 1;
	}	}
	if (iLevel >= 2) {
		discriminator.Setup(sFileset);
		if (discriminator.Failed()) {
			cout << "Discriminator setup failed\n";
			if (iLevel > 2) {
				cout << "Returning\n";
				return 1;
			} else {
				cout << "Continuing\n";
				iLevel = 1;
	}	}	}

	t_start = steady_clock::now();
	if (iLevel <= 2) processor.BusinessTime();
	if (iLevel >= 2) discriminator.Discriminate();
	t_end = steady_clock::now();
	t_elapsed = duration_cast<duration<double>>(t_end-t_start);
	iElapsed = t_elapsed.count();
	cout << "Total time elapsed: " << iElapsed/3600 << 'h' << (iElapsed%3600)/60 << 'm' << iElapsed%60 << "s\n";
	return 0;
}
