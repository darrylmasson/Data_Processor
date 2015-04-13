#include "Digitizer.h"
#include <cstring>
#include <iostream>

Digitizer::Digitizer(char in[], int special_in) {
	strcpy(cName,in);
	iSpecial = special_in;
	if (strcmp(cName, "DT5730") == 0) {
		dSamplerate = 5E8;
		sResolution = (1 << 14);
		if (iSpecial > 0) sResolution = (1 << (14 - iSpecial));
		dVpp = 2;
		iBaselength = 20;
		iFailed = 0;
		id = dt5730;
	} else if (strcmp(cName, "DT5751") == 0) {
		dSamplerate = (iSpecial == 0) ? 5E8 : 1E9;
		sResolution = (1 << 10);
		dVpp = 1;
		iBaselength = (iSpecial == 0) ? 20 : 40;
		iFailed = 0;
		id = dt5751;
	} else if (strcmp(cName, "DT5751DES") == 0) {
		dSamplerate = 2E9;
		sResolution = (1 << 10);
		dVpp = 1;
		iBaselength = 80;
		iFailed = 0;
		id = dt5751des;
	} else if (strcmp(cName, "V1724") == 0) {
		dSamplerate = 1E8;
		sResolution = (1 << 14); 
		dVpp = 2.25;
		iBaselength = 15;
		iFailed = 0;
		id = v1724;
	} else {
		dSamplerate = -1;
		sResolution = -1;
		dVpp = -1;
		iBaselength = 1;
		iFailed |= (1 << dig_error);
		id = invalid_dig;
	}
	
	dScaleV = dVpp/(double)sResolution; // volts/bin
	dScaleT = 1E9/dSamplerate; // ns
}

Digitizer::~Digitizer() {
	if (g_verbose) std::cout << " dig d'tor ";
}
