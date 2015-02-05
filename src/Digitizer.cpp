#include "Digitizer.h"
#include <cstring>
#include <iostream>

Digitizer::Digitizer(char in[], int special_in) {
	strcpy(name,in);
	special = special_in;
	if (strcmp(name, "DT5730") == 0) {
		samplerate = 5E8;
		resolution = (1 << 14);
		if (special > 0) resolution = (1 << (14 - special));
		V_pp = 2;
		baselength = 20;
		failed = 0;
		id = dt5730;
	} else if (strcmp(name, "DT5751") == 0) {
		samplerate = (special == 0) ? 5E8 : 1E9;
		resolution = (1 << 10);
		V_pp = 1;
		baselength = (special == 0) ? 20 : 40;
		failed = 0;
		id = dt5751;
	} else if (strcmp(name, "DT5751DES") == 0) {
		samplerate = 2E9;
		resolution = (1 << 10);
		V_pp = 1;
		baselength = 80;
		failed = 0;
		id = dt5751des;
	} else if (strcmp(name, "V1724") == 0) {
		samplerate = 1E8;
		resolution = (1 << 14); 
		V_pp = 2.25;
		baselength = 15;
		failed = 0;
		id = v1724;
	} else {
		samplerate = -1;
		resolution = -1;
		V_pp = -1;
		baselength = 1;
		failed = 2; // 2 = digitizer error
		id = invalid_dig;
	}
	
	scaleV = V_pp/(double)resolution;
	scaleT = 1E9/samplerate; // ns
}

Digitizer::~Digitizer() {
	std::cout << " dig d'tor ";
}
