#include "Config.h"
#include <cstdlib>

void GetFileHeader(ifstream* fin, config_t* config, f_header_t* f_header) {
	long filesize(0l);
	int i(0);
	char buffer[sizeof_f_header];
	fin->seekg(0, fin->end);
	filesize = fin->tellg();
	fin->seekg(0, fin->beg);
	fin->read(buffer,sizeof_f_header);
	strncpy(f_header->dig_name, buffer, sizeof(f_header->dig_name)); // digitizer name
	
	memcpy(&f_header->mask, buffer + 12, sizeof(f_header->mask)); // channel mask
	config->mask = f_header->mask;
	
	memcpy(&f_header->ev_len, buffer + 14, sizeof(f_header->ev_len)); // eventlength
	config->eventlength = f_header->ev_len;
		
	memcpy(&f_header->trig_post, buffer + 18, sizeof(f_header->trig_post)); // post-trigger
	config->trig_post = f_header->trig_post;
	
	memcpy(f_header->dc_off, buffer + 22, sizeof(f_header->dc_off)); // dc offsets
	memcpy(config->dc_offset, f_header->dc_off, sizeof(f_header->dc_off));
	
	memcpy(f_header->threshold, buffer + 54, sizeof(f_header->threshold)); // trigger thresholds
	memcpy(config->threshold, f_header->threshold, sizeof(f_header->threshold));
	
	config->nchan = 0;
	for (i = 0; i < MAX_CH; i++) if (f_header->mask & (1<<i)) config->chan[config->nchan++] = i;
	config->eventsize = sizeof_ev_header + config->nchan*f_header->ev_len*2; // 2 bytes per sample
	config->numEvents = (filesize - sizeof_f_header)/config->eventsize;
	
}

int ParseConfigFile(string& filename, config_t* config, char dig_name[12]) {
	char file[64];
	sprintf(file, "%s/config/%s", path, filename.c_str());
	ifstream fin(file,ios::in);
	if (!fin.is_open()) return file_error;

	char buffer[64] = {'\0'}, temp[32] = {'\0'};
	int ch(-1), code(0);
	while (!fin.eof()) {
		fin.getline(buffer, 64, '\n');
		if (buffer[0] == '#') continue;
		if (strcmp(buffer, "METHODS") == 0) { // checks for active processing methods
			fin.getline(buffer, 64, '\n');
			while (strstr(buffer, "END") == NULL) {
				sscanf(buffer, "%s %i", temp, &code);
				for (int i = 0; i < NUM_METHODS; i++) if (strcmp(temp, method_names[i]) == 0) config->method_active[i] = code;
				fin.getline(buffer, 64, '\n');
			} // end of while
		} // end of if method
		if (strcmp(buffer + 10, dig_name) == 0) {
			fin.getline(buffer, 64, '\n');
			while (strstr(buffer, "END") == NULL) {
				if (strstr(buffer, "CHANNEL") != NULL) { // loads processing parameters
					ch = atoi(&buffer[8]);
					if ((ch >= MAX_CH) || (ch < 0)) {fin.close(); return config_file_error;}
					fin.getline(buffer, 64, '\n');
					sscanf(buffer, "SLOW %i FAST %i PGA %i GAIN_N %f GAIN_Y %f", &config->slowTime[ch], &config->fastTime[ch], &config->pga_samples[ch], &config->gain[ch][0], &config->gain[ch][1]);
					}
				fin.getline(buffer, 64, '\n');
			} // end of while
		} // end of if dig
	} // end of file
	fin.close();

	return no_error;
}
