#ifndef CONFIG_H
#define CONFIG_H

#include <fstream>
#include <cstring>

#ifndef NGDP_TYPES_H
#include "NGDP_types.h"
#endif

// Gets info (config_t and f_header_t structs) from the file.  Leaves the get pointer at the start of the 0th event (Processor() does not seek before reading data)
void GetFileHeader(ifstream* fin, config_t* config, f_header_t* f_header);

// Gets info for processing run
int ParseConfigFile(string& filename, config_t* config, char dig_name[12]);

#endif // CONFIG_H
