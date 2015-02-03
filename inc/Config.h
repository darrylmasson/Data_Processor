#ifndef CONFIG_H
#define CONFIG_H

#include <fstream>
#include <cstring>

#ifndef METHOD_H
#include "Method.h"
#endif

// Gets info (config_t and f_header_t structs) from the file.  Leaves the get pointer at the start of the 0th event (Processor() does not seek before reading data)
void GetFileHeader(ifstream* fin, config_t* config, f_header_t* f_header);

// Gets info for processing run
int ParseConfigFile(string& filename, config_t* config, const shared_ptr<Digitizer> digitizer);

#endif // CONFIG_H
