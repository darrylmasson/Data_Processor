#ifndef PROCESSOR_H
#define PROCESSOR_H

#include <fstream>
#ifndef CCM_H
#include "CCM.h"
#endif
#ifndef DFT_H
#include "DFT.h"
#endif
#ifndef XSQ_H
#include "XSQ.h"
#endif

// Argument is thread_data_t* passed as void*
void* Process(void* arg);

// Main processor function.  Returns error code
int Processor(config_t* config, ifstream* fin, TFile* file, Digitizer* dig, bool verbose);

#endif // PROCESSOR_H
