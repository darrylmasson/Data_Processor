#ifndef DIGITIZER_H
#define DIGITIZER_H

#ifndef NGDP_TYPES_H
#include "NGDP_types.h"
#endif

class Digitizer {
	private:
		int 		iFailed;
		char		cName[12];
		double		dSamplerate;
		short		sResolution;
		double		dVpp;
		double		dScaleT;
		double		dScaleV;
		int			iBaselength;
		int			iSpecial;
		dig_id_t	id;
		
	public:
		Digitizer(char in[], int special_in);
		~Digitizer();
		int& Failed()			{return iFailed;}
		char* Name()			{return cName;}
		double& ScaleT()		{return dScaleT;}
		double& ScaleV()		{return dScaleV;}
		double& Samplerate()	{return dSamplerate;}
		short& Resolution()		{return sResolution;}
		double& Vpp()			{return dVpp;}
		int& Baselength()		{return iBaselength;}
		int& Special()			{return iSpecial;}
		dig_id_t& ID()			{return id;}
};

#endif // DIGITIZER_H
