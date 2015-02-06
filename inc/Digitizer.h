#ifndef DIGITIZER_H
#define DIGITIZER_H

extern bool g_verbose;

enum dig_id_t {
  dt5751 = 0,
  dt5751des,
  dt5730,
  v1724,
  invalid_dig
};

class Digitizer {
	private:
		int 		failed;
		char		name[12];
		double		samplerate;
		short		resolution;
		double		V_pp;
		double		scaleT;
		double		scaleV;
		int			baselength;
		int			special;
		dig_id_t	id; // an easier way of identifying digitizers
		
	public:
		Digitizer(char in[], int special_in);
		~Digitizer();
		int& Failed()			{return failed;}
		char* Name()			{return name;}
		double& ScaleT()		{return scaleT;}
		double& ScaleV()		{return scaleV;}
		double& Samplerate()	{return samplerate;}
		short& Resolution()		{return resolution;}
		double& Vpp()			{return V_pp;}
		int& Baselength()		{return baselength;}
		int& Special()			{return special;}
		dig_id_t& ID()			{return id;}
		
		const int& Failed() const			{return failed;}
		const char* Name() const			{return name;}
		const double& ScaleT() const		{return scaleT;}
		const double& ScaleV() const		{return scaleV;}
		const double& Samplerate() const	{return samplerate;}
		const short& Resolution() const		{return resolution;}
		const double& Vpp() const			{return V_pp;}
		const int& Baselength() const		{return baselength;}
		const int& Special() const			{return special;}
		const dig_id_t& ID() const			{return id;}
};

#endif // DIGITIZER_H
