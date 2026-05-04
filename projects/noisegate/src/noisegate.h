
#pragma once

//--------------------------------------------------
// Defines the module's audio processing behavior
//--------------------------------------------------

#include "cephean-lv2.h"
#include "cephean-dsp.h"

using namespace cephean;

//======================================================

//config definitions
#define PLUG_URI "https://github.com/creesullivan/cephean/noisegate"
#define PLUG_CLASS noisegate

#define PLUG_CONTROLS (2)
#define PLUG_INPUTS (1)
#define PLUG_OUTPUTS (1)

//======================================================

//class definition
class noisegate : public plugin
{
public:
	noisegate(float setFs);
	~noisegate();

	void init() override;
	void step(int len) override;
	void deactivate() override;

private:
	//ports -----------------
	fcontrol low{ this, 0 };
	fcontrol high{ this, 1 };

	input in{ this, 2 };

	output out{ this, 3 };

	//constants -------------
	const float slew = 0.1f;			//parameter slew rate

	const int maxlen = 128;				//maximum block size

	const float freqlow = 200.0f;		//Hz crossover between low and mid
	const float freqhigh = 2000.0f;		//Hz crossover between mid and high

	const float latency = 0.0005f;		//latency of the system for better onset performance
	const float rmsdur = 0.001f;		//RMS duration of the system
	const float holddur = 0.04f;		//peak hold duration
	const float reldur = 0.04f;			//release duration

	const int Nlat;						//latency in samples

	//state etc.-------------
	rebuffer rebuff;
	fvec inscr, scr1, scr2, scr3;

	//*note* Wiener gains are used here with spectral power subtraction: max{X - G}/X

	sofcasc<2> lowfilt;					//fixed filter cascade for low freq band - 1x LPF and 1x APF
	flatbuffer lowdelay;				//buffer for low band delay to align with power averaging
	lma lowave;							//low band moving power average
	hold4 lowhold;						//low band moving max
	attrel lowrel;						//low band envelope release
	float lowthreshmult = 1.0f;			//low band filter power multiplier for thresh computation
	slewed<float> lowthresh;			//low band smoothed power threshold

	sofcasc<2> midfilt;					//fixed filter cascade for mid freq band - 1x LPF and 1x HPF
	flatbuffer middelay;				//buffer for mid band delay to align with power averaging
	lma midave;							//mid band moving power average
	hold4 midhold;						//mid band moving max
	attrel midrel;						//mid band envelope release
	float midthreshmult = 1.0f;			//mid band filter power multiplier for thresh computation
	slewed<float> midthresh;			//mid band smoothed power threshold

	sofcasc<2> highfilt;				//fixed filter cascade for high freq band - 2x HPF
	flatbuffer highdelay;				//buffer for high band delay to align with power averaging
	lma highave;						//high band moving power average
	hold4 highhold;						//high band moving max
	attrel highrel;						//high band envelope release
	float highthreshmult = 1.0f;		//high band filter power multiplier for thresh computation
	slewed<float> highthresh;			//high band smoothed power threshold

	//helper functions ------
	void update(int len = 0);
	void clear();
};
