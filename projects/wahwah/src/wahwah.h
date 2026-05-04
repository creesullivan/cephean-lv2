
#pragma once

//--------------------------------------------------
// Defines the module's audio processing behavior
//--------------------------------------------------

#include "cephean-lv2.h"
#include "cephean-dsp.h"

using namespace cephean;

//======================================================

//config definitions
#define PLUG_URI "https://github.com/creesullivan/cephean/wahwah"
#define PLUG_CLASS wahwah

#define PLUG_CONTROLS (7)
#define PLUG_INPUTS (1)
#define PLUG_OUTPUTS (1)

//======================================================

//class definition
class wahwah : public plugin
{
public:
	wahwah(float setFs);
	~wahwah();

	void init() override;
	void step(int len) override;
	void deactivate() override;

private:
	//ports -----------------
	fcontrol ccontrol{ this, 0 };
	fcontrol cfreq1{ this, 1 };
	fcontrol cfreq2{ this, 2 };
	fcontrol cres1{ this, 3 };
	fcontrol cres2{ this, 4 };
	fcontrol crolloff{ this, 5 };
	fcontrol clevel{ this, 6 };

	input in{ this, 7 };

	output out{ this, 8 };

	//constants -------------
	const float slewdur = 0.1f;			//seconds parameter slew

	const float bypatt = 0.1f;			//attack time out of an auto bypass
	const float byprel = 0.5f;			//release time into an auto bypass
	const float bypdur = 1.0f;			//seconds we wait to bypass automatically when the control sits below controlthresh
	const float bypthresh = 1.0f;		//control percent value that starts counting down to bypass

	const int maxlen = 128;				//maximum block size
	const int Nbypass = 96000;			//samples corresponding to bypassdur

	//state etc.-------------
	rebuffer rebuff;
	fvec inscr, scr1;

	//references to the current value of filter shape definition controls for recomputation speed
	float vfreq1 = 0.0f;
	float vfreq2 = 0.0f;
	float vres1 = 0.0f;
	float vres2 = 0.0f;

	slewed<float> slewFreq; //slewed filter center frequency in log units
	slewed<float> slewRes; //slewed filter resonance in log units
	slewed<float> slewRolloff; //slewed rolloff control

	int nbyp = 0; //counter for auto-bypass

	sof resonator; //core resonator filter for the classic wah
	sof lpf; //LR LPF for rolloff management
	sof hpf; //LR HPF for rolloff management
	gain outGain;

	float byptarg = 1.0f; //bypass gain target (1 for applied, 0 for bypassed)
	attrel bypgain;	//asymmetric bypass gain

	//helper functions ------
	void update(int len = 0);
	void clear();
};
