
#pragma once

//--------------------------------------------------
// Defines the module's audio processing behavior
//--------------------------------------------------

#include "cephean-lv2.h"
#include "cephean-dsp.h"

using namespace cephean;

//======================================================

//config definitions
#define PLUG_URI "https://github.com/creesullivan/cephean/selfmod"
#define PLUG_CLASS selfmod

#define PLUG_CONTROLS (5)
#define PLUG_INPUTS (1)
#define PLUG_OUTPUTS (1)

//======================================================

//class definition
class selfmod : public plugin
{
public:
	selfmod(float setFs);
	~selfmod();

	void init() override;
	void step(int len) override;
	void deactivate() override;

private:
	//ports -----------------
	fcontrol depth{ this, 0 };
	fcontrol bandwidth{ this, 1 };
	fcontrol rolloff{ this, 2 };
	fcontrol wet{ this, 3 };
	fcontrol dry{ this, 4 };

	input in{ this, 5 };

	output out{ this, 6 };

	//constants -------------
	const float slewdur = 0.1f;			//seconds

	const float ctrlfreq = 80.0f;		//Hz ctrl frequency peak, -6 dB above and below
	const float attdur = 0.001f;		//seconds
	const float reldur = 0.01f;			//seconds
	const float maxboost = 10000.0f;	//+80 dB available control signal gain

	const int maxlen = 128;				//maximum block size

	//state etc.-------------
	rebuffer rebuff;
	fvec inscr, scr1, scr2, scr3, scr4;

	sofcasc<2> filt0; //allpass cascade for 0deg path and dry signal
	sofcasc<4> filt90; //cascade for 90deg path
	sof filtdir0; //0deg post filter to form modulatee signal
	sof filtdir90; //90deg post filter to form modulatee signal
	fof filtapf; //0deg post filter to form dry signal
	sof filtctrl; //0deg post filter to form control signal

	attrel detect; //primary attack/release envelope detector
	saturator limit; //safety limiting clipper on the control envelope

	slewed<float> slewDepth; //depth multiplier in rads vs a unity amplitude control signal

	slewed<float> slewRolloff;
	sof lpf; //rolloff lowpass
	sof hpf0; //rolloff highpass on 0 deg path
	sof hpf90; //rolloff highpass on 90 deg path

	gain wetGain;
	gain dryGain;

	//helper functions ------
	void update(int len = 0);
	void clear();
};
