
#pragma once

//--------------------------------------------------
// Defines the module's audio processing behavior
//--------------------------------------------------

#include "cephean-lv2.h"
#include "cephean-dsp.h"

using namespace cephean;

//======================================================

//config definitions
#define PLUG_URI "https://github.com/creesullivan/cephean/flanger"
#define PLUG_CLASS flanger

#define PLUG_CONTROLS (6)
#define PLUG_INPUTS (1)
#define PLUG_OUTPUTS (1)

//======================================================

//class definition
class flanger : public plugin
{
public:
	flanger(float setFs);
	~flanger();

	void init() override;
	void step(int len) override;
	void deactivate() override;

private:
	//ports -----------------
	fcontrol depth{ this, 0 };
	fcontrol stop1{ this, 1 };
	fcontrol stop2{ this, 2 };
	fcontrol rate{ this, 3 };
	fcontrol tone{ this, 4 };
	fcontrol level{ this, 5 };

	input in{ this, 6 };

	output out{ this, 7 };

	//constants -------------
	const float slewdur = 0.1f;			//seconds
	const int maxlen = 128;				//maximum block size

	const float maxFreq = 8000.0f;		//max frequency in Hz
	const float minFreq = 400.0f;		//min frequency in Hz

	//state etc.-------------
	rebuffer rebuff;

	flangercore fcore; //core delay line of the flanger

	slewed<float> dph; //current oscillator radian frequency (ie phase step)
	slewed<float> del1, del2; //slewed stop delays
	gain dryGain, wetGain; //tone & level management

	fastcircfloat ph; //oscillator phase on 0 - 2pi

	fvec inscr;

	//helper functions ------
	void update(int len = 0);
	void clear();
};
