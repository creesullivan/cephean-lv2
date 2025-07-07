
#pragma once

//--------------------------------------------------
// Defines the module's audio processing behavior
//--------------------------------------------------

#include "cephean-lv2.h"
#include "cephean-dsp.h"

using namespace cephean;

//======================================================

//config definitions
#define PLUG_URI "https://github.com/creesullivan/cephean/tremvibe"
#define PLUG_CLASS tremvibe

#define PLUG_CONTROLS (4)
#define PLUG_INPUTS (1)
#define PLUG_OUTPUTS (1)

//======================================================

//class definition
class tremvibe : public plugin
{
public:
	tremvibe(float setFs);
	~tremvibe();

	void init() override;
	void step(int len) override;
	void deactivate() override;

private:
	//ports -----------------
	fcontrol tremelo{ this, 0 };
	fcontrol vibrato{ this, 1 };
	fcontrol rate{ this, 2 };
	fcontrol random{ this, 3 };

	input in{ this, 4 };

	output out{ this, 5 };

	//constants -------------
	const float slewdur = 0.1f;			//seconds

	const float maxdel = 0.02f; //seconds, under-supporting wide vibrato at low rates
	const int maxdelN; //samples

	const int maxlen = 256;				//maximum block size

	//state etc.-------------
	rebuffer rebuff;

	float tdepth = 0.0f; //linear tremelo depth multiplier (+/-)
	float vdepth = 0.0f; //linear vibrato peak step multiplier
	float dph = 0.0f; //rate as a phase differential
	float rmult = 0.0f; //how much randomness is dialed in each period (+/-)

	float gT = 0.0f;
	float gV = 0.0f;
	float curdph = 0.0f;

	fastcircfloat ph{ 2.0f * constants.pi, 0.0f }; //phase for oscillator

	flatbuffer vbuff; //vibrato delay buffer

	//helper functions ------
	void update(int len = 0);
	
	void newperiod();
};
