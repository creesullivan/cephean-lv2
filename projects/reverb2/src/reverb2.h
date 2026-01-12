
#pragma once

//--------------------------------------------------
// Defines the module's audio processing behavior
//--------------------------------------------------

#include "cephean-lv2.h"
#include "cephean-dsp.h"

using namespace cephean;

//======================================================

//config definitions
#define PLUG_URI "https://github.com/creesullivan/cephean/reverb2"
#define PLUG_CLASS reverb2

#define PLUG_CONTROLS (6)
#define PLUG_INPUTS (1)
#define PLUG_OUTPUTS (1)

//======================================================

//class definition
class reverb2 : public plugin
{
public:
	reverb2(double setFs);
	~reverb2();

	void init() override;
	void step(int len) override;
	void deactivate() override;

private:
	//ports -----------------
	fcontrol tail{ this, 0 };
	fcontrol color{ this, 1 };
	fcontrol rolloff{ this, 2 };
	fcontrol energy{ this, 3 };
	fcontrol wet{ this, 4 };
	fcontrol dry{ this, 5 };

	input in{ this, 6 };

	output out{ this, 7 };

	//constants -------------
	const float slewdur = 0.1f;			//seconds
	const int maxlen = 128;				//samples

	//state etc.-------------
	rebuffer rebuff;
	
	slewed<float> scolor;
	slewed<float> srolloff;

	waveverb2<32> verb;
	sofcasc<2> bpf; //HPF for rolloff, LPF for color
	gain gainWet;
	gain gainDry;

	fvec inscr;

	//helper functions ------
	void update(int len = 0);
	void clear();
};
