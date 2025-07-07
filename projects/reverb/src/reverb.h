
#pragma once

//--------------------------------------------------
// Defines the module's audio processing behavior
//--------------------------------------------------

#include "cephean-lv2.h"
#include "cephean-dsp.h"

using namespace cephean;

//======================================================

//config definitions
#define PLUG_URI "https://github.com/creesullivan/cephean/reverb"
#define PLUG_CLASS reverb

#define PLUG_CONTROLS (5)
#define PLUG_INPUTS (1)
#define PLUG_OUTPUTS (1)

//======================================================

//class definition
class reverb : public plugin
{
public:
	reverb(double setFs, int setSeed = -1);
	~reverb();

	void init() override;
	void step(int len) override;
	void deactivate() override;

private:
	//ports -----------------
	icontrol seed{ this, 0 };
	fcontrol duration{ this, 1 };
	fcontrol color{ this, 2 };
	fcontrol wet{ this, 3 };
	fcontrol dry{ this, 4 };

	input in{ this, 5 };

	output out{ this, 6 };

	//constants -------------
	const float slewdur = 0.1f;			//seconds
	const int maxlen = 128;				//samples

	//state etc.-------------
	rebuffer rebuff;
	
	slewed<float> scolor;

	waveverb verb;
	sofcasc<2> bpf; //LPF/HPF
	gain gainWet;
	gain gainDry;

	fvec dryscr; //dry path scratch

	//helper functions ------
	void update(int len = 0);
};
