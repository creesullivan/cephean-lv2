
#pragma once

//--------------------------------------------------
// Defines the module's audio processing behavior
//--------------------------------------------------

#include "cephean-lv2.h"
#include "cephean-dsp.h"

using namespace cephean;

//======================================================

//config definitions
#define PLUG_URI "https://github.com/creesullivan/cephean/diffuse"
#define PLUG_CLASS diffuse

#define PLUG_CONTROLS (5)
#define PLUG_INPUTS (1)
#define PLUG_OUTPUTS (1)

//======================================================

//class definition
class diffuse : public plugin
{
public:
	diffuse(float setFs);
	~diffuse();

	void init() override;
	void step(int len) override;
	void deactivate() override;

private:
	//ports -----------------
	icontrol seed{ this, 0 };
	fcontrol spread{ this, 1 };
	fcontrol color{ this, 2 };
	fcontrol wet{ this, 3 };
	fcontrol dry{ this, 4 };

	input in{ this, 5 };

	output out{ this, 6 };

	//constants -------------
	const float slewdur = 0.1f;			//seconds

	const int maxlen = 256;				//maximum block size

	//state etc.-------------
	rebuffer rebuff;

	fvec dryscr;				//dry path scratch
	slewed<float> slewColor;	//manual block slew for color 

	diffuser alg;				//core diffusion short verb algorithm
	sofcasc<2> lpf;				//color lowpass filter
	sof apf;					//phase alignment allpass for dry path
	gain wetGain;				//mixing gains
	gain dryGain;

	//helper functions ------
	void update(int len = 0);
};
