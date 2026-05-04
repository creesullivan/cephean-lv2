
#pragma once

//--------------------------------------------------
// Defines the module's audio processing behavior
//--------------------------------------------------

#include "cephean-lv2.h"
#include "cephean-dsp.h"

using namespace cephean;

//======================================================

//config definitions
#define PLUG_URI "https://github.com/creesullivan/cephean/ditherdist"
#define PLUG_CLASS ditherdist

#define PLUG_CONTROLS (8)
#define PLUG_INPUTS (1)
#define PLUG_OUTPUTS (1)

//======================================================

//class definition
class ditherdist : public plugin
{
public:
	ditherdist(float setFs);
	~ditherdist();

	void init() override;
	void step(int len) override;
	void deactivate() override;

private:
	//ports -----------------
	fcontrol drive{ this, 0 };
	fcontrol bias{ this, 1 };
	fcontrol gate{ this, 2 };
	fcontrol low{ this, 3 };
	fcontrol high{ this, 4 };
	fcontrol tone{ this, 5 };
	fcontrol dry{ this, 6 };
	fcontrol wet{ this, 7 };

	input in{ this, 8 };

	output out{ this, 9 };

	//constants -------------
	const float slewdur = 0.1f;			//seconds

	const int maxlen = 128;				//maximum block size

	//state etc.-------------
	rebuffer rebuff;
	fvec inscr;

	slewedvec<float> slewfilt; //slewed filter parameters: low, high, tone
	sofcasc<2> prefilt; //cascade of FOF-HPF + FOF-tone, and SOF-LPF

	timedither dither;
	transistor waveshaper;
	
	sofcasc<2> postfilt; //cascade of FOF-HPF + FOF-itone, and SOF-LPF

	gain dryGain, wetGain; //slewed dry and wet gain

	//helper functions ------
	void update(int len = 0);
	void clear();
};
