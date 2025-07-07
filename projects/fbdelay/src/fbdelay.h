
#pragma once

//--------------------------------------------------
// Defines the module's audio processing behavior
//--------------------------------------------------

#include "cephean-lv2.h"
#include "cephean-dsp.h"

using namespace cephean;

//======================================================

//config definitions
#define PLUG_URI "https://github.com/creesullivan/cephean/fbdelay"
#define PLUG_CLASS fbdelay

#define PLUG_CONTROLS (6)
#define PLUG_INPUTS (1)
#define PLUG_OUTPUTS (1)

#define PLUG_TESTING_ECHO (false)

//======================================================

//class definition
class fbdelay : public plugin
{
public:
	fbdelay(float setFs, int setSeed=-1);
	~fbdelay();

	void init() override;
	void step(int len) override;
	void deactivate() override;

private:
	//ports -----------------
	fcontrol time{ this, 0 };
	fcontrol feedback{ this, 1 };
	fcontrol wet{ this, 2 };
	fcontrol spread{ this, 3 };
	fcontrol low{ this, 4 };
	fcontrol high{ this, 5 };

	input in{ this, 6 };

	output out{ this, 7 };

	//constants -------------
	const float slewdur = 0.1f;			//seconds for most parameters
	const float tslewdur = 0.2f;		//seconds crossfade for time change

	const float maxdur = 1.0f;			//seconds max time duration

	const int minecho = 480;			//minimum sample echo delay @ 48k
	const int maxecho = 4800;			//maximum sample echo delay @ 48k
	const float dryecho = 0.1f;			//echo dry mix gain for stability
	const float recho = 1.6f;			//empirical echo power compensation parameter #1
	const float gecho = 2.1f;			//empirical echo power compensation parameter #2

	const int maxlen = 128;				//maximum block size

	//state etc.-------------
	rebuffer rebuff;

	circbuffer buff; //feedback buffer

	gain fbGain;
	faded<int> delay; //delay encoding time difference
	sof hpf; //highpass filter
	diffuser echo; //echo filter
	sofcasc<2> lpf; //feedback lowpass filter

	gain wetGain;

	//helper functions ------
	void update(int len = 0);
};
