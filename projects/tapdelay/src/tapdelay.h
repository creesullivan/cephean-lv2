
#pragma once

//--------------------------------------------------
// Defines the module's audio processing behavior
//--------------------------------------------------

#include "cephean-lv2.h"
#include "cephean-dsp.h"

using namespace cephean;

//======================================================

//config definitions
#define PLUG_URI "https://github.com/creesullivan/cephean/tapdelay"
#define PLUG_CLASS tapdelay

#define PLUG_CONTROLS (7)
#define PLUG_INPUTS (1)
#define PLUG_OUTPUTS (1)

//======================================================

//class definition
class tapdelay : public plugin
{
public:
	tapdelay(float setFs);
	~tapdelay();

	void init() override;
	void step(int len) override;
	void deactivate() override;

private:
	//ports -----------------
	icontrol mode{ this, 0 }; //0 - 1, 1 - 3/4, 2 - 2/3, 3 - 1/2, 4 - 1/3, 5 - 1/4
	fcontrol tempo{ this, 1 };
	fcontrol feedback{ this, 2 };
	fcontrol color{ this, 3 };
	fcontrol rolloff{ this, 4 };
	fcontrol dry{ this, 5 };
	fcontrol wet{ this, 6 };

	input in{ this, 7 };

	output out{ this, 8 };

	//constants -------------
	const float slewdur = 0.1f;			//seconds
	const float tslewdur = 0.2f;		//seconds crossfade for time change

	const float maxdur = 1.0f;			//seconds max time duration in seconds

	const int maxlen = 128;				//maximum block size

	//state etc.-------------
	rebuffer rebuff;
	fvec inscr;
	fvec fbscr;

	slewed<float> slewRolloff;

	pow2buffer buff; //feedback buffer

	gain fbGain;
	faded<int> delay; //delay encoding time difference
	onepolelpf lpf; //color echo lowpass filter
	sof hpf; //rolloff highpass filter

	gain dryGain; //mixing gains
	gain wetGain;

	//helper functions ------
	void update(int len = 0);
	void clear();
};
