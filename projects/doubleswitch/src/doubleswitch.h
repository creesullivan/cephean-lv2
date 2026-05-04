
#pragma once

//--------------------------------------------------
// Defines the module's audio processing behavior
//--------------------------------------------------

#include "cephean-lv2.h"
#include "cephean-dsp.h"

using namespace cephean;

//======================================================

//config definitions
#define PLUG_URI "https://github.com/creesullivan/cephean/doubleswitch"
#define PLUG_CLASS doubleswitch

#define PLUG_CONTROLS (2)
#define PLUG_INPUTS (4)
#define PLUG_OUTPUTS (1)

//======================================================

//class definition
class doubleswitch : public plugin
{
public:
	doubleswitch(float setFs);
	~doubleswitch();

	void init() override;
	void step(int len) override;
	void deactivate() override;

private:
	//ports -----------------
	fcontrol switch1{ this, 0 };
	fcontrol switch2{ this, 1 };

	input in1{ this, 2 };
	input in2{ this, 3 };
	input in3{ this, 4 };
	input in4{ this, 5 };

	output out{ this, 6 };

	//constants -------------
	const float slewdur = 0.03f;		//seconds

	const int maxlen = 128;				//maximum block size

	//state etc.-------------
	rebuffer rebuff;
	fvec inscr1, inscr2, inscr3, inscr4;

	gain gain1, gain2, gain3, gain4;

	//helper functions ------
	void update(int len = 0);
	void clear();
};
