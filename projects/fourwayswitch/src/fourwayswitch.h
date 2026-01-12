
#pragma once

//--------------------------------------------------
// Defines the module's audio processing behavior
//--------------------------------------------------

#include "cephean-lv2.h"
#include "cephean-dsp.h"

using namespace cephean;

//======================================================

//config definitions
#define PLUG_URI "https://github.com/creesullivan/cephean/fourwayswitch"
#define PLUG_CLASS fourwayswitch

#define PLUG_CONTROLS (1)
#define PLUG_INPUTS (4)
#define PLUG_OUTPUTS (1)

//======================================================

//class definition
class fourwayswitch : public plugin
{
public:
	fourwayswitch(float setFs);
	~fourwayswitch();

	void init() override;
	void step(int len) override;
	void deactivate() override;

private:
	//ports -----------------
	fcontrol position{ this, 0 };

	input in1{ this, 1 };
	input in2{ this, 2 };
	input in3{ this, 3 };
	input in4{ this, 4 };

	output out{ this, 5 };

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
