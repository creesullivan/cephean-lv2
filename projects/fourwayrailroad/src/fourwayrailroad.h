
#pragma once

//--------------------------------------------------
// Defines the module's audio processing behavior
//--------------------------------------------------

#include "cephean-lv2.h"
#include "cephean-dsp.h"

using namespace cephean;

//======================================================

//config definitions
#define PLUG_URI "https://github.com/creesullivan/cephean/fourwayrailroad"
#define PLUG_CLASS fourwayrailroad

#define PLUG_CONTROLS (1)
#define PLUG_INPUTS (1)
#define PLUG_OUTPUTS (4)

//======================================================

//class definition
class fourwayrailroad : public plugin
{
public:
	fourwayrailroad(float setFs);
	~fourwayrailroad();

	void init() override;
	void step(int len) override;
	void deactivate() override;

private:
	//ports -----------------
	fcontrol position{ this, 0 };

	input in{ this, 1 };

	output out1{ this, 2 };
	output out2{ this, 3 };
	output out3{ this, 4 };
	output out4{ this, 5 };

	//constants -------------
	const float slewdur = 0.03f;		//seconds

	const int maxlen = 128;				//maximum block size

	//state etc.-------------
	rebuffer rebuff;

	gain gain1, gain2, gain3, gain4;

	//helper functions ------
	void update(int len = 0);
	void clear();
};
