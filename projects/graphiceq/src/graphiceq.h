
#pragma once

//--------------------------------------------------
// Defines the module's audio processing behavior
//--------------------------------------------------

#include "cephean-lv2.h"
#include "cephean-dsp.h"

using namespace cephean;

//======================================================

//config definitions
#define PLUG_URI "https://github.com/creesullivan/cephean/graphiceq"
#define PLUG_CLASS graphiceq

#define PLUG_CONTROLS (9)
#define PLUG_INPUTS (1)
#define PLUG_OUTPUTS (1)

//======================================================

//class definition
class graphiceq : public plugin
{
public:
	graphiceq(float setFs);
	~graphiceq();

	void init() override;
	void step(int len) override;
	void deactivate() override;

private:
	//ports -----------------
	fcontrol band1{ this, 0 };
	fcontrol band2{ this, 1 };
	fcontrol band3{ this, 2 };
	fcontrol band4{ this, 3 };
	fcontrol band5{ this, 4 };
	fcontrol band6{ this, 5 };
	fcontrol band7{ this, 6 };
	fcontrol band8{ this, 7 };
	fcontrol band9{ this, 8 };

	input in{ this, 9 };

	output out{ this, 10 };

	//constants -------------
	const float slewdur = 0.1f;			//seconds
	const int maxlen = 128;				//maximum block size

	const float fx[8] = { 80.0f, 160.0f, 320.0f, 640.0f, 1280.0f, 2560.0f, 5120.0f, 10240.0f };	//crossover frequencies
	const float f0[9] = { 56.6f, 113.0f, 226.0f, 452.0f, 904.0f, 1808.0f, 3616.0f, 7232.0f, 14464.0f }; //center frequencies

	//state etc.-------------
	rebuffer rebuff;
	fmat sbscr; //subband scratch block <-- LEFT OFF HERE, fmat is untested, but let's probably write the cascade first

	slewedvec<float> gain; //gains that apply the graphic EQ
	
	//helper functions ------
	void update(int len = 0);
};
