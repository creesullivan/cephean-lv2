
#pragma once

//--------------------------------------------------
// Defines the module's audio processing behavior
//--------------------------------------------------

#include "cephean-lv2.h"
#include "cephean-dsp.h"

using namespace cephean;

//======================================================

//config definitions
#define PLUG_URI "https://github.com/creesullivan/cephean/linkwitz"
#define PLUG_CLASS linkwitz

#define PLUG_CONTROLS (1)
#define PLUG_INPUTS (2)
#define PLUG_OUTPUTS (1)

//======================================================

//class definition
class linkwitz : public plugin
{
public:
	linkwitz(float setFs);
	~linkwitz();

	void init() override;
	void step(int len) override;
	void deactivate() override;

private:
	//ports -----------------
	fcontrol frequency{ this, 0 };

	input inLow{ this, 1 };
	input inHigh{ this, 2 };

	output out{ this, 3 };

	//constants -------------
	const float slewdur = 0.1f;			//seconds

	const int maxlen = 128;				//maximum block size

	//state etc.-------------
	rebuffer rebuff;
	fvec inscr1, inscr2;

	slewed<float> slewFreq;
	sof lpf, hpf;

	//helper functions ------
	void update(int len = 0);
	void clear();
};
