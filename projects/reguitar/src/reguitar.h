
#pragma once

//--------------------------------------------------
// Defines the module's audio processing behavior
//--------------------------------------------------

#include "cephean-lv2.h"
#include "cephean-dsp.h"

using namespace cephean;

//======================================================

//config definitions
#define PLUG_URI "https://github.com/creesullivan/cephean/reguitar"
#define PLUG_CLASS reguitar

#define PLUG_CONTROLS (5)
#define PLUG_INPUTS (1)
#define PLUG_OUTPUTS (1)

//======================================================

//class definition
class reguitar : public plugin
{
public:
	reguitar(float setFs);
	~reguitar();

	void init() override;
	void step(int len) override;
	void deactivate() override;

private:
	//ports -----------------
	fcontrol freq1{ this, 0 };
	fcontrol res1{ this, 1 };
	fcontrol freq2{ this, 2 };
	fcontrol res2{ this, 3 };
	fcontrol gain{ this, 4 };

	input in{ this, 5 };

	output out{ this, 6 };

	//constants -------------
	const float slewdur = 0.1f;			//seconds

	const int maxlen = 128;				//maximum block size

	//state etc.-------------
	rebuffer rebuff;
	fvec inscr;

	sof filt; //core filter that applies the entire effect lol

	//helper functions ------
	void update(int len = 0);
	void clear();
};
