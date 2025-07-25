
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
	fcontrol* bands[9]; //for convenience

	input in{ this, 9 };

	output out{ this, 10 };

	//constants -------------
	const float slewdur = 0.1f;						//seconds
	const int maxlen = 128;							//maximum block size
	const int Nband = 9;							//number of bands

	const float Qp = 1.0f;				//peaking filter Q
	const float Qs = 0.707f;				//shelf filter Q

	const float fx[8] = { 50.0f, 106.6f, 227.2f, 484.3f, 1032.4f, 2200.7f, 4691.2f, 10000.0f };	//crossover frequencies
	const float f0[9] = { 34.2f, 73.0f, 155.6f, 331.7f, 707.1f, 1507.3f, 3213.1f, 6849.2f, 14600.2f }; //center frequencies

	const float lev0 = powf(10.0f, -12.0f / 20.0f); //maps to 0% band gain
	const float lev25 = powf(10.0f, -6.0f / 20.0f); //maps to 25% band gain
	const float lev75 = powf(10.0f, 6.0f / 20.0f); //maps to 75% band gain
	const float lev100 = powf(10.0f, 12.0f / 20.0f); //maps to 100% band gain

	//state etc -------------
	rebuffer rebuff;

	sofcasc<9> casc; //9-band cascade
	
	//helper functions ------
	void update(int len = 0);
};
