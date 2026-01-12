
#pragma once

//--------------------------------------------------
// Defines the module's audio processing behavior
//--------------------------------------------------

#include "cephean-lv2.h"
#include "cephean-dsp.h"

using namespace cephean;

//======================================================

//config definitions
#define PLUG_URI "https://github.com/creesullivan/cephean/distortion2"
#define PLUG_CLASS distortion2

#define PLUG_CONTROLS (6)
#define PLUG_INPUTS (1)
#define PLUG_OUTPUTS (1)

//======================================================

//class definition
class distortion2 : public plugin
{
public:
	distortion2(float setFs);
	~distortion2();

	void init() override;
	void step(int len) override;
	void deactivate() override;

private:
	//ports -----------------
	fcontrol drive{ this, 0 };
	fcontrol shape{ this, 1 };
	fcontrol tone { this, 2 };
	fcontrol dirty{ this, 3 };
	fcontrol color{ this, 4 };
	fcontrol level{ this, 5 };

	input in{ this, 6 };

	output out{ this, 7 };

	//constants -------------
	const float slewdur = 0.1f;				//seconds duration for parameter slews
	const int maxlen = 128;					//maximum block size

	const float attdur = 0.002f;			//attack envelope time in seconds
	const float holddur = 0.02f;			//hold envelope time in seconds
	const float reldur = 0.02f;				//release envelope time in seconds

	const float ftone = 1000.0f;			//fulcrum frequency for tone control

	//state etc.-------------
	rebuffer rebuff;
	fvec inscr;

	slewed<float> slewDrive;
	slewed<float> slewShape;
	slewed<float> slewDirty;
	slewed<float> slewColor;
	slewed<float> slewTone;

	gain driveGain; //drive gain into the nonlinearity
	fof preTone; //preemphasis tone filter

	attrel envAttack; //envelope detection for input compressor
	hold4 envHold;
	attrel envRelease;
	
	fastratio compress; //compression curve
	simshaper distort; //distortion curve with lowpass filtering
	gain makeupGain; //makeup gain after the nonlinearity
	sof makeupFilter; //rolloff frequency makeup boost
	fof postTone; //deemphasis tone filter

	gain outGain; //final output level gain

	fvec scr; //general purpose scratch vectors

	//helper functions ------
	void update(int len = 0);
	void clear();
};
