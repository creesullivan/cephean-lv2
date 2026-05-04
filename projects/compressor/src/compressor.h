
#pragma once

//--------------------------------------------------
// Defines the module's audio processing behavior
//--------------------------------------------------

#include "cephean-lv2.h"
#include "cephean-dsp.h"

using namespace cephean;

//======================================================

//config definitions
#define PLUG_URI "https://github.com/creesullivan/cephean/compressor"
#define PLUG_CLASS compressor

#define PLUG_CONTROLS (4)
#define PLUG_INPUTS (1)
#define PLUG_OUTPUTS (1)

//======================================================

//class definition
class compressor : public plugin
{
public:
	compressor(float setFs);
	~compressor();

	void init() override;
	void step(int len) override;
	void deactivate() override;

private:
	//ports -----------------
	fcontrol drive{ this, 0 };
	fcontrol tone{ this, 1 };
	fcontrol dry{ this, 2 };
	fcontrol wet{ this, 3 };

	input in{ this, 4 };

	output out{ this, 5 };

	//constants -------------
	const float slewdur = 0.1f;			//seconds

	const int maxlen = 128;				//maximum block size

	const float Trms = 0.05f;			//RMS time constant in seconds
	const float Tatt = 0.0001f;			//limiter attack time in seconds
	const float Trel = 0.005f;			//limiter release time in seconds
	const float crest = 0.15f;			//short-term crest factor target as inverse power multiplier

	const float knee = 0.5f;			//amplitude knee width of compression curve

	const float tonefreq = 1000.0f;		//Hz fulcrum point of tone filter

	//state etc.-------------
	rebuffer rebuff;
	fvec inscr, scr1, scr2, scr3;

	slewed<float> slewTone;				//slewed tone adjustment parameter
	fof prefilt, postfilt;				//pre- and post-filter for color processing
	lma rmsdetect;						//core RMS detector
	attrel release;						//post asym release for limit smoothing etc

	softclip curve;						//compression curve with soft knee

	gain driveGain, dryGain, wetGain;	//gains for drive, dry, and wet parameters

	//helper functions ------
	void update(int len = 0);
	void clear();
};
