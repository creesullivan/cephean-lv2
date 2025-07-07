
#pragma once

//--------------------------------------------------
// Defines the module's audio processing behavior
//--------------------------------------------------

#include "cephean-lv2.h"
#include "cephean-dsp.h"

using namespace cephean;

//======================================================

//config definitions
#define PLUG_URI "https://github.com/creesullivan/cephean/distortion"
#define PLUG_CLASS distortion

#define PLUG_CONTROLS (6)
#define PLUG_INPUTS (1)
#define PLUG_OUTPUTS (1)

//======================================================

//class definition
class distortion : public plugin
{
public:
	distortion(float setFs);
	~distortion();

	void init() override;
	void step(int len) override;
	void deactivate() override;

private:
	//ports -----------------
	fcontrol drive{ this, 0 };
	fcontrol shape{ this, 1 };
	fcontrol dirty{ this, 2 };
	fcontrol detail{ this, 3 };
	fcontrol color{ this, 4 };
	fcontrol level{ this, 5 };

	input in{ this, 6 };

	output out{ this, 7 };

	//constants -------------
	const float slewdur = 0.1f;			//seconds

	const int maxlen = 128;				//maximum block size
	const int rate = 2;					//canonical rate to apply nonlinearity at

	const float fpre = 1000.0f;			//preemphasis frequency in Hz
	const float gpre = 0.3f;			//preemphasis depth multiplier

	const float attdur = 0.001f;		//attack envelope time in seconds
	const float holddur = 0.01f;		//hold envelope time in seconds
	const float reldur = 0.02f;			//release envelope time in seconds

	//state etc.-------------
	rebuffer rebuff;

	slewed<float> slewdetail;
	slewed<float> slewcolor;

	int usr; //upsampling rate
	float Fsu; //upsampled sampling rate
	int maxlenu; //upsampled maximum block size

	fof preemph;
	sofcasc<2> filtDetail; //detail LPF

	upsampler up;
	sofcasc<4> aau; //upsampling antialiasing filter

	attrel envAttack; //envelope detection for compressor
	hold4 envHold;
	attrel envRelease;

	fastratio curve; //compression/distortion gain curve
	gain mixComp; //mixing gains for compression vs distortion
	gain mixDist;
	sofcasc<2> filtColor; //color LPF
	sof phaseColor; //phase matching APF
	
	sofcasc<4> aad; //downsampling antialiasing filter
	downsampler down;

	fof deemph;
	gain outLevel;

	fvec uscr1, uscr2, uscr3; //upsampled scratch vectors

	//helper functions ------
	void update(int len = 0);
};
