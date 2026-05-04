
#pragma once

//--------------------------------------------------
// Defines the module's audio processing behavior
//--------------------------------------------------

#include "cephean-lv2.h"
#include "cephean-dsp.h"

using namespace cephean;

//======================================================

//config definitions
#define PLUG_URI "https://github.com/creesullivan/cephean/autotune"
#define PLUG_CLASS autotune

#define PLUG_CONTROLS (6)
#define PLUG_INPUTS (1)
#define PLUG_OUTPUTS (1)

//======================================================

//class definition
class autotune : public plugin
{
public:
	autotune(float setFs);
	~autotune();

	void init() override;
	void step(int len) override;
	void deactivate() override;

private:
	//ports -----------------
	icontrol fundamental{ this, 0 };	//0-11, cast to scale::fundamental
	icontrol mode{ this, 1 };			//0-8, cast to scale::mode
	fcontrol tuning{ this, 2 };			//A4 in Hz, 440 typical
	icontrol pitch1{ this, 3 };			//pitch minimum as a note number, C2 - E6
	icontrol pitch2{ this, 4 };			//pitch maximum as a note number, C2 - E6
	fcontrol smooth{ this, 5 };			//smoothing time as %

	input in{ this, 6 };

	output out{ this, 7 };

	//constants -------------
	const int maxlen = 128;				//maximum block size

	const float thresh = 1e-8f;				//linear power thresh for correlation safety (-80 dB)

	const int bufflen = 2048;				//length of buffer including chunk at nominal rate
	const int chunklen = 448;				//length of chunk at nominal rate
	const int blendlen = 512;				//length of crossfading at nominal rate
	const int minper = 36;					//minimum period length at nominal rate
	const int DSR = 8;						//downsampling rate of correlation at nominal rate
	const int Rlen = (2048 - 448) / 8 + 1;	//length of correlation sequence at downsampled rate
	const int maxper = ((2048 - 448) / 8 + 1) * 8; //absolute maximum period length at nominal rate

	const float Tblend = 0.95f;				//blending allowance
	const float Ttrans = 0.5f;				//correlation threshold indicating a transient

	const float maxShift = powf(2.0f, 200.0f / 1200.0f); //theoretical worst case shift size

	//state etc.-------------
	rebuffer rebuff;
	fvec inscr;

	fastcircfloat ph; //test object

	corrbuffer buff; //buffer for computing correlation for pitch detection and shifting
	varshift shift; //core variable pitch shifter

	int minPitchPeriod = 0; //corresponding to pitch2
	int maxPitchPeriod = 0; //corresponding to pitch1
	scale curscale; //defines the current tuning, fundamental, and scale mode

	float curdsamp = 1.0f; //smoothed pitch shifting step size in samples/sample
	float salpha = 0.0f; //manual one pole smoother
	float sbeta = 1.0f;

	fvec R; //scratch storage for autocorrelation sequence
	fvec scr;

	//helper functions ------
	void update(int len = 0);
	void clear();
};
