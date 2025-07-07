
#pragma once

//--------------------------------------------------
// Defines the module's audio processing behavior
//--------------------------------------------------

#include "cephean-lv2.h"
#include "cephean-dsp.h"
#include "cephean-fft.h"

//DEVELOPMENT PLAN:
//	-linear phase timbre shifting <-- IN PROGRESS
//	-min phase timbre shifting
//	-transient manip
//	-victory lap with noise gating and level

using namespace cephean;

//======================================================

//config definitions
#define PLUG_URI "https://github.com/creesullivan/cephean/respectrum"
#define PLUG_CLASS respectrum

#define PLUG_CONTROLS (6)
#define PLUG_INPUTS (1)
#define PLUG_OUTPUTS (1)

//======================================================

//class definition
class respectrum : public plugin
{
public:
	respectrum(float setFs);
	~respectrum();

	void init() override;
	void step(int len) override;
	void deactivate() override;

private:
	//ports -----------------
	fcontrol gateLow{ this, 0 };
	fcontrol gateHigh{ this, 1 };
	fcontrol timbre{ this, 2 };
	fcontrol attack{ this, 3 };
	fcontrol time{ this, 4 };
	fcontrol level{ this, 5 };

	input in{ this, 6 };

	output out{ this, 7 };

	//constants -------------
	const float slewdur = 0.1f;			//seconds

	const float fmin = 20.0f;			//Hz analysis DC reject
	const float fpre = 2000.0f;			//Hz preemphasis
	const float gpre = -20.0f;			//dB preemphasis

	const float fnoisemin = 200.0f;		//Hz noise curve low end
	const float fnoisemax = 4000.0f;	//Hz noise curve high end

	const float Noct = 1.0f / 1.0f;		//octave smoothing
	const float fsmomin = 200.0f;		//Hz min smoothing frequency
	const float fsmomax = 500.0f;		//Hz max smoothing frequency	
	const float sbdet = 0.01f;			//seconds subband env detect time
	const float sbhold = 0.04f;			//seconds subband hold time
	const float sbrel = 0.01f;			//seconds subband release time

	const float mincut = 0.0f;			//cut depth
	const float maxboost = 10.0f;		//boost depth

	const int maxlen = 128;				//maximum block size
	const int npts0 = 512;				//npts at 48k nominal

	//state etc.-------------
	rebuffer rebuff;

	fof preemph; //preemphasis filter
	fof deemph; //deemphasis filter

	const int npts; //effective npts at the current sampling rate
	const int K; //effective num subbands at the current sampling rate

	fvec wdw; //FFT window
	fft F; //FFT object
	vrel fenv; //initial spectral envelope detection
	vhold4 fhold; //spectral envelope hold
	vrel frel; //spectral envelope release
	specsmooth fsmo; //spectral log smoothing

	float vshift = 1.0f; //current shift multiplier value
	interpolator shift; //timbre shifter

	slewed<float> noiseLF; //low frequency noise floor in dB
	slewed<float> noiseHF; //high frequency noise floor in dB

	slewed<float> gattsus; //linear gain on net attack + sustain envelope
	slewed<float> gsus; //linear gain on sustain only envelope
	vatt fsustain1, fsustain2; //sustain split filter

	vsmo fhsmo; //spectral filter smoothing
	fastfir fir; //filter application

	flatbuffer ybuff; //filter buffer
	flatbuffer xbuff; //analysis buffer
	fvec X; //input spectral magnitude
	fvec Hmag; //spectral filter magnitude
	cfvec Hang; //spectral filter magnitude
	fvec Hdeemph; //deemphasis magnitude response
	fvec Hsupp; //suppression noise floor

	fvec fscratch1, fscratch2, fscratch3; //subband scratch vectors

	gain outLevel;

	//helper functions ------
	void update(int len = 0);
};
