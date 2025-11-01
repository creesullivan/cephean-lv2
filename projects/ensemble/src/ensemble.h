
#pragma once

//--------------------------------------------------
// Defines the module's audio processing behavior
//--------------------------------------------------

#include "cephean-lv2.h"
#include "cephean-dsp.h"

using namespace cephean;

//======================================================

//config definitions
#define PLUG_URI "https://github.com/creesullivan/cephean/ensemble"
#define PLUG_CLASS ensemble

#define PLUG_CONTROLS (7)
#define PLUG_INPUTS (1)
#define PLUG_OUTPUTS (1)

//======================================================

//class definition
class ensemble : public plugin
{
public:
	ensemble(float setFs);
	~ensemble();

	void init() override;
	void step(int len) override;
	void deactivate() override;

private:
	//ports -----------------
	icontrol mode{ this, 0 };
	fcontrol time{ this, 1 };
	fcontrol pitch{ this, 2 };
	fcontrol level{ this, 3 };
	fcontrol rate{ this, 4 };
	fcontrol dry{ this, 5 };
	fcontrol wet{ this, 6 };

	input in{ this, 7 };

	output out{ this, 8 };

	//constants -------------
	const float slewdur = 0.1f;			//seconds
	const int maxlen = 128;				//maximum block size

	const float minDelay = 0.01f;		//smallest maximum delay in sec
	const float maxDelay = 0.1f;		//greatest maximum delay in sec

	const float minShift = 0.1f;		//smallest maximum shift in +/- semitones
	const float maxShift = 1.0f;		//greatest maximum shift in semitones

	const float minLevel = 0.0f;		//smallest maximum level adjust in +/- dB
	const float maxLevel = 6.0f;		//greatest maximum level adjust in +/- dB

	const float maxPeriod = 4.0f;		//maximum period (rate inverse) in sec
	const float minPeriod = 0.5f;		//minimum period (rate inverse) in sec

	//state etc.-------------
	rebuffer rebuff;

	fvec inscr;
	pow2buffer delay; //TESTING the fast pow2 circular buffer with a simple time delay
	//circbuffer delay;

	//helper functions ------
	void update(int len = 0);
	void clear();

	//Object that generates the delay and normalized gain trajectory for
	//a single chorus voice based on the rate and time/pitch/level boundaries
	class voice
	{
	public:
		voice();
		~voice();

		//pitchMax - samples/sample, [-pitchMax, pitchMax]
		//delayMax - samples, [0, delayMax]
		//levelMax - dB >= 0, [-levelMax, levelMax]
		void setRange(float pitchMax, uint16_t delayMax, float levelMax);

		//note that the g trajectory should be normalized downstream to ensure
		//steady power across multiple voices, ie gnorm = 1/sqrt(sum_n{g_n})
		void stepBlock(uint16_t* del, float* g, int len);

		//LEFT OFF HERE <------- get this bit working, at least a little, then
		// use it to test the speed of a getNN() with pow2buffer vs circbuffer

	private:

	};
};
