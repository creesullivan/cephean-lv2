
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
	icontrol cmode{ this, 0 }; //0 - byp, 1 - 1x, 2 - 2x, 3 - 4x, 4 - 8x
	fcontrol ctime{ this, 1 };
	fcontrol cpitch{ this, 2 };
	fcontrol crate{ this, 3 };
	fcontrol camount{ this, 4 };
	fcontrol crolloff{ this, 5 };
	fcontrol clevel{ this, 6 };

	input in{ this, 7 };

	output out{ this, 8 };

	//constants -------------
	const float slewdur = 0.1f;			//seconds
	const float bypassdur = 0.01f;		//seconds
	const int maxlen = 128;				//maximum block size

	const float sepDelay = 0.005f;		//fixed minimum delay in sec

	const float minDelay = 0.01f;		//smallest maximum delay range in sec
	const float maxDelay = 0.1f;		//greatest maximum delay range in sec

	const float minShift = 0.1f;		//smallest maximum shift in +/- semitones
	const float maxShift = 0.5f;		//greatest maximum shift in semitones

	const float minLevel = 0.0f;		//smallest maximum level adjust in +/- dB
	const float maxLevel = 6.0f;		//greatest maximum level adjust in +/- dB

	const float maxPeriod = 4.0f;		//maximum period (rate inverse) in sec
	const float minPeriod = 0.5f;		//minimum period (rate inverse) in sec

	//state etc.-------------
	rebuffer rebuff;
	bypass byp;

	fvec inscr;
	fvec delscr;
	fvec vscr;
	fvec gscr;
	fvec gtotscr;

	pow2buffer buff;
	const float del0;

	slewed<float> crossoverFreq;
	sof lpf, hpf; //bass crossover

	slewed<float> amountGain; //not a gain object so we can multiply by makeup vector
	gain outputGain;

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
		//periodMin/Max - sample range of each period
		void set(float pitchMax, float delayMax, float levelMax,
			unsigned int periodMin, unsigned int periodMax);

		//note that the g trajectory should be normalized downstream to ensure
		//steady power across multiple voices, ie gnorm = 1/sqrt(sum_n{g_n})
		void stepBlock(float* del, float* g, int len);

		//clears current state to random settings within the current range
		void clear();

	private:
		//state ---------------
		float n; //current fractional delay position
		int l; //current interpolation step index
		int L; //current interpolation sample length
		float pnext; //next pitch shift step target
		float plast; //previous pitch shift step target
		float gnext; //next gain multiplier target
		float glast; //previous gain multiplier target

		//bounds --------------
		float pmax; //absolute max +/- pitch shift applied on update
		const float nmin = 0.0f; //absolute min fractional sample delay
		float nmax; //absolute max fractional sample delay
		float Gmax; //absolute max +/- gain in dB
		int Lmin; //minimum possible interpolation period
		int Lmax; //maximum possible interpolation period

		//helper functions --------------
		void update();
	};

	int Nvoice = 8; //active number of voices not including dry
	voice v[8];

};
