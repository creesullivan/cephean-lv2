
//--------------------------------------------------
// Defines the module's audio processing behavior
//--------------------------------------------------

#include "fbdelay.h"

fbdelay::fbdelay(float setFs, int setSeed) :
	plugin(PLUG_CONTROLS, PLUG_INPUTS, PLUG_OUTPUTS, setFs),
	rebuff(maxlen),
	buff((int)ceilf(maxdur* Fs + maxlen)),
	fbGain((int)roundf(slewdur* setFs)),
	delay((int)roundf(tslewdur* Fs)),
	hpf((int)roundf(slewdur* Fs)),
	echo((int)roundf(slewdur * Fs), maxecho * getNominalUSR(Fs)),
	lpf((int)roundf(slewdur* Fs)),
	wetGain((int)roundf(slewdur* setFs))
{
	echo.setSeed(15); //fixed
	/* only valid when using miniverb, not diffuser!
	//set the echo delay vector
	int* pdel = echo.getDelayPointer();
	if (setSeed < 0) { //load the preset vector that sounds great
		echo.setSeed(8);
	}
	else { //randomize directly from seed
		srand((unsigned int)setSeed);
		int viewdels[9];
		for (int m = 0; m < 9; ++m) {
			viewdels[m] = (int)(minecho + (maxecho - minecho) * randf());
			pdel[m] = viewdels[m];
		}
		float bpdummy = 0.0f;
	}
	vmult(pdel, getNominalUSR(Fs), pdel, 9); //adjust for sampling rate
	*/
}
fbdelay::~fbdelay() {}

void fbdelay::init()
{
	rebuff.reset();
	buff.clear();
	echo.clear();

	commit();
	update();
}

void fbdelay::deactivate() {}

void fbdelay::update(int samples)
{
	bool force = (samples == 0); //if no processing samples, we are initializing

	//time --------------------------------------
	if (time.isnew() || force) {
		float vtime = 0.0f;
		if (time <= 25.0f) { //slapback territory
			vtime = boundlinmap(time, 0.0f, 25.0f, 0.03f, 0.15f);
		}
		else if(time <= 50.0f) {
			vtime = boundlogmap(time, 25.0f, 50.0f, 0.15f, 0.3f);
		}
		else {
			vtime = boundlogmap(time, 50.0f, 100.0f, 0.3f, 1.0f);
		}
		delay.target(max((int)roundf(vtime * Fs), maxlen), force);
		time.clearnew();
	}

	//feedback --------------------------------------
	if (feedback.isnew() || force) {
		float vfeedback = boundlinmap(feedback, 0.0f, 100.0f, 0.0f, 0.99f);
		fbGain.setLevel(vfeedback, force);
		feedback.clearnew();
	}

	// spread  --------------------------------------
	if(spread.isnew() || force){
		float vspread = boundlinmap(spread, 0.0f, 100.0f, 0.0f, 0.25f);
		echo.setDecay(vspread * Fs, force);
		spread.clearnew();
	}

	// wet -------------------------------------------
	if (wet.isnew() || force) {
		float vwet = 0.0f;
		if (wet <= 50.0f) {
			vwet = boundlinmap(wet, 0.0f, 50.0f, 0.0f, 0.5f);
		}
		else {
			vwet = boundlogmap(wet, 50.0f, 100.0f, 0.5f, 2.0f);
		}
		wetGain.setLevel(vwet, force);
		wet.clearnew();
	}

	//low --------------------------------------
	if (low.isnew() || force) {
		float vlow = 0.0f;
		if (low <= 50.0f) {
			vlow = boundlogmap(low, 0.0f, 50.0f, 20.0f, 400.0f);
		}
		else if(low <= 75.0f) {
			vlow = boundlogmap(low, 50.0f, 75.0f, 400.0f, 1500.0f);
		}
		else {
			vlow = boundlogmap(low, 75.0f, 100.0f, 1500.0f, 20000.0f);
		}
		hpf.setCoefs(filt::highpass2(vlow / (Fs / 2.0f), 0.707f,
			20.0f / (Fs / 2.0f), 20000.0f / (Fs / 2.0f)), force);
		low.clearnew();
	}

	//high --------------------------------------
	if (high.isnew() || force) {
		float vhigh = 0.0f;
		if (high <= 25.0f) {
			vhigh = boundlogmap(high, 0.0f, 25.0f, 20.0f, 1500.0f);
		}
		else if(high <= 50.0f) {
			vhigh = boundlogmap(high, 25.0f, 50.0f, 1500.0f, 4000.0f);
		}
		else {
			vhigh = boundlogmap(high, 50.0f, 100.0f, 4000.0f, 20000.0f);
		}
		sofcasc<2>::coefs setc;
		lpf.setCoefs(filt::lowpass4(setc, vhigh / (Fs / 2.0f), 0.707f, 0.0f,
			20.0f / (Fs / 2.0f), 20000.0f / (Fs / 2.0f)), force);
		high.clearnew();
	}
}

void fbdelay::step(int len)
{
	//grab I/O pointers
	const float* pin = in;
	float* pout = out;

	//commit control values at the top of the block
	commit();

	//run processing loop
	int curlen = 0;
	while (len > 0) {
		curlen = rebuff.next(len); //current sub-block length

		//update to slewed parameters
		update(curlen);

#if PLUG_TESTING_ECHO

		//apply echo with wet gain only for seed determination and debugging/tuning
		echo.stepBlock(pin, pout, curlen);
		//vmultaccum(pin, dryecho, pout, curlen);
		wetGain.stepBlock(pout, pout, curlen);

#else

		//retrieve from feedback delay
		delay.swap();
		if (delay.check()) {
			float g = 0.0f;
			float y = 0.0f;
			for (int n = 0; n < curlen; ++n) {
				delay.fade();
				g = delay.gain();
				y = g * buff.get(delay.current() - curlen, curlen, n);
				y += (1.0f - g) * buff.get(delay.last() - curlen, curlen, n);
				pout[n] = y;
			}
		}
		else { //no crossfade
			buff.get(pout, curlen, delay.current() - curlen);
		}

		//apply feedback gain
		fbGain.stepBlock(pout, pout, curlen);

		//put back to feedback buffer
		vadd(pout, pin, pout, curlen);
		buff.put(pout, curlen);

		//wet path processing
		vsub(pout, pin, pout, curlen); //remove dry component
		hpf.stepBlock(pout, pout, curlen); //highpass
		lpf.stepBlock(pout, pout, curlen); //lowpass
		echo.stepBlock(pout, pout, curlen); //add echo spread
		wetGain.stepBlock(pout, pout, curlen); //apply wet gain

		//final mix
		vadd(pout, pin, pout, curlen); //mix dry with wet

#endif

		//step to next sub-block
		pin += curlen;
		pout += curlen;
		len -= curlen;
	}
}
