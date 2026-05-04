
//--------------------------------------------------
// Defines the module's audio processing behavior
//--------------------------------------------------

#include "tapdelay.h"

tapdelay::tapdelay(float setFs) :
	plugin(PLUG_CONTROLS, PLUG_INPUTS, PLUG_OUTPUTS, setFs),
	rebuff(maxlen),
	inscr(maxlen),
	fbscr(maxlen),

	slewRolloff(maxlen),

	//buff((int)ceilf(maxdur * Fs + maxlen)),
	buff(maxlen),

	fbGain((int)roundf(slewdur * Fs)),
	delay((int)roundf(tslewdur * Fs)),
	hpf((int)roundf(slewdur * Fs)),
	lpf((int)roundf(slewdur * Fs)),

	dryGain((int)roundf(slewdur * setFs)),
	wetGain((int)roundf(slewdur * setFs))
{}
tapdelay::~tapdelay() {}

void tapdelay::init()
{
	rebuff.reset();

	clear();
	commit();
	update();
}

void tapdelay::deactivate() {}

void tapdelay::update(int samples)
{
	bool force = (samples == 0); //if no processing samples, we are initializing

	//tempo & mode --------------------------------------
	if (tempo.isnew() || mode.isnew() || force) {
		int vmode = mode.get();
		float modemult = 1.0f;
		switch (vmode)
		{
		case 0:
			modemult = 1.0f;
			break;
		case 1:
			modemult = 3.0f / 4.0f;
			break;
		case 2:
			modemult = 2.0f / 3.0f;
			break;
		case 3:
			modemult = 1.0f / 2.0f;
			break;
		case 4:
			modemult = 1.0f / 3.0f;
			break;
		case 5:
			modemult = 1.0f / 4.0f;
			break;
		}
		int vtime = bound((int)roundf((60.0f / tempo.get()) * modemult * Fs),
			maxlen, buff.size() - maxlen);
		delay.target(vtime, force);

		tempo.clearnew();
		mode.clearnew();
	}

	//feedback --------------------------------------
	if (feedback.isnew() || force) {
		float vfeedback = boundlinmap(feedback, 0.0f, 100.0f, 0.0f, 0.99f);
		fbGain.setLevel(vfeedback, force);

		feedback.clearnew();
	}

	// color  --------------------------------------
	if (color.isnew() || force) {
		float vcolor = boundlogmap(color, 0.0f, 100.0f, 400.0f, 20000.0f);
		lpf.setLowpass(vcolor / (Fs / 2.0f), force);

		color.clearnew();
	}

	//rolloff --------------------------------------
	if (rolloff.isnew() || force) {
		float vrolloff = 0.0f;
		if (rolloff.get() < 25.0f) {
			vrolloff = boundlinmap(rolloff, 0.0f, 25.0f, 20.0f, 200.0f);
		}
		else {
			vrolloff = boundlogmap(rolloff, 25.0f, 100.0f, 200.0f, 2000.0f);
		}
		slewRolloff.target(vrolloff, force);

		rolloff.clearnew();
	}
	if (slewRolloff.slew(samples) || force) {
		hpf.setCoefs(filt::highpass2(slewRolloff.get() / (Fs / 2.0f), 0.707f,
			20.0f / (Fs / 2.0f), 20000.0f / (Fs / 2.0f)), force);
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

	// dry -------------------------------------------
	if (dry.isnew() || force) {
		float vdry = boundlinmap(dry, 0.0f, 100.0f, 0.0f, 2.0f);
		dryGain.setLevel(vdry, force);

		dry.clearnew();
	}

	// wet -------------------------------------------
	if (wet.isnew() || force) {
		float vwet = boundlinmap(wet, 0.0f, 100.0f, 0.0f, 2.0f);
		wetGain.setLevel(vwet, force);

		wet.clearnew();
	}
}

void tapdelay::clear() {}

void tapdelay::step(int len)
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
		
		safeInput(pin, inscr, curlen);

		//retrieve from feedback delay
		delay.swap();
		if (delay.check()) { //crossfading
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

		//apply feedback gain and color
		fbGain.stepBlock(pout, pout, curlen);
		lpf.stepBlock(pout, pout, curlen);

		//put back to feedback buffer
		wetGain.stepBlock(inscr, fbscr, curlen); //wet gain here to fade out delays
		vadd(fbscr.ptr(), pout, fbscr.ptr(), curlen);
		buff.put(fbscr.ptr(), curlen);

		//final mix
		dryGain.stepBlock(inscr, inscr, curlen);
		hpf.stepBlock(pout, pout, curlen); //highpass the echoes
		vadd(pout, inscr.ptr(), pout, curlen); //mix dry with wet
		
		if (safeOutput(pout, curlen)) {
			clear();
		}

		//step to next sub-block
		pin += curlen;
		pout += curlen;
		len -= curlen;
	}
}
