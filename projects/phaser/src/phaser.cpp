
//--------------------------------------------------
// Defines the module's audio processing behavior
//--------------------------------------------------

#include "phaser.h"

phaser::phaser(float setFs) :
	plugin(PLUG_CONTROLS, PLUG_INPUTS, PLUG_OUTPUTS, setFs),
	rebuff(maxlen),
	ph(2.0f*constants.pi),
	pcore((int)roundf(slewdur* Fs), maxlen),

	dph((int)roundf(slewdur* Fs)),
	del1((int)roundf(slewdur * Fs)),
	del2((int)roundf(slewdur * Fs)),
	dryGain((int)roundf(slewdur* Fs)),
	wetGain((int)roundf(slewdur* Fs)),

	inscr(maxlen)
{}
phaser::~phaser() {}

void phaser::init()
{
	rebuff.reset();

	clear();
	commit();
	update();
}

void phaser::deactivate() {}

void phaser::update(int samples)
{
	bool force = (samples == 0); //if no processing samples, we are initializing

	//depth ------------------------------
	if (depth.isnew() || force) {
		float fbval = 0.0f;
		float depval = 0.0f;
		if (depth <= 25.0f) {
			fbval = 0.0f;
			depval = boundlinmap(depth, 0.0f, 25.0f, 0.0f, 0.5f);
		}
		else if (depth <= 50.0f) {
			fbval = boundlinmap(depth, 25.0f, 50.0f, 0.0f, 0.5f);
			depval = boundlinmap(depth, 25.0f, 50.0f, 0.5f, 1.0f);
		}
		else {
			fbval = boundlogmap(depth, 50.0f, 100.0f, 0.5f, 0.9f);
			depval = 1.0f;
		}
		pcore.setDepth(depval, -fbval, force);

		depth.clearnew();
	}

	//stop1 ------------------------------
	if (stop1.isnew() || force) {
		float val = boundlogmap(stop1, 0.0f, 100.0f, 400.0f, 8000.0f);
		del1.target((Fs / val) / (float)order, force);

		stop1.clearnew();
	}

	//stop2 ------------------------------
	if (stop2.isnew() || force) {
		float val = boundlogmap(stop2, 0.0f, 100.0f, 400.0f, 8000.0f);
		del2.target((Fs / val) / (float)order, force);

		stop2.clearnew();
	}

	//rate ------------------------------
	if (rate.isnew() || force) {
		float val = boundlogmap(rate, 0.0f, 100.0f, 0.1f, 10.0f);
		dph.target((2.0f * constants.pi) / (Fs / val), force);

		rate.clearnew();
	}

	//tone & level ----------------------
	if (tone.isnew() || level.isnew() || force) {
		float tval = boundlinmap(tone, 0.0f, 100.0f, 0.0f, -0.5f);
		float lval = boundlinmap(level, 0.0f, 100.0f, 0.0f, 2.0f);
		dryGain.setLevel(lval * tval, force);
		wetGain.setLevel(lval * (1.0f - fabsf(tval)), force);

		tone.clearnew();
		level.clearnew();
	}
}

void phaser::clear()
{
	pcore.clear();
	ph = 0.0f;
}

void phaser::step(int len)
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

		//step the oscillator
		dph.slew(curlen);
		del1.slew(curlen);
		del2.slew(curlen);
		dph.slew(curlen);
		ph += (dph * curlen);
		float curalpha = boundlogmap(triwave(ph), -1.0f, 1.0f, del1, del2);
		curalpha = pcore.getAlphaForDelay(curalpha);

		//apply the warped buffer with feedback
		pcore.stepBlock(inscr, curalpha, pout, curlen);

		//mix dry/wet to generate comb filtered response
		dryGain.stepBlock(inscr, inscr, curlen);
		wetGain.stepBlock(pout, pout, curlen);
		vadd(inscr.ptr(), pout, pout, curlen);

		if (safeOutput(pout, curlen)) {
			clear();
		}

		//step to next sub-block
		pin += curlen;
		pout += curlen;
		len -= curlen;
	}

	
}
