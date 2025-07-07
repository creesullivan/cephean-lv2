
//--------------------------------------------------
// Defines the module's audio processing behavior
//--------------------------------------------------

#include "reverb.h"

reverb::reverb(double setFs, int setSeed) :
	plugin(PLUG_CONTROLS, PLUG_INPUTS, PLUG_OUTPUTS, (float)setFs),
	rebuff(maxlen),
	scolor((int)roundf(slewdur * (float)setFs)),

	gainWet((int)roundf(slewdur* (float)setFs)),
	gainDry((int)roundf(slewdur* (float)setFs)),
	verb((int)roundf(slewdur * Fs), maxlen, getNominalUSR(Fs)),
	bpf(maxlen),

	dryscr(maxlen)
{
	
	if(setSeed >= 0) { //randomize directly from seed when testing and ignoring seed parameter
		int* pdel = verb.getDelayPointer();
		srand((unsigned int)setSeed);
		for (int m = 0; m < 16; ++m) {
			pdel[m] = (int)(480 + (4800 - 480) * randf());
		}
		int del[16];
		vcopy(pdel, del, 16);
		
		float bpdummy = 0.0f;

		vmult(pdel, getNominalUSR(Fs), pdel, 16); //adjust for sampling rate
	}
	
	//design fixed highpass for rumble rejection
	bpf[0].setCoefs(filt::highpass2(60.0f / (Fs / 2.0f), 0.707f), true);
}
reverb::~reverb() {}

void reverb::init()
{
	rebuff.reset();

	verb.clear();
	bpf.clear();

	commit();
	update();
}

void reverb::deactivate() {}

void reverb::update(int samples)
{
	bool force = (samples == 0); //if no processing samples, we are initializing

	//seed ------------------------------------------
	if (seed.isnew() || force) {
		unsigned int vseed = (unsigned int)bound(seed.get(), 0, 15);
		verb.setSeed(vseed);
		seed.clearnew();
	}

	//duration --------------------------------------
	if (duration.isnew() || force) {
		float vdur = 0.0f;
		if (duration <= 25.0f) {
			vdur = boundlinmap(duration, 0.0f, 25.0f, 0.0f, 0.5f);
		}
		else if (duration <= 75.0f) {
			vdur = boundlogmap(duration, 25.0f, 75.0f, 0.5f, 2.0f);
		}
		else {
			vdur = boundlogmap(duration, 75.0f, 100.0f, 2.0f, 10.0f);
		}
		verb.setDecay(vdur * Fs, force);

		duration.clearnew();
	}

	//color -----------------------------------------
	if (color.isnew() || force) {
		float vcolor = 0.0f;
		if (color <= 75.0f) {
			vcolor = boundlogmap(color, 0.0f, 75.0f, 1000.0f, 4000.0f);
		}
		else {
			vcolor = boundlogmap(color, 75.0f, 100.0f, 4000.0f, 20000.0f);
		}
		scolor.target(vcolor / (Fs / 2.0f), force);

		color.clearnew();
	}
	if (scolor.slew(samples) || force) {
		bpf[1].setCoefs(filt::lowpass2(scolor.get(), 0.707f,
			20.0f / (Fs / 2.0f), 20000.0f / (Fs / 2.0f)), force);
	}

	//wet ----------------------------------------
	if (wet.isnew() || force) {
		float vwet = 0.0f;
		if (wet <= 50.0f) {
			vwet = boundlinmap(wet, 0.0f, 50.0f, 0.0f, 1.0f);
		}
		else {
			vwet = boundlogmap(wet, 50.0f, 100.0f, 1.0f, 10.0f);
		}
		gainWet.setLevel(vwet, force);

		wet.clearnew();
	}

	//dry ----------------------------------------
	if (dry.isnew() || force) {
		float vdry = 0.0f;
		if (dry <= 50.0f) {
			vdry = boundlinmap(dry, 0.0f, 50.0f, 0.0f, 1.0f);
		}
		else {
			vdry = boundlogmap(dry, 50.0f, 100.0f, 1.0f, 10.0f);
		}
		gainDry.setLevel(vdry, force);

		dry.clearnew();
	}

}

void reverb::step(int len)
{
	//grab I/O pointers
	const float* pin = in;
	float* pout = out;

	//commit control values at the top of the block
	commit();

	//run processing loop
	int curlen = 0;
	while (len > 0) {
		curlen = rebuff.next(len);

		//update to slewed parameters
		update(curlen);

		//apply the algorithm, forming the wet path
		gainWet.stepBlock(pin, pout, curlen);
		bpf.stepBlock(pout, pout, curlen);
		verb.stepBlock(pout, pout, curlen); //actual verb last for smoothness
		
		//mix dry back in
		gainDry.stepBlock(pin, dryscr.ptr(), curlen);
		vadd(pout, dryscr.ptr(), pout, curlen);

		pin += curlen;
		pout += curlen;
		len -= curlen;
	}
}
