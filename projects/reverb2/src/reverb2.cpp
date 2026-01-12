
//--------------------------------------------------
// Defines the module's audio processing behavior
//--------------------------------------------------

#include "reverb2.h"

reverb2::reverb2(double setFs) :
	plugin(PLUG_CONTROLS, PLUG_INPUTS, PLUG_OUTPUTS, (float)setFs),
	rebuff(maxlen),
	scolor((int)roundf(slewdur * Fs)),
	srolloff((int)roundf(slewdur * Fs)),

	verb((int)roundf(slewdur * Fs), maxlen, getNominalUSR(Fs)),
	bpf(maxlen),
	gainWet((int)roundf(slewdur * Fs)),
	gainDry((int)roundf(slewdur * Fs)),

	inscr(maxlen)
{}
reverb2::~reverb2() {}

void reverb2::init()
{
	rebuff.reset();

	clear();
	commit();
	update();
}

void reverb2::deactivate() {}

void reverb2::update(int samples)
{
	bool force = (samples == 0); //if no processing samples, we are initializing

	//tail --------------------------------------
	if (tail.isnew() || force) {
		float vtail = 0.0f;
		if (tail <= 50.0f) {
			vtail = boundlinmap(tail, 0.0f, 50.0f, 0.0f, 0.5f);
		}
		else {
			vtail = boundlogmap(tail, 50.0f, 100.0f, 0.5f, 5.0f);
		}
		verb.setDecay(vtail * Fs, force);

		tail.clearnew();
	}

	//color -----------------------------------------
	if (color.isnew() || force) {
		float vcolor = 0.0f;
		if (color <= 75.0f) {
			vcolor = boundlinmap(color, 0.0f, 75.0f, 1000.0f, 5000.0f);
		}
		else {
			vcolor = boundlogmap(color, 75.0f, 100.0f, 5000.0f, 20000.0f);
		}
		scolor.target(vcolor / (Fs / 2.0f), force);

		color.clearnew();
	}
	if (scolor.slew(samples) || force) {
		bpf[1].setCoefs(filt::lowpass2(scolor.get(), 0.707f,
			20.0f / (Fs / 2.0f), 20000.0f / (Fs / 2.0f)), force);
	}

	//rolloff -----------------------------------------
	if (rolloff.isnew() || force) {
		float vrolloff = 0.0f;
		if (rolloff <= 25.0f) {
			vrolloff = boundlogmap(rolloff, 0.0f, 25.0f, 60.0f, 200.0f);
		}
		else{
			vrolloff = boundlogmap(rolloff, 25.0f, 100.0f, 200.0f, 2000.0f);
		}
		srolloff.target(vrolloff / (Fs / 2.0f), force);

		rolloff.clearnew();
	}
	if (srolloff.slew(samples) || force) {
		bpf[0].setCoefs(filt::highpass2(srolloff.get(), 0.707f), force);
	}

	//energy -----------------------------------------
	if (energy.isnew() || force) {
		float venergy = boundlinmap(energy, 0.0f, 100.0f, 0.0f, 12.0f); //max cents
		verb.setRate(expf(constants.two2nats * venergy / 1200.0f) - 1.0f, force);

		energy.clearnew();
	}

	//wet ----------------------------------------
	if (wet.isnew() || force) {
		float vwet = 0.0f;
		if (wet <= 25.0f) {
			vwet = boundlinmap(wet, 0.0f, 25.0f, 0.0f, 0.1f);
		}
		else {
			vwet = boundlogmap(wet, 25.0f, 100.0f, 0.1f, 10.0f);
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

void reverb2::clear()
{
	verb.clear();
	bpf.clear();
}

void reverb2::step(int len)
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

		safeInput(pin, inscr, curlen);

		//apply the algorithm, forming the wet path
		gainWet.stepBlock(inscr, pout, curlen);
		verb.stepBlock(pout, pout, curlen);
		bpf.stepBlock(pout, pout, curlen); //bandpass last for smoothness and correctness
		
		//mix dry back in
		gainDry.stepBlock(inscr, inscr, curlen);
		vadd(pout, inscr.ptr(), pout, curlen);

		if (safeOutput(pout, curlen)) {
			clear();
		}

		pin += curlen;
		pout += curlen;
		len -= curlen;
	}
}
