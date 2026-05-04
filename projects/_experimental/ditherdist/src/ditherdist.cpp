
//--------------------------------------------------
// Defines the module's audio processing behavior
//--------------------------------------------------

#include "ditherdist.h"

ditherdist::ditherdist(float setFs) :
	plugin(PLUG_CONTROLS, PLUG_INPUTS, PLUG_OUTPUTS, setFs),
	rebuff(maxlen),
	inscr(maxlen),

	slewfilt(3, (int)roundf(slewdur * Fs)),
	prefilt(maxlen),

	dither(maxlen),
	waveshaper((int)roundf(slewdur * Fs)),

	postfilt(maxlen),
	dryGain((int)roundf(slewdur * Fs)),
	wetGain((int)roundf(slewdur * Fs))
{
	/*
	//TESTING with harcoded parameters
	waveshaper.set(100.0f, 0.0f, 0.0f, true);

	fof::coefs hpf1 = filt::highpass1(150.0f / (Fs / 2.0f));
	fof::coefs tone1 = filt::balance1(1000.0f / (Fs / 2.0f), 1.7f);
	sof::coefs lpf2 = filt::lowpass2(10000.0f / (Fs / 2.0f), 0.707f);
	prefilt[0].setCoefs(filt::merge(hpf1, tone1), true);
	prefilt[1].setCoefs(lpf2, true);
	postfilt[0].setCoefs(filt::merge(hpf1, filt::inverse(tone1)), true);
	postfilt[1].setCoefs(lpf2, true);

	dryGain.setLevel(1.0f, true);
	wetGain.setLevel(0.25f, true); //for debugging safety
	*/
}
ditherdist::~ditherdist() {}

void ditherdist::init()
{
	rebuff.reset();

	clear();
	commit();
	update();
}

void ditherdist::deactivate() {}

void ditherdist::update(int samples)
{
	bool force = (samples == 0); //if no processing samples, we are initializing

	//drive, bias, gate --------------------------------------
	if (drive.isnew() || bias.isnew() || gate.isnew() || force) {
		float vdrive = 0.0f;
		if (drive <= 25.0f) {
			vdrive = boundlogmap(drive, 0.0f, 25.0f, 1.0f, 3.0f);
		}
		else if(drive <= 50.0f){
			vdrive = boundlogmap(drive, 25.0f, 50.0f, 3.0f, 10.0f);
		}
		else if (drive <= 75.0f) {
			vdrive = boundlogmap(drive, 50.0f, 75.0f, 10.0f, 100.0f);
		}
		else {
			vdrive = boundlogmap(drive, 75.0f, 100.0f, 100.0f, 10000.0f);
		}

		float vbias = boundlinmap(bias, -100.0f, 100.0f, -1.0f, 1.0f);

		float vgate = 0.0f;
		if (gate <= 25.0f) {
			vgate = boundlinmap(gate, 0.0f, 25.0f, 0.0f, 0.0001f);
		}
		else {
			vgate = boundlogmap(gate, 25.0f, 100.0f, 0.0001f, 1.0f);
		}

		waveshaper.set(vdrive, vbias, vgate, force);

		drive.clearnew();
		bias.clearnew();
		gate.clearnew();
	}

	//low, high, tone ---------------------------------------
	if (low.isnew() || high.isnew() || tone.isnew() || force) {
		float vlow = boundlogmap(low, 0.0f, 100.0f, 30.0f, 480.0f);
		float vhigh = boundlinmap(high, 0.0f, 100.0f, 1000.0f, 10000.0f);
		float vtone = boundlogmap(tone, -100.0f, 100.0f, 0.5f, 2.0f);

		float* filtparams = slewfilt.target();
		filtparams[0] = vlow;
		filtparams[1] = vhigh;
		filtparams[2] = vtone;
		if (force) {
			slewfilt.converge();
		}

		low.clearnew();
		high.clearnew();
		tone.clearnew();
	}
	if (slewfilt.slew(samples) || force) { //design the pre and post filters given the slewed parameters
		fof::coefs hpf1 = filt::highpass1(slewfilt.get()[0] / (Fs / 2.0f));
		sof::coefs lpf2 = filt::lowpass2(slewfilt.get()[1] / (Fs / 2.0f), 0.707f);
		fof::coefs tone1 = filt::balance1(1000.0f / (Fs / 2.0f), slewfilt.get()[2]);

		prefilt[0].setCoefs(filt::merge(hpf1, filt::inverse(tone1)), force);
		prefilt[1].setCoefs(lpf2, force);
		postfilt[0].setCoefs(filt::merge(hpf1, tone1), force);
		postfilt[1].setCoefs(lpf2, force);
	}

	//dry --------------------------------------
	if(dry.isnew() || force){
		float vdry = boundlinmap(dry, 0.0f, 100.0f, 0.0f, 2.0f);
		dryGain.setLevel(vdry, force);
		
		dry.clearnew();
	}

	//wet --------------------------------------
	if (wet.isnew() || force) {
		float vwet = 0.0f;
		if (wet <= 25.0f) {
			vwet = boundlinmap(wet, 0.0f, 25.0f, 0.0f, 0.03f);
		}
		else {
			vwet = boundlinmap(wet, 25.0f, 100.0f, 0.03f, 1.0f);
		}
		wetGain.setLevel(vwet, force);

		wet.clearnew();
	}
}

void ditherdist::clear()
{
	dither.clear();
}

void ditherdist::step(int len)
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

		prefilt.stepBlock(inscr, inscr, curlen);

		//dither.apply(inscr, pout, curlen);
		vcopy(inscr.ptr(), pout, curlen);
		waveshaper.stepBlock(pout, pout, curlen);
		//dither.invert(pout, pout, curlen);

		dryGain.stepBlock(inscr, inscr, curlen); //dry path totally avoids noise from the dithering operation
		wetGain.stepBlock(pout, pout, curlen);
		vadd(inscr.ptr(), pout, pout, curlen); //mix dry/wet

		postfilt.stepBlock(pout, pout, curlen);
		
		if (safeOutput(pout, curlen)) {
			clear();
		}

		//step to next sub-block
		pin += curlen;
		pout += curlen;
		len -= curlen;
	}
}
