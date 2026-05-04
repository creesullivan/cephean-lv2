
//--------------------------------------------------
// Defines the module's audio processing behavior
//--------------------------------------------------

#include "compressor.h"

compressor::compressor(float setFs) :
	plugin(PLUG_CONTROLS, PLUG_INPUTS, PLUG_OUTPUTS, setFs),
	rebuff(maxlen),
	inscr(maxlen), scr1(maxlen), scr2(maxlen), scr3(maxlen),

	slewTone((int)roundf(slewdur * Fs)),
	prefilt(maxlen),
	postfilt(maxlen),

	rmsdetect(1, (double)Fs, (int)roundf(Trms * Fs), maxlen),
	release(Tatt * Fs, Trel * Fs),

	curve(),

	driveGain((int)roundf(slewdur * setFs)),
	dryGain((int)roundf(slewdur * setFs)),
	wetGain((int)roundf(slewdur * setFs))
{
	//design fixed LMA, release is pre-designed from constructor
	rmsdetect.setLength((int)roundf(Trms * Fs), true);

	//design fixed curve
	curve.setKneeWidth(knee, true);
}
compressor::~compressor() {}

void compressor::init()
{
	rebuff.reset();

	clear();
	commit();
	update();
}

void compressor::deactivate() {}

void compressor::update(int samples)
{
	bool force = (samples == 0); //if no processing samples, we are initializing

	//drive --------------------------------------
	if(drive.isnew() || force){
		float vdrive = boundlogmap(drive, 0.0f, 100.0f, 1.0f, 100.0f);
		driveGain.setLevel(vdrive, force);
		
		drive.clearnew();
	}

	//tone --------------------------------------
	if (tone.isnew() || force) {
		float vtone = boundlogmap(tone, -100.0f, 100.0f, 0.5f, 2.0f);
		slewTone.target(vtone, force);
		
		tone.clearnew();
	}
	if (slewTone.slew(samples) || force) {
		fof::coefs ctone = filt::balance1(tonefreq / (Fs / 2.0f), slewTone.get());
		prefilt.setCoefs(filt::inverse(ctone), force);
		postfilt.setCoefs(ctone, force);
	}

	//dry --------------------------------------
	if (dry.isnew() || force) {
		float vdry = boundlinmap(dry, 0.0f, 100.0f, 0.0f, 2.0f);
		dryGain.setLevel(vdry, force);

		dry.clearnew();
	}

	//wet --------------------------------------
	if (wet.isnew() || force) {
		float vwet = 0.0f;
		if (wet < 25.0f) {
			vwet = boundlinmap(wet, 0.0f, 25.0f, 0.0f, 0.125f);
		}
		else {
			vwet = boundlogmap(wet, 25.0f, 100.0f, 0.125f, 1.0f);
		}
		wetGain.setLevel(vwet, force);

		wet.clearnew();
	}
}

void compressor::clear()
{
	prefilt.clear();
	postfilt.clear();
}

void compressor::step(int len)
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

		//preemph and drive gain
		prefilt.stepBlock(inscr, inscr, curlen);
		driveGain.stepBlock(inscr, scr1, curlen);

		//envelope detection
		vmult(scr1.ptr(), scr1.ptr(), scr2.ptr(), curlen); //x*x
		rmsdetect.stepBlock(scr2, scr3, curlen); //power
		vmult(scr2.ptr(), crest, scr2.ptr(), curlen); //scale by power crest factor
		vmax(scr2.ptr(), scr3.ptr(), scr2.ptr(), curlen); //power with transient correction
		vsqrt(scr2.ptr(), scr2.ptr(), curlen); //RMS envelope
		release.stepBlock(scr2, scr2, curlen); //final asym smoothed envelope

		//gain computation and application
		curve.stepGainBlock(scr2, scr2, curlen);
		vmult(scr1.ptr(), scr2.ptr(), scr1.ptr(), curlen);

		//dry/wet mix and deemph
		dryGain.stepBlock(inscr.ptr(), inscr.ptr(), curlen);
		wetGain.stepBlock(scr1.ptr(), scr1.ptr(), curlen);
		vadd(inscr.ptr(), scr1.ptr(), pout, curlen);
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
