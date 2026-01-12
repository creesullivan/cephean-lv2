
//--------------------------------------------------
// Defines the module's audio processing behavior
//--------------------------------------------------

#include "distortion3.h"

distortion3::distortion3(float setFs) :
	plugin(PLUG_CONTROLS, PLUG_INPUTS, PLUG_OUTPUTS, setFs),
	rebuff(maxlen),
	inscr(maxlen),

	slewTone((int)roundf(Fs * slewdur)),

	preTone(maxlen),

	envAttack(),
	envHold(0.0f),
	envRelease(),

	compress(maxlen),
	distort(maxlen, maxlen),

	postTone(maxlen),

	outGain((int)roundf(slewdur* Fs)),

	scr(maxlen)
{
	//envelope detection
	envAttack.setAttack(attdur * Fs);
	envAttack.setRelease(attdur * Fs);
	envHold.setQuarterLength((int)roundf(holddur * Fs / 4.0f));
	envRelease.setAttack(0.0f);
	envRelease.setRelease(reldur * Fs);
}
distortion3::~distortion3() {}

void distortion3::init()
{
	rebuff.reset();

	clear();
	commit();
	update();
}

void distortion3::deactivate() {}

void distortion3::update(int samples)
{
	bool force = (samples == 0); //if no processing samples, we are initializing

	//drive & shape --------------------------------------
	if (cdrive.isnew() || cshape.isnew() || force) {
		float vdrive = 0.0f;
		if (cdrive <= 25.0f) {
			vdrive = boundlinmap(cdrive, 0.0f, 25.0f, 0.01f, 0.1f); //true 0 drive unsupported
		}
		else {
			vdrive = boundlogmap(cdrive, 25.0f, 100.0f, 0.1f, 10.0f);
		}
		
		float vshape = 0.0f;
		if (cshape <= 25.0f) {
			vshape = boundlinmap(cshape, 0.0f, 25.0f, 0.999f, 0.5f);
		}
		else if (cshape <= 75.0f) {
			vshape = boundlinmap(cshape, 25.0f, 75.0f, 0.5f, 0.25f);
		}
		else {
			vshape = boundlinmap(cshape, 75.0f, 100.0f, 0.25f, 0.1f);
		}

		compress.set(vshape, vdrive, force);

		cdrive.clearnew();
		cshape.clearnew();
	}

	//gain & color --------------------------------------
	if (cgain.isnew() || ccolor.isnew() || force) {
		float vgain = 0.0f;
		if (cgain <= 25.0f) {
			vgain = boundlinmap(cgain, 0.0f, 25.0f, 0.025f, 0.5f); //true 0 gain unsupported
		}
		else {
			vgain = boundlogmap(cgain, 25.0f, 100.0f, 0.5f, 100.0f);
		}

		float vcolor = 0.0f;
		if (ccolor <= 25.0f) {
			vcolor = boundlinmap(ccolor, 0.0f, 25.0f, 2000.0f, 4000.0f);
			vcolor = expf(-2.0f * constants.pi * vcolor / Fs); //to alpha
		}
		else if (ccolor <= 50.0f) {
			vcolor = boundlinmap(ccolor, 25.0f, 50.0f, 4000.0f, 6000.0f);
			vcolor = expf(-2.0f * constants.pi * vcolor / Fs); //to alpha
		}
		else if (ccolor <= 75.0f) {
			vcolor = boundlinmap(ccolor, 50.0f, 75.0f, 6000.0f, 10000.0f);
			vcolor = expf(-2.0f * constants.pi * vcolor / Fs); //to alpha
		}
		else {
			vcolor = boundlinmap(ccolor, 75.0f, 100.0f, expf(-2.0f * constants.pi * 10000.0f / Fs), 0.0f);
		}

		distort.set(vgain, vcolor, force);

		cgain.clearnew();
		ccolor.clearnew();
	}

	//tone --------------------------------------
	if (ctone.isnew() || force) {
		ctone.clearnew();

		float vtone = 0.0f;
		if (ctone <= 50.0f) {
			vtone = boundlogmap(ctone, 0.0f, 50.0f, 0.25f, 1.0f);
		}
		else {
			vtone = boundlogmap(ctone, 50.0f, 100.0f, 1.0f, 4.0f);
		}
		slewTone.target(vtone, force);
	}
	if (slewTone.check() || force) { //sub-slew update
		slewTone.slew(samples);
		if (slewTone.get() <= 1.0f) {
			preTone.setCoefs(filt::lowshelf1(ftone / (0.5f * Fs), slewTone), force);
		}
		else {
			preTone.setCoefs(filt::highshelf1(ftone / (0.5f * Fs), 1.0f / slewTone), force);
		}
		postTone.setCoefs(filt::inverse(preTone.getCoefs()), force);
	}

	//level --------------------------------------
	if (clevel.isnew() || force) {
		clevel.clearnew();

		float vlev = 0.0f;
		if (clevel <= 50.0f) {
			vlev = boundlinmap(clevel, 0.0f, 50.0f, 0.0f, 1.0f);
		}
		else{
			vlev = boundlogmap(clevel, 50.0f, 100.0f, 1.0f, 40.0f);
		}
		outGain.setLevel(vlev, force);
	}
}

void distortion3::clear()
{
	slewTone.converge();

	preTone.clear();
	
	envAttack.clear();
	envHold.clear();
	envRelease.clear();

	distort.clear();

	postTone.clear();
}

void distortion3::step(int len)
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
		safeInput(pin, inscr, curlen);

		//update to slewed parameters
		update(curlen);
		
		//distortion preprocessing
		preTone.stepBlock(inscr.ptr(), inscr.ptr(), curlen);

		//envelope detection
		vabs(inscr.ptr(), pout, curlen);
		envAttack.stepBlock(pout, pout, curlen);
		envHold.stepBlock(pout, pout, curlen);
		envRelease.stepBlock(pout, pout, curlen);

		//form, process, and apply gain
		compress.stepGainBlock(pout, scr.ptr(), curlen);
		vmult(inscr.ptr(), scr.ptr(), inscr.ptr(), curlen); //apply compression gain on signal
		vmult(pout, scr.ptr(), pout, curlen); //apply compression gain on envelope
		distort.stepBlock(inscr, pout, scr, curlen); //saturate to the target envelope

		//post-processing
		postTone.stepBlock(scr, pout, curlen);
		outGain.stepBlock(pout, pout, curlen);
		
		if (safeOutput(pout, curlen)) {
			clear();
		}

		//step to next sub-block
		pin += curlen;
		pout += curlen;
		len -= curlen;
	}
}
