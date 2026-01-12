
//--------------------------------------------------
// Defines the module's audio processing behavior
//--------------------------------------------------

#include "distortion2.h"

distortion2::distortion2(float setFs) :
	plugin(PLUG_CONTROLS, PLUG_INPUTS, PLUG_OUTPUTS, setFs),
	rebuff(maxlen),
	inscr(maxlen),

	slewDrive((int)roundf(Fs* slewdur)),
	slewShape((int)roundf(Fs * slewdur)),
	slewDirty((int)roundf(Fs * slewdur)),
	slewColor((int)roundf(Fs* slewdur)),
	slewTone((int)roundf(Fs * slewdur)),

	preTone(maxlen),

	envAttack(),
	envHold(0.0f),
	envRelease(),

	driveGain(maxlen),
	compress(maxlen),
	distort(maxlen),
	makeupGain(maxlen),
	makeupFilter(maxlen),

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
distortion2::~distortion2() {}

void distortion2::init()
{
	rebuff.reset();

	clear();
	commit();
	update();
}

void distortion2::deactivate() {}

void distortion2::update(int samples)
{
	bool force = (samples == 0); //if no processing samples, we are initializing

	//drive --------------------------------------
	if (drive.isnew() || force) {
		float vdrive = 0.0f;
		if (drive <= 25.0f) {
			vdrive = boundlinmap(drive, 0.0f, 25.0f, 0.01f, 0.1f); //true 0 drive unsupported
		}
		else {
			vdrive = boundlogmap(drive, 25.0f, 100.0f, 0.1f, 10.0f);
		}
		slewDrive.target(vdrive, force);

		drive.clearnew();
	}

	//shape --------------------------------------
	if (shape.isnew() || force) {
		float vshape = 0.0f;
		if (shape <= 25.0f) {
			vshape = boundlinmap(shape, 0.0f, 25.0f, 0.999f, 0.5f);
		}
		else if (shape <= 75.0f) {
			vshape = boundlinmap(shape, 25.0f, 75.0f, 0.5f, 0.25f);
		}
		else {
			vshape = boundlinmap(shape, 75.0f, 100.0f, 0.25f, 0.1f);
		}
		slewShape.target(vshape, force);

		shape.clearnew();
	}

	//dirty --------------------------------------
	if (dirty.isnew() || force) {
		float vdirty = 0.0f;
		if (dirty <= 25.0f) {
			vdirty = boundlinmap(dirty, 0.0f, 25.0f, 0.0f, 0.1f);
		}
		else {
			vdirty = boundlogmap(dirty, 25.0f, 100.0f, 0.1f, 1.0f);
		}
		slewDirty.target(vdirty, force);

		dirty.clearnew();
	}
	
	//color  --------------------------------------
	if (color.isnew() || force) {
		float vcolor = 0.0f;
		if (color <= 25.0f) {
			vcolor = boundlinmap(color, 0.0f, 25.0f, 2000.0f, 4000.0f);
		}
		else if (color <= 50.0f) {
			vcolor = boundlinmap(color, 25.0f, 50.0f, 4000.0f, 6000.0f);
		}
		else if (color <= 75.0f) {
			vcolor = boundlinmap(color, 50.0f, 75.0f, 6000.0f, 8000.0f);
		}
		else {
			vcolor = boundlogmap(color, 75.0f, 100.0f, 8000.0f, 22000.0f);
		}
		slewColor.target(vcolor, force);

		color.clearnew();
	}

	//shared update from drive, shape, dirty, and color sub-slew
	if (slewDrive.check() || slewShape.check() || slewDirty.check() || slewColor.check() || force) {
		slewDrive.slew(samples);
		slewShape.slew(samples);
		slewDirty.slew(samples);
		slewColor.slew(samples);

		const float x0 = 0.01f;	//x input value used to derive the gain rolloff
		const float rat1 = powf(slewShape, 1.0f - slewDirty);
		const float rat2 = powf(slewShape, slewDirty);
		const float subdrive = powf(x0, (1.0f - rat1) * (1.0f - rat2) / (2.0f - rat1 - rat2));
		const float makeup = powf(slewDrive, -slewShape);
		const float rfreq = min(slewColor.get() / (0.5f * Fs), 0.95f);

		driveGain.setLevel(slewDrive, force);
		compress.set(rat1, subdrive, force);
		distort.set(rat2, subdrive, rfreq, force);
		makeupFilter.setCoefs(filt::peaking2(rfreq, 0.707f, 4.0f), force);
		makeupGain.setLevel(makeup, force);
	}

	//tone --------------------------------------
	if (tone.isnew() || force) {
		tone.clearnew();

		float vtone = 0.0f;
		if (tone <= 50.0f) {
			vtone = boundlogmap(tone, 0.0f, 50.0f, 0.25f, 1.0f);
		}
		else {
			vtone = boundlogmap(tone, 50.0f, 100.0f, 1.0f, 4.0f);
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
	if (level.isnew() || force) {
		level.clearnew();

		float vlev = 0.0f;
		if (level <= 20.0f) {
			vlev = boundlinmap(level, 0.0f, 20.0f, 0.0f, powf(10.0f, -40.0f / 20.0f));
		}
		else {
			vlev = boundlogmap(level, 20.0f, 100.0f, powf(10.0f, -40.0f / 20.0f), 1.0f);
		}
		outGain.setLevel(vlev, force);
	}
}

void distortion2::clear()
{
	slewDrive.converge();
	slewShape.converge();
	slewDirty.converge();
	slewColor.converge();
	slewTone.converge();

	preTone.clear();
	
	envAttack.clear();
	envHold.clear();
	envRelease.clear();

	distort.clear();

	makeupFilter.clear();
	postTone.clear();
}

void distortion2::step(int len)
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
		driveGain.stepBlock(inscr.ptr(), inscr.ptr(), curlen);
		preTone.stepBlock(inscr.ptr(), inscr.ptr(), curlen);

		//envelope detection
		vabs(inscr.ptr(), pout, curlen);
		envAttack.stepBlock(pout, pout, curlen);
		envHold.stepBlock(pout, pout, curlen);
		envRelease.stepBlock(pout, pout, curlen);

		//form, process, and apply gain
		compress.stepGainBlock(pout, pout, curlen);
		vmult(inscr.ptr(), pout, pout, curlen); //apply compression gain
		distort.stepBlock(pout, scr, curlen);
		makeupGain.stepBlock(scr, pout, curlen);
		makeupFilter.stepBlock(pout, pout, curlen);

		//post-processing
		postTone.stepBlock(pout, pout, curlen);
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
