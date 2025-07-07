
//--------------------------------------------------
// Defines the module's audio processing behavior
//--------------------------------------------------

#include "distortion.h"

distortion::distortion(float setFs) :
	plugin(PLUG_CONTROLS, PLUG_INPUTS, PLUG_OUTPUTS, setFs),
	rebuff(maxlen),

	slewdetail((int)roundf(slewdur * Fs)),
	slewcolor((int)roundf(slewdur * Fs)),

	usr(rate / getNominalUSR(Fs)), //1, 2, or 4
	Fsu(Fs * usr),
	maxlenu(usr * maxlen),

	filtDetail(2, maxlen),

	up(usr),
	aau(),

	envAttack(attdur * Fsu, 0.0f),
	envHold(0.0f),
	envRelease(0.0f, reldur * Fsu),

	curve((int)roundf(slewdur* Fsu)),
	mixComp((int)roundf(slewdur* Fsu)),
	mixDist((int)roundf(slewdur* Fsu)),
	filtColor(maxlenu),
	phaseColor(maxlenu),
	
	aad(),
	down(usr),

	outLevel((int)roundf(slewdur * Fs)),

	uscr1(maxlenu),
	uscr2(maxlenu),
	uscr3(maxlenu)
{
	//preemphasis/deemphasis filters
	fof::coefs emphtemp = filt::highshelf1(fpre / (0.5f * Fs), gpre);
	preemph.setCoefs(emphtemp, true);
	deemph.setCoefs(filt::inverse(emphtemp), true);

	//antialiasing filter
	sofcasc<4>::coefs aatemp;
	filt::antialias8(aatemp, usr);
	aau.setCoefs(aatemp, true);
	aad.setCoefs(aatemp, true);
	/*
	fvec h(65536);
	fvec H(getKForNpts(h.size()));
	fft F(h.size());
	impz<4>(aatemp, h.ptr(), h.size());
	F.rmag(h.ptr(), H.ptr());
	vadd(H.ptr(), 0.000001f, H.ptr(), H.size());
	vmagdB(H.ptr(), H.ptr(), H.size());
	*/
	//upsampled envelope detection
	envAttack.setAttack(attdur * Fsu);
	envAttack.setRelease(attdur * Fsu);
	envHold.setQuarterLength((int)roundf(holddur * Fsu / 4.0f));
	envRelease.setAttack(0.0f);
	envRelease.setRelease(reldur * Fsu);
}
distortion::~distortion() {}

void distortion::init()
{
	rebuff.reset();

	commit();
	update();
}

void distortion::deactivate() {}

void distortion::update(int samples)
{
	bool force = (samples == 0); //if no processing samples, we are initializing

	//drive & shape --------------------------------------
	if (drive.isnew() || shape.isnew() || force) {
		drive.clearnew();
		shape.clearnew();

		float vdrive = 0.0f;
		if (drive <= 25.0f) {
			vdrive = boundlinmap(drive, 0.0f, 25.0f, 0.0f, 0.01f);
		}
		else if (drive <= 75.0f) {
			vdrive = boundlogmap(drive, 25.0f, 75.0f, 0.01f, 1.0f);
		}
		else {
			vdrive = boundlogmap(drive, 75.0f, 100.0f, 1.0f, 10.0f);
		}
		float vshape = 0.0f;
		if (shape <= 25.0f) {
			vshape = boundlinmap(shape, 0.0f, 25.0f, 1.0f, 0.5f);
		}
		else if (shape <= 75.0f) {
			vshape = boundlinmap(shape, 25.0f, 75.0f, 0.5f, 0.25f);
		}
		else {
			vshape = boundlinmap(shape, 75.0f, 100.0f, 0.25f, 0.1f);
		}
		curve.set(vshape, vdrive, force);
	}

	//dirty --------------------------------------
	if (dirty.isnew() || force) {
		dirty.clearnew();

		float vdirt = boundlinmap(dirty, 0.0f, 100.0f, 0.0f, 1.0f);
		mixDist.setLevel(vdirt, force); //linear easy mix
		mixComp.setLevel(1.0f - vdirt, force);
	}

	//detail --------------------------------------
	if (detail.isnew() || force) {
		detail.clearnew();

		float vdetail = 0.0f;
		if (detail <= 25.0f) {
			vdetail = boundlinmap(detail, 0.0f, 25.0f, 1000.0f, 3000.0f);
		}
		else if (detail <= 50.0f) {
			vdetail = boundlogmap(detail, 25.0f, 50.0f, 3000.0f, 7000.0f);
		}
		else if (detail <= 75.0f) {
			vdetail = boundlogmap(detail, 50.0f, 75.0f, 7000.0f, 11000.0f);
		}
		else {
			vdetail = boundlogmap(detail, 75.0f, 100.0f, 11000.0f, 20000.0f);
		}
		
		slewdetail.target(vdetail / (0.5f * Fs), force);
	}
	if (slewdetail.slew(samples) || force) {
		sofcasc<2>::coefs detailtemp;
		filtDetail.setCoefs(filt::lowpass4(detailtemp, slewdetail.get(), 0.707f, 0.0f,
			20.0f / (0.5f * Fs), 20000.0f / (0.5f * Fs)), force);
	}
	
	//color --------------------------------------
	if (color.isnew() || force) {
		color.clearnew();

		float vcolor = 0.0f;
		if (color <= 25.0f) {
			vcolor = boundlinmap(color, 0.0f, 25.0f, 100.0f, 1000.0f);
		}
		else if (color <= 50.0f) {
			vcolor = boundlogmap(color, 25.0f, 50.0f, 1000.0f, 4000.0f);
		}
		else if (color <= 75.0f) {
			vcolor = boundlogmap(color, 50.0f, 75.0f, 4000.0f, 8000.0f);
		}
		else {
			vcolor = boundlogmap(color, 75.0f, 100.0f, 8000.0f, 20000.0f);
		}

		slewcolor.target(vcolor / (0.5f * Fsu), force);
	}
	if(slewcolor.slew(samples) || force){
		sof::coefs colortemp = filt::lowpass2(slewcolor.get(), 0.707f);
		filtColor[0].setCoefs(colortemp, force);
		filtColor[1].setCoefs(filt::mp2mp(colortemp), force);
		phaseColor.setCoefs(filt::matchphase(colortemp), force);
	}

	//level --------------------------------------
	if (level.isnew() || force) {
		level.clearnew();

		float vlev = 0.0f;
		if (level <= 20.0f) {
			vlev = boundlinmap(level, 0.0f, 20.0f, 0.0f, powf(10.0f, -40.0f / 20.0f));
		}
		else{
			vlev = boundlogmap(level, 20.0f, 100.0f, powf(10.0f, -40.0f / 20.0f), 1.0f);
		}
		outLevel.setLevel(vlev, force);
	}
}

void distortion::step(int len)
{
	//grab I/O pointers
	const float* pin = in;
	float* pout = out;

	//commit control values at the top of the block
	commit();

	//run processing loop
	int curlen = 0;
	int curlenu = 0;
	while (len > 0) {
		curlen = rebuff.next(len); //current sub-block length

		//update to slewed parameters
		update(curlen);

		//preprocessing
		preemph.stepBlock(pin, pout, curlen);
		filtDetail.stepBlock(pout, pout, curlen);
		
		//upsample
		curlenu = curlen * usr;
		up.stepBlock(pout, uscr1.ptr(), curlen, curlenu);
		aau.stepBlock(uscr1.ptr(), uscr1.ptr(), curlenu);
		
		//envelope detection
		vabs(uscr1.ptr(), uscr2.ptr(), curlenu);
		envAttack.stepBlock(uscr2.ptr(), uscr2.ptr(), curlenu);
		envHold.stepBlock(uscr2.ptr(), uscr2.ptr(), curlenu);
		envRelease.stepBlock(uscr2.ptr(), uscr2.ptr(), curlenu);

		//form, process, and apply gain
		curve.stepBlock(uscr1.ptr(), uscr2.ptr(), uscr3.ptr(), uscr2.ptr(), curlenu);
		mixComp.stepBlock(uscr2.ptr(), uscr2.ptr(), curlenu); //compression gain
		mixDist.stepBlock(uscr3.ptr(), uscr3.ptr(), curlenu); //distortion gain
		vadd(uscr2.ptr(), uscr3.ptr(), uscr2.ptr(), curlenu); //net gain
		filtColor.stepBlock(uscr2.ptr(), uscr2.ptr(), curlenu);
		phaseColor.stepBlock(uscr1.ptr(), uscr1.ptr(), curlenu); //phase match
		vmult(uscr1.ptr(), uscr2.ptr(), uscr1.ptr(), curlenu); //apply
		
		//downsample
		aad.stepBlock(uscr1.ptr(), uscr1.ptr(), curlenu);
		down.stepBlock(uscr1.ptr(), pout, curlenu);
		
		//postprocessing
		deemph.stepBlock(pout, pout, curlen);
		outLevel.stepBlock(pout, pout, curlen);

		//step to next sub-block
		pin += curlen;
		pout += curlen;
		len -= curlen;
	}
}
