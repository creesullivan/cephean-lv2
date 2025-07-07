
//--------------------------------------------------
// Defines the module's audio processing behavior
//--------------------------------------------------

#include "diffuse.h"

diffuse::diffuse(float setFs) :
	plugin(PLUG_CONTROLS, PLUG_INPUTS, PLUG_OUTPUTS, setFs),
	rebuff(maxlen),

	dryscr(maxlen),
	slewColor((int)roundf(slewdur * Fs)),

	alg((int)roundf(slewdur * setFs), maxlen, getNominalUSR(Fs)),
	lpf(maxlen),
	apf(maxlen),
	wetGain((int)roundf(slewdur * setFs)),
	dryGain((int)roundf(slewdur* setFs))
{}
diffuse::~diffuse() {}

void diffuse::init()
{
	rebuff.reset();

	commit();
	update();
}

void diffuse::deactivate() {}

void diffuse::update(int samples)
{
	bool force = (samples == 0); //if no processing samples, we are initializing

	// seed  -----------------------------------
	if (seed.isnew() || force) {
		unsigned int vseed = (unsigned int)bound(seed.get(), 0, 15);
		alg.setSeed(vseed);
		seed.clearnew();
	}

	// spread  ---------------------------------
	if (spread.isnew() || force) {
		float vspread = boundlinmap(spread, 0.0f, 100.0f, 0.0f, 0.25f);
		alg.setDecay(vspread * Fs, force);

		spread.clearnew();
	}

	//filter color -----------------------------
	if (color.isnew() || force) {
		float vcolor;
		if (color <= 25.0f) {
			vcolor = boundlogmap(color, 0.0f, 25.0f, 1000.0f, 3000.0f);
		}
		else if (color <= 50.0f) {
			vcolor = boundlogmap(color, 25.0f, 50.0f, 3000.0f, 7000.0f);
		}
		else if (color <= 75.0f) {
			vcolor = boundlogmap(color, 50.0f, 75.0f, 7000.0f, 11000.0f);
		}
		else {
			vcolor = boundlogmap(color, 75.0f, 100.0f, 11000.0f, 20000.0f);
		}

		slewColor.target(vcolor / (Fs / 2.0f), force);
		color.clearnew();
	}
	if (slewColor.slew(samples) || force) { //assign phase aligned lowpass
		sof::coefs ctemp = filt::lowpass2(slewColor.get(), 0.707f,
			20.0f / (Fs / 2.0f), 20000.0f / (Fs / 2.0f));
		lpf[0].setCoefs(ctemp, force);
		lpf[1].setCoefs(filt::mp2mp(ctemp), force);
		apf.setCoefs(filt::matchphase(ctemp), force);
	}

	//wet --------------------------------------
	if (wet.isnew() || force) {
		float vwet = boundlinmap(wet, 0.0f, 100.0f, 0.0f, 2.0f);
		wetGain.setLevel(vwet, force);

		wet.clearnew();
	}

	//dry --------------------------------------
	if (dry.isnew() || force) {
		float vdry = boundlinmap(dry, 0.0f, 100.0f, 0.0f, 2.0f);
		dryGain.setLevel(vdry, force);

		dry.clearnew();
	}
}

void diffuse::step(int len)
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

		//dry path processing
		vcopy(pin, dryscr.ptr(), curlen);
		apf.stepBlock(dryscr.ptr(), dryscr.ptr(), curlen);
		dryGain.stepBlock(dryscr.ptr(), dryscr.ptr(), curlen);

		//wet path processing
		alg.stepBlock(pin, pout, curlen);
		lpf.stepBlock(pout, pout, curlen);
		wetGain.stepBlock(pout, pout, curlen);

		//mix
		vadd(pout, dryscr.ptr(), pout, curlen);

		//step to next sub-block
		pin += curlen;
		pout += curlen;
		len -= curlen;
	}
}
