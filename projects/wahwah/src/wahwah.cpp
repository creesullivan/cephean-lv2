
//--------------------------------------------------
// Defines the module's audio processing behavior
//--------------------------------------------------

#include "wahwah.h"

wahwah::wahwah(float setFs) :
	plugin(PLUG_CONTROLS, PLUG_INPUTS, PLUG_OUTPUTS, setFs),
	rebuff(maxlen),
	inscr(maxlen),
	scr1(maxlen),

	Nbypass((int)roundf(Fs * bypdur)),

	slewFreq((int)roundf(slewdur * Fs)),
	slewRes((int)roundf(slewdur* Fs)),
	slewRolloff((int)roundf(slewdur* Fs)),

	resonator(maxlen),
	lpf(maxlen),
	hpf(maxlen),

	outGain((int)roundf(slewdur * Fs)),

	bypgain(bypatt * Fs, byprel * Fs)
{}
wahwah::~wahwah() {}

void wahwah::init()
{
	rebuff.reset();

	commit();
	update();
	clear();
}

void wahwah::deactivate() {}

void wahwah::update(int samples)
{
	bool force = (samples == 0); //if no processing samples, we are initializing

	//freq1, freq2, res1, res2 --------------------------------------
	bool filterShapesChanged = false;
	if (cfreq1.isnew() || cfreq2.isnew() || cres1.isnew() || cres2.isnew() || force) {
		vfreq1 = logf(boundlinmap(cfreq1, 0.0f, 100.0f, 150.0f, 550.0f) / (Fs / 2.0f));
		vfreq2 = logf(boundlinmap(cfreq2, 0.0f, 100.0f, 1000.0f, 3400.0f) / (Fs / 2.0f));
		vres1 = logf(boundlogmap(cres1, 0.0f, 100.0f, 1.0f, 16.0f));
		vres2 = logf(boundlogmap(cres2, 0.0f, 100.0f, 1.0f, 16.0f));

		filterShapesChanged = true;

		cfreq1.clearnew();
		cfreq2.clearnew();
		cres1.clearnew();
		cres2.clearnew();
	}

	//control -------------------------------------------
	if (ccontrol.isnew() || filterShapesChanged || force) {
		float vcontrol = ccontrol.get() / 100.0f;
		slewFreq.target((1.0f - vcontrol) * vfreq1 + vcontrol * vfreq2, force);
		slewRes.target((1.0f - vcontrol) * vres1 + vcontrol * vres2, force);

		ccontrol.clearnew();
	}
	if (slewFreq.check() || slewRes.check() || force) {
		slewFreq.slew(samples);
		slewRes.slew(samples);
		float freq = expf(slewFreq.get());
		float res = expf(0.5f * slewRes.get()) * constants.root2;
		resonator.setCoefs(filt::regain(filt::resonator2(freq, res), res / 2.0f), force);
	}
	if (nbyp < Nbypass) { //not bypassed
		if (ccontrol.get() <= bypthresh) { //count up time until we hit the bypass duration
			nbyp += samples;
		}
		else { //clear bypass counter
			nbyp = 0;
		}
		if (nbyp >= Nbypass) { //change to bypass
			byptarg = 0.0f;
		}
	}
	else { //currently bypassed
		if (ccontrol.get() >= bypthresh) { //change to not bypass
			nbyp = 0;
			byptarg = 1.0f;
		}
	}

	//rolloff --------------------------------------
	if (crolloff.isnew() || force) {
		float vrolloff = 0.0f;
		if (crolloff.get() < 25.0f) {
			vrolloff = boundlinmap(crolloff, 0.0f, 25.0f, 20.0f, 100.0f);
		}
		else {
			vrolloff = boundlinmap(crolloff, 25.0f, 100.0f, 100.0f, 400.0f);
		}
		slewRolloff.target(vrolloff / (0.5f * Fs), force);

		crolloff.clearnew();
	}
	if (slewRolloff.slew(samples) || force) {
		sof::coefs clpf, chpf;
		filt::crossover2(clpf, chpf, slewRolloff, 20.0f / (0.5f * Fs), 20000.0f / (0.5f * Fs));
		lpf.setCoefs(clpf, force);
		hpf.setCoefs(chpf, force);
	}

	//level --------------------------------------
	if(clevel.isnew() || force){
		float vlev = boundlinmap(clevel, 0.0f, 100.0f, 0.0f, 2.0f);
		outGain.setLevel(vlev, force);
		
		clevel.clearnew();
	}
}

void wahwah::clear()
{
	resonator.clear();
	lpf.clear();
	hpf.clear();

	if (ccontrol.get() <= bypthresh) { //init to bypass
		nbyp = Nbypass;
		byptarg = 0.0f;
	}
	else { //init to not bypassed
		nbyp = 0;
		byptarg = 1.0f;
	}
	bypgain.clear(byptarg);
}

void wahwah::step(int len)
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

		lpf.stepBlock(inscr, pout, curlen); //preserve lowpassed dry
		resonator.stepBlock(inscr, scr1, curlen); //apply wah filter
		hpf.stepBlock(scr1, scr1, curlen); //highpass wet signal
		vadd(scr1.ptr(), pout, pout, curlen); //mix dry/wet for rolloff control
		outGain.stepBlock(pout, pout, curlen); //output gain

		vset(scr1.ptr(), byptarg, curlen); //manage smoothed autobypass
		bypgain.stepBlock(scr1, scr1, curlen);
		vmult(pout, scr1.ptr(), pout, curlen);
		vsub(1.0f, scr1.ptr(), scr1.ptr(), curlen);
		vmultaccum(inscr.ptr(), scr1.ptr(), pout, curlen);
		
		if (safeOutput(pout, curlen)) {
			clear();
		}

		//step to next sub-block
		pin += curlen;
		pout += curlen;
		len -= curlen;
	}
}
