
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

	freqRolloff((int)roundf(slewdur* Fs)),
	lpf(maxlen), hpf(maxlen),

	outGain((int)roundf(slewdur* Fs)),

	inscr(maxlen)
{
	//TESTING THE DEPTH I WANT TO USE -- IMPULSE RESPONSE
	int npts = 8192;
	float alpha = 0.0f;// pcore.getAlphaForDelay(96.0f);
	fvec x(npts);
	vset(x.ptr(), 0.0f, npts);
	x[0] = 1.0f;
	pcore.setDepth(0.9f, 0.45f, true); //let fb = dep/2 as an easy rule!
	pcore.clear(alpha);
	pcore.stepBlock(x, alpha, x, npts);

	int K = getKForNpts(npts);
	fvec X(K);
	fvec XdB(K);
	fft F(npts);
	F.rmag(x, X, npts);
	vmagdB(X, XdB, K);
	vmax(XdB.ptr(), -20.0f, XdB.ptr(), K);

 	float bpdummy = 0.0f;

	//if:
	//H(w) = gd + gw*e^(-1*jw)/(1 - fb*e^(-1j*w)), resp @w=pi?
	//H(pi) = gd - gw/(1 + fb) = 1 + depth <-- for positive depth only
	//H(0) = gd + gw/(1 - fb) = 1 - depth <--
	//gd = (1-depth) - gw/(1 - fb)
	//(1-depth) - gw(1/(1-fb) + 1/(1+fb)) = (1+depth),
	//-gw((1 + fb + 1 - fb = 2)/((1-fb)(1+fb)=1-fb*fb)) = 2*depth
	// gw = depth*(fb*fb-1)
}
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
		float depval = 0.0f;
		float fbval = 0.0f;
		if (depth.get() <= 75.0f) {
			depval = boundlinmap(depth, 0.0f, 75.0f, 0.0f, 1.0f);
			fbval = 0.5f * depval;
		}
		else {
			depval = 1.0f;
			fbval = boundlinmap(depth, 75.0f, 100.0f, 0.5f, 0.9f);
		}
		pcore.setDepth(depval, fbval, force);

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

	//rolloff --------------------------------------
	if (rolloff.isnew() || force) {
		float vrolloff = 0.0f;
		if (rolloff.get() < 25.0f) {
			vrolloff = boundlinmap(rolloff, 0.0f, 25.0f, 20.0f, 100.0f);
		}
		else {
			vrolloff = boundlogmap(rolloff, 25.0f, 100.0f, 100.0f, 1000.0f);
		}
		freqRolloff.target(vrolloff / (0.5f * Fs), force);

		rolloff.clearnew();
	}
	if (freqRolloff.slew(samples) || force) {
		sof::coefs clpf, chpf;
		filt::crossover2(clpf, chpf, freqRolloff, 20.0f / (0.5f * Fs), 20000.0f / (0.5f * Fs));
		lpf.setCoefs(clpf, force);
		hpf.setCoefs(chpf, force);
	}

	//level ----------------------
	if (level.isnew() || force) {
		//float tval = boundlinmap(tone, 0.0f, 100.0f, 0.0f, -0.5f);
		float lval = boundlinmap(level, 0.0f, 100.0f, 0.0f, 2.0f);
		//dryGain.setLevel(lval * tval, force);
		//wetGain.setLevel(lval * (1.0f - fabsf(tval)), force); //TODO <--- embed dry/wet gain in phasercore
		outGain.setLevel(lval, force);

		level.clearnew();
	}
}

void phaser::clear()
{
	pcore.clear();
	ph = 0.0f;
	lpf.clear();
	hpf.clear();
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
		ph += (dph * curlen);
		float curalpha = boundlogmap(triwave(ph), -1.0f, 1.0f, del1, del2);
		curalpha = pcore.getAlphaForDelay(curalpha);

		//apply the warped buffer with feedback and fixed dry/wet mix
		pcore.stepBlock(inscr, curalpha, pout, curlen);

		//crossover and mixing
		lpf.stepBlock(inscr, inscr, curlen);
		hpf.stepBlock(pout, pout, curlen);
		vadd(inscr.ptr(), pout, pout, curlen);

		//output gain
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
