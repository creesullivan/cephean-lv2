
//--------------------------------------------------
// Defines the module's audio processing behavior
//--------------------------------------------------

#include "respectrum.h"

respectrum::respectrum(float setFs) :
	plugin(PLUG_CONTROLS, PLUG_INPUTS, PLUG_OUTPUTS, setFs),
	rebuff(maxlen),
	npts(npts0* getNominalUSR(setFs)),
	K(getKForNpts(npts)),

	wdw(npts),
	F(npts),
	fenv(K),
	fhold(K),
	frel(K),
	fsmo(K),
	shift(K, (int)roundf(slewdur* setFs / maxlen)),

	noiseLF((int)roundf(slewdur* setFs / maxlen)),
	noiseHF((int)roundf(slewdur* setFs / maxlen)),

	gattsus((int)roundf(slewdur* setFs / maxlen)),
	gsus((int)roundf(slewdur * setFs / maxlen)),
	fsustain1(K), fsustain2(K),

	fhsmo(K),
	fir(npts - maxlen, F, maxlen - 1),

	ybuff(npts),
	xbuff(npts),
	X(K),
	Hmag(K),
	Hdeemph(K),
	Hsupp(K),
	Hang(K),

	fscratch1(K),
	fscratch2(K),
	fscratch3(K),

	outLevel((int)roundf(slewdur * setFs))
{
	//design window
	mphann(wdw.ptr(), npts, 0);
	float wdwDC = sum(wdw.ptr(), npts);
	vdiv(wdw.ptr(), wdwDC, wdw.ptr(), npts); //normalize levels

	//design fixed spectral envelope detection
	fenv.setRelease(sbdet * Fs / maxlen);
	fhold.setQuarterLength((int)roundf((sbhold * Fs) / (4.0f * maxlen)));
	frel.setRelease(sbrel * Fs / maxlen);
	fsmo.setResolution(Noct, fsmomin / (Fs / 2.0f), fsmomax / (Fs / 2.0f));
	fsmo.setUnityWeighting();
	fhsmo.setSmooth(2.0f);

	//design filters
	preemph.setCoefs(filt::lowshelf1(fpre / (Fs / 2.0f), powf(10.0f, gpre / 20.0f)), true);
	deemph.setCoefs(filt::inverse(preemph.getCoefs()), true);
	
	//compute deemphasis spectra
	impz(deemph.getCoefs(), xbuff.get(npts), npts);
	F.rmag(xbuff.get(npts), Hdeemph.ptr()); //deemph spectra

	//compute linear phase angle spectra
	vset(Hmag.ptr(), 1.0f, K);
	F.linphase(Hmag.ptr(), Hang.ptr(), (npts * 3) / 8);
}
respectrum::~respectrum(){}

void respectrum::init()
{
	rebuff.reset();

	//clear state
	ybuff.clear();
	xbuff.clear();
	fhold.clear();
	frel.clear();
	fsustain1.clear();
	fsustain2.clear();
	fhsmo.clear();

	commit();
	update();
}

void respectrum::deactivate() {}

void respectrum::update(int samples)
{
	bool force = (samples == 0); //if no processing samples, we are initializing

	//gate low/high --------------------------------------
	float vgatelow = boundlinmap(gateLow, 0.0f, 100.0f, -140.0f, 0.0f);
	noiseLF.target(vgatelow, force);

	float vgatehigh = boundlinmap(gateHigh, 0.0f, 100.0f, -140.0f, 0.0f);
	noiseHF.target(vgatehigh, force);

	bool updateGate = noiseLF.slew();
	updateGate = (noiseHF.slew() || updateGate);
	if (updateGate || force) { //recompute suppression gate power
		float kmin = fnoisemin * npts / Fs;
		float kmax = fnoisemax * npts / Fs;
		float nselow = expf(2.0f * constants.db2nats * noiseLF.get());
		float nsehigh = expf(2.0f * constants.db2nats * noiseHF.get());
		float r = logf(nselow / nsehigh) / logf(kmin / kmax);
		float g = nselow / powf(kmin, r);

		Hsupp[0] = nselow;
		for (int k = 1; k < K; ++k) {
			Hsupp[k] = g * powf((float)k, r);
		}
		vbound(Hsupp.ptr(), nselow, nsehigh, Hsupp.ptr(), K);

		float bpdummy = 0.0f;
	}

	//timbre --------------------------------------
	float newvshift = boundlogmap(timbre, 0.0f, 100.0f, 4.0f, 0.25f);
	bool updateShift = (newvshift != vshift);
	if (updateShift || force){ //don't set if we don't have to
		vshift = newvshift;
		float* shiftind = shift.set();
		vincspace(shiftind, 0.0f, vshift, K);
		vmin(shiftind, K - 1.0f, shiftind, K);
		if (force) {
			shift.converge();
		}
	}

	//attack --------------------------------------
	float vattack;
	float vsustain;
	if (attack.get() <= 50.0f) {
		vattack = boundlinmap(attack, 0.0f, 50.0f, -0.3f, 1.0f);
		vsustain = 1.0f;
	}
	else {
		vattack = 1.0f;
		vsustain = boundlinmap(attack, 50.0f, 100.0f, 1.0f, -0.3f);
	}
	gattsus.target(vattack, force);
	gattsus.slew();
	gsus.target(vsustain - vattack, force);
	gsus.slew();

	//time --------------------------------------
	float vtime = boundlogmap(time, 0.0f, 100.0f, 0.03f, 1.0f);
	fsustain1.setAttack(0.5f * vtime * Fs / maxlen);
	fsustain2.setAttack(0.5f * vtime * Fs / maxlen);

	//level --------------------------------------
	float vlev;
	if (level.get() <= 50.0f) {
		vlev = boundlinmap(level, 0.0f, 50.0f, 0.0f, 1.0f);
	}
	else {
		vlev = boundlogmap(level, 50.0f, 100.0f, 1.0f, 10.0f);
	} 
	outLevel.setLevel(vlev, force);
}

void respectrum::step(int len)
{
	disableDenormals d; //scoped denormal disabler

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
		
		//prefilter
		preemph.stepBlock(pin, pout, curlen);

		//buffering and FFT analysis
		xbuff.put(pout, curlen);
		F.rmag(xbuff.get(npts), X.ptr(), wdw.ptr());
		vmult(X.ptr(), Hdeemph.ptr(), X.ptr(), K);
		X[K - 1] = 0.0f; //clear NQ for computational simplicity

		//run envelope detection
		fenv.step(X.ptr(), X.ptr(), K); //envelope for timbre
		fhold.step(X.ptr(), fscratch1.ptr(), K);
		frel.step(fscratch1.ptr(), fscratch1.ptr(), K); //envelope for noise/trans

		//noise suppression
		vmult(fscratch1.ptr(), fscratch1.ptr(), fscratch3.ptr(), K); //square sig + nse
		vsub(fscratch3.ptr(), Hsupp.ptr(), fscratch2.ptr(), K); //sig = sig + nse - nse
		vmax(fscratch2.ptr(), 0.0f, fscratch2.ptr(), K);
		vadd(fscratch3.ptr(), constants.eps, fscratch3.ptr(), K);
		vdiv(fscratch2.ptr(), fscratch3.ptr(), Hmag.ptr(), K); //W = sig/(sig + nse)

		//transient manipulation
		vmult(fscratch1.ptr(), Hmag.ptr(), fscratch1.ptr(), K); //preapply suppression
		fsustain1.step(fscratch1.ptr(), fscratch2.ptr(), K); //form sustain only
		fsustain2.step(fscratch2.ptr(), fscratch2.ptr(), K);
		vmult(fscratch1.ptr(), gattsus.get(), fscratch3.ptr(), K); //gain modifiers
		vmult(fscratch2.ptr(), gsus.get(), fscratch2.ptr(), K);
		vadd(fscratch3.ptr(), fscratch2.ptr(), fscratch3.ptr(), K);
		vmax(fscratch3.ptr(), 0.0f, fscratch3.ptr(), K); //target
		vadd(fscratch1.ptr(), constants.eps, fscratch1.ptr(), K); //orig
		vdiv(fscratch3.ptr(), fscratch1.ptr(), fscratch3.ptr(), K); //filter
		vmult(fscratch3.ptr(), Hmag.ptr(), Hmag.ptr(), K); //embed manip filter

		//preapply suppression & manipulation filter before timbre shifting
		vbound(Hmag.ptr(), mincut, maxboost, fscratch1.ptr(), K);
		hannsmooth(fscratch1.ptr(), fscratch2.ptr(), K);
		vmult(X.ptr(), fscratch2.ptr(), X.ptr(), K);
		
		//timbre shifting
		shift.apply(X.ptr(), fscratch1.ptr(), K, K);
		fsmo.maxsmooth(X.ptr(), fscratch2.ptr(), K); 
		fsmo.maxsmooth(fscratch1.ptr(), fscratch3.ptr(), K);
		vadd(fscratch2.ptr(), constants.eps, fscratch2.ptr(), K);
		vdiv(fscratch3.ptr(), fscratch2.ptr(), fscratch1.ptr(), K); //timbre weights
		vmult(fscratch1.ptr(), Hmag.ptr(), Hmag.ptr(), K); //include in net

		//smooth the net filter
		vbound(Hmag.ptr(), mincut, maxboost, Hmag.ptr(), K);
		hannsmooth(Hmag.ptr(), fscratch1.ptr(), K);
		fhsmo.step(fscratch1.ptr(), Hmag.ptr(), K);

		//cast magnitude weights to complex spectra
		complex<float>* H = fir.set();
		for (int k = 0; k < K; ++k) {
			H[k] = Hmag[k] * Hang[k];
		}

		//apply frequency domain FIR
		ybuff.put(pout, curlen);
		fir.apply(ybuff.get(npts), pout, npts, curlen);

		//postfilter
		deemph.stepBlock(pout, pout, curlen);

		//post-gain
		outLevel.stepBlock(pout, pout, curlen);

		//step to next sub-block
		pin += curlen;
		pout += curlen;
		len -= curlen;
	}
}
