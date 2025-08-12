
//--------------------------------------------------
// Defines the module's audio processing behavior
//--------------------------------------------------

#include "octaver.h"

octaver::octaver(float setFs) :
	plugin(PLUG_CONTROLS, PLUG_INPUTS, PLUG_OUTPUTS, setFs),
	rebuff(maxlen),

	buff(bufflen * getNominalUSR(Fs), chunklen * getNominalUSR(Fs), DSR * getNominalUSR(Fs), maxlen),
	downoct(0.5f, blendlen * getNominalUSR(Fs), buff),
	upfifth(1.5f, blendlen * getNominalUSR(Fs), buff),
	upoct(2.0f, (blendlen / 2) * getNominalUSR(Fs), buff),
	up2oct(4.0f, (blendlen / 4) * getNominalUSR(Fs), buff),

	R(Rlen),

	gainDownOne((int)roundf(slewdur * Fs)),
	gainDry((int)roundf(slewdur* Fs)),
	gainUpHalf((int)roundf(slewdur* Fs)),
	gainUpOne((int)roundf(slewdur* Fs)),
	gainUpTwo((int)roundf(slewdur* Fs)),

	slewcolor((int)round(slewdur * Fs)),
	filtColor(maxlen),

	inscr(maxlen),
	vscr(maxlen)
{
	//set global correlation parameters
	buff.setThresh(thresh);
	buff.setMinPeriod(minper * getNominalUSR(Fs));

	//set shifter parameters
	downoct.setBlendThresh(Tblend);
	downoct.setSnapThresh(Tsnap);
	upfifth.setBlendThresh(Tblend);
	upfifth.setSnapThresh(Tsnap);
	upfifth.setMinPeriod(minper * getNominalUSR(Fs));
	upoct.setBlendThresh(Tblend);
	upoct.setSnapThresh(Tsnap);
	upoct.setMinPeriod(minper * getNominalUSR(Fs));
	up2oct.setBlendThresh(Tblend);
	up2oct.setSnapThresh(Tsnap);
	up2oct.setMinPeriod(minper * getNominalUSR(Fs));

	/* TESTING LMA FILTER 
	int niter = 128;
	int L = niter * maxlen;
	fvec h(L);
	vset(h.ptr(), 0.0f, h.size());
	h[0] = 1.0f;
	lma filt(1, 48000.0, 4096, maxlen);
	filt.setLength(4096, true);
	float* hptr = h.ptr();
	for (int it = 0; it < niter; ++it) {
		filt.stepBlock(hptr, hptr, maxlen);
		hptr += maxlen;
	}
	float hsum = sum(h.ptr(), h.size());
	float bpdummy = 0.0f;
	*/

	/* TESTING CORRELATION BUFFER
	int niter = 16;
	int L = niter * maxlen; //=4096, 2x the buff len to fill up the whole thing w/o edge effects

	fvec x(L); //input signal
	float ph = 0.0f; //sine wave testing
	float dph = 2.0f * constants.pi / 497.311f;
	//vset(x.ptr(), 0.0f, x.size());
	//x[L - maxlen] = 1.0f; //IR testing to find the beginning of R

	float* xptr = x.ptr();
	for (int it = 0; it < niter; ++it) {
		for (int n = 0; n < maxlen; ++n) {
			xptr[n] = sinf(ph); //generate sine wave
			ph += dph;
		}
		buff.put(xptr, maxlen);
		xptr += maxlen;
	}
	
	fvec R(buff.Rsize());
	buff.corr(R); //check the autocorrelation sequence

	float bpdummy = 0.0f;
	*/
}
octaver::~octaver() {}

void octaver::clear()
{
	buff.clear();
	downoct.clear();
	upfifth.clear();
	upoct.clear();
	up2oct.clear();
	filtColor.clear();
}

void octaver::init()
{
	rebuff.reset();

	clear();
	commit();
	update();
}

void octaver::deactivate() {}

void octaver::update(int samples)
{
	bool force = (samples == 0); //if no processing samples, we are initializing

	//down 1-octave level --------------------------------------
	if (downone.isnew() || force) {
		float vgain = boundlinmap(downone, 0.0f, 100.0f, 0.0f, 2.0f);
		gainDownOne.setLevel(vgain, force);
		downone.clearnew();
	}

	//dry level --------------------------------------
	if (dry.isnew() || force) {
		float vgain = boundlinmap(dry, 0.0f, 100.0f, 0.0f, 2.0f);
		gainDry.setLevel(vgain, force);
		dry.clearnew();
	}

	//up fifth level --------------------------------------
	if (uphalf.isnew() || force) {
		float vgain = boundlinmap(uphalf, 0.0f, 100.0f, 0.0f, 2.0f);
		gainUpHalf.setLevel(vgain, force);
		uphalf.clearnew();
	}

	//up 1-octave level --------------------------------------
	if (upone.isnew() || force) {
		float vgain = boundlinmap(upone, 0.0f, 100.0f, 0.0f, 2.0f);
		gainUpOne.setLevel(vgain, force);
		upone.clearnew();
	}

	//up 2-octave level --------------------------------------
	if (uptwo.isnew() || force) {
		float vgain = boundlinmap(uptwo, 0.0f, 100.0f, 0.0f, 2.0f);
		gainUpTwo.setLevel(vgain, force);
		uptwo.clearnew();
	}

	//color --------------------------------------
	if (color.isnew() || force) {
		float vcolor = 0.0f;
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
		slewcolor.target(vcolor / (0.5f * Fs), force);
		color.clearnew();
	}
	if (slewcolor.slew(samples) || force) {
		sofcasc<2>::coefs colortemp;
		filtColor.setCoefs(filt::lowpass4(colortemp, slewcolor.get(), 0.707f, 0.0f,
			20.0f / (0.5f * Fs), 20000.0f / (0.5f * Fs)), force);
	}
}

void octaver::step(int len)
{
	//grab I/O pointers
	const float* pin = in;
	float* pout = out;

	//commit control values at the top of the block
	commit();

	//run processing loop
	while (len > 0) {
		int curlen = rebuff.next(len); //current sub-block length
		safeInput(pin, inscr, curlen); //block NaN and large input values

		//update to slewed parameters
		update(curlen);

		//manage autocorrelation sequence
		buff.put(inscr, curlen);
		buff.corr(R);
		int Rdel = 0;
		float Rpeak = buff.corrpeak(R, Rdel);

		//preserving input for dry mix at the end
		vset(pout, 0.0f, curlen);

		//down octave
		downoct.stepBlock(inscr, vscr, curlen, R, Rpeak, Rdel);
		gainDownOne.stepBlock(vscr, vscr, curlen);
		vadd(pout, vscr.ptr(), pout, curlen);

		//up fifth
		upfifth.stepBlock(inscr, vscr, curlen, R, Rpeak, Rdel);
		gainUpHalf.stepBlock(vscr, vscr, curlen);
		vadd(pout, vscr.ptr(), pout, curlen);

		//up octave
		upoct.stepBlock(inscr, vscr, curlen, R, Rpeak, Rdel);
		gainUpOne.stepBlock(vscr, vscr, curlen);
		vadd(pout, vscr.ptr(), pout, curlen);

		//up 2x octave
		up2oct.stepBlock(inscr, vscr, curlen, R, Rpeak, Rdel);
		gainUpTwo.stepBlock(vscr, vscr, curlen);
		vadd(pout, vscr.ptr(), pout, curlen);

		//color filter applied to total of all shifts (not dry)
		filtColor.stepBlock(pout, pout, curlen);

		//add dry mix
		gainDry.stepBlock(pin, vscr, curlen);
		vadd(pout, vscr.ptr(), pout, curlen);

		//step to next sub-block
		pin += curlen;
		pout += curlen;
		len -= curlen;
	}

	//output safety and auto-reset for glitch recovery
	bool badState = safeOutput(pout, len);
	if (badState) {
		clear();
	}
}
