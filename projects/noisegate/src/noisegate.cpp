
//--------------------------------------------------
// Defines the module's audio processing behavior
//--------------------------------------------------

#include "noisegate.h"

noisegate::noisegate(float setFs) :
	plugin(PLUG_CONTROLS, PLUG_INPUTS, PLUG_OUTPUTS, setFs),
	rebuff(maxlen),
	inscr(maxlen), scr1(maxlen), scr2(maxlen), scr3(maxlen),
	Nlat((int)roundf(latency * Fs)),

	lowave(1, 1.0 * Fs, (int)ceilf(rmsdur * Fs) + 1, maxlen),
	lowdelay(Nlat + maxlen),
	lowrel(0.0f, reldur * Fs),
	lowthresh((int)roundf(slew*Fs)),

	midave(1, 1.0 * Fs, (int)ceilf(rmsdur* Fs) + 1, maxlen),
	middelay(Nlat + maxlen),
	midrel(0.0f, reldur * Fs),
	midthresh((int)roundf(slew* Fs)),

	highave(1, 1.0 * Fs, (int)ceilf(rmsdur* Fs) + 1, maxlen),
	highdelay(Nlat + maxlen),
	highrel(0.0f, reldur * Fs),
	highthresh((int)roundf(slew* Fs))
{
	//design filters
	lowfilt[0].setCoefs(filt::lowpass2(freqlow / (Fs / 2.0f), 0.5f), true);
	lowfilt[1].setCoefs(filt::merge(filt::allpass1(freqhigh/ (Fs / 2.0f)), filt::unity1()), true);
	midfilt[0].setCoefs(filt::highpass2(freqlow / (Fs / 2.0f), 0.5f), true);
	midfilt[1].setCoefs(filt::invert(filt::lowpass2(freqhigh / (Fs / 2.0f), 0.5f)), true); //invert to align with low and high
	highfilt[0].setCoefs(filt::highpass2(freqlow / (Fs / 2.0f), 0.5f), true);
	highfilt[1].setCoefs(filt::highpass2(freqhigh / (Fs / 2.0f), 0.5f), true);

	//compute normalization power multipliers
	const int L = 4096 * getNominalUSR(Fs);
	fvec hlow(L);
	fvec hmid(L);
	fvec hhigh(L);
	fvec h0(L);
	vset(h0.ptr(), 0.0f, L); h0[0] = 1.0f;

	lowfilt.stepBlock(h0, hlow, L);
	lowthreshmult = ssq(hlow.ptr(), L);
	lowfilt.clear();

	midfilt.stepBlock(h0, hmid, L);
	midthreshmult = ssq(hmid.ptr(), L);
	midfilt.clear();

	highfilt.stepBlock(h0, hhigh, L);
	highthreshmult = ssq(hhigh.ptr(), L);
	highfilt.clear();
	/*
	vadd(hlow.ptr(), hmid.ptr(), h0.ptr(), L);
	vadd(h0.ptr(), hhigh.ptr(), h0.ptr(), L); //DEBUGGING to check the filters sum to flat
	const int K = getKForNpts(L);
	fvec Hnet(K);
	fft F(L);
	F.rmag(h0, Hnet);
	*/

	//design envelope detectors
	lowave.setLength((int)roundf(rmsdur * Fs), true);
	lowhold.setQuarterLength((int)roundf(holddur * Fs / 4.0f));
	midave.setLength((int)roundf(rmsdur * Fs), true);
	midhold.setQuarterLength((int)roundf(holddur * Fs / 4.0f));
	highave.setLength((int)roundf(rmsdur * Fs), true);
	highhold.setQuarterLength((int)roundf(holddur * Fs / 4.0f));
}
noisegate::~noisegate() {}

void noisegate::init()
{
	rebuff.reset();

	clear();
	commit();
	update();
}

void noisegate::deactivate() {}

void noisegate::update(int samples)
{
	bool force = (samples == 0); //if no processing samples, we are initializing

	//low & high --------------------------------------
	if(low.isnew() || high.isnew() || force){
		float vlow = 0.0f;
		if (low.get() < 10.0f) {
			vlow = boundlinmap(low, 0.0f, 10.0f, 0.0f, 0.00001f);
		}
		else {
			vlow = boundlogmap(low, 10.0f, 100.0f, 0.00001f, 1.0f);
		}

		float vhigh = 0.0f;
		if (high.get() < 10.0f) {
			vhigh = boundlinmap(high, 0.0f, 10.0f, 0.0f, 0.00001f);
		}
		else {
			vhigh = boundlogmap(high, 10.0f, 100.0f, 0.00001f, 1.0f);
		}
		
		float vmid = sqrtf(vlow * vhigh); //dB average of low and high thresholds

		lowthresh.target(vlow * vlow * lowthreshmult, force);
		midthresh.target(vmid * vmid * midthreshmult, force);
		highthresh.target(vhigh * vhigh * highthreshmult, force);
		
		low.clearnew();
		high.clearnew();
	}
}

void noisegate::clear()
{
	lowave.clear();
	lowdelay.clear();
	lowfilt.clear();
	lowhold.clear();
	lowrel.clear();

	midave.clear();
	middelay.clear();
	midfilt.clear();
	midhold.clear();
	midrel.clear();

	highave.clear();
	highdelay.clear();
	highfilt.clear();
	highhold.clear();
	highrel.clear();
}

void noisegate::step(int len)
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

		//low band -------------------------------------
		lowfilt.stepBlock(inscr, scr1, curlen); //band signal
		vmult(scr1.ptr(), scr1.ptr(), scr2.ptr(), curlen); //to power
		lowave.stepBlock(scr2, scr2, curlen);
		lowhold.stepBlock(scr2, scr2, curlen);
		lowrel.stepBlock(scr2, scr2, curlen); //envelope
		lowthresh.slewBlock(scr3, curlen); //threshold
		vsub(scr2.ptr(), scr3.ptr(), scr3.ptr(), curlen);
		vmax(scr3.ptr(), 0.0f, scr3.ptr(), curlen); //numerator
		vadd(scr2.ptr(), constants.eps, scr2.ptr(), curlen); //denominator
		vdiv(scr3.ptr(), scr2.ptr(), scr3.ptr(), curlen); //Wiener gain
		lowdelay.put(scr1, curlen); //delay alignment
		vmult(lowdelay.get(curlen, Nlat), scr3.ptr(), pout, curlen); //apply into output

		//mid band -------------------------------------
		midfilt.stepBlock(inscr, scr1, curlen); //band signal
		vmult(scr1.ptr(), scr1.ptr(), scr2.ptr(), curlen); //to power
		midave.stepBlock(scr2, scr2, curlen);
		midhold.stepBlock(scr2, scr2, curlen);
		midrel.stepBlock(scr2, scr2, curlen); //envelope
		midthresh.slewBlock(scr3, curlen); //threshold
		vsub(scr2.ptr(), scr3.ptr(), scr3.ptr(), curlen);
		vmax(scr3.ptr(), 0.0f, scr3.ptr(), curlen); //numerator
		vadd(scr2.ptr(), constants.eps, scr2.ptr(), curlen); //denominator
		vdiv(scr3.ptr(), scr2.ptr(), scr3.ptr(), curlen); //Wiener gain
		middelay.put(scr1, curlen); //delay alignment
		vmultaccum(middelay.get(curlen, Nlat), scr3.ptr(), pout, curlen); //apply into output

		//high band -------------------------------------
		highfilt.stepBlock(inscr, scr1, curlen); //band signal
		vmult(scr1.ptr(), scr1.ptr(), scr2.ptr(), curlen); //to power
		highave.stepBlock(scr2, scr2, curlen);
		highhold.stepBlock(scr2, scr2, curlen);
		highrel.stepBlock(scr2, scr2, curlen); //envelope
		highthresh.slewBlock(scr3, curlen); //threshold
		vsub(scr2.ptr(), scr3.ptr(), scr3.ptr(), curlen);
		vmax(scr3.ptr(), 0.0f, scr3.ptr(), curlen); //numerator
		vadd(scr2.ptr(), constants.eps, scr2.ptr(), curlen); //denominator
		vdiv(scr3.ptr(), scr2.ptr(), scr3.ptr(), curlen); //Wiener gain
		highdelay.put(scr1, curlen); //delay alignment
		vmultaccum(highdelay.get(curlen, Nlat), scr3.ptr(), pout, curlen); //apply into output

		if (safeOutput(pout, curlen)) {
			clear();
		}

		//step to next sub-block
		pin += curlen;
		pout += curlen;
		len -= curlen;
	}
}
