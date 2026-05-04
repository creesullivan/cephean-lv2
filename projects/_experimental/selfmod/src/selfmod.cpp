
//--------------------------------------------------
// Defines the module's audio processing behavior
//--------------------------------------------------

#include "selfmod.h"

selfmod::selfmod(float setFs) :
	plugin(PLUG_CONTROLS, PLUG_INPUTS, PLUG_OUTPUTS, setFs),
	rebuff(maxlen),
	inscr(maxlen), scr1(maxlen), scr2(maxlen), scr3(maxlen), scr4(maxlen),

	filt0(), filt90(), filtctrl(),

	slewDepth((int)roundf(slewdur * Fs)),

	slewRolloff((int)roundf(slewdur * Fs)),
	lpf(maxlen), hpf0(maxlen), hpf90(maxlen),

	wetGain((int)roundf(slewdur * Fs)),
	dryGain((int)roundf(slewdur * Fs))
{
	//design fixed Hilbert transform filters
	sofcasc<2>::coefs coef0;
	sofcasc<4>::coefs coef90;
	filt::hilbert7(coef0, coef90);
	filt0.setCoefs(coef0, true);
	filt90.setCoefs(coef90, true);

	//TESTING with fixed bandwidth filters
	const float fbw = 1000.0f;
	filtdir0.setCoefs(filt::lowpass2(fbw / (Fs / 2.0f), 0.5f), true);
	filtdir90.setCoefs(filt::lowpass2(fbw / (Fs / 2.0f), 0.5f), true);
	filtapf.setCoefs(filt::allpass1(fbw / (Fs / 2.0f)), true);

	//design fixed control signal path
	filtctrl.setCoefs(filt::resonator2(ctrlfreq / (Fs / 2.0f), 1.0f), true);
	detect.setAttack(attdur * Fs);
	detect.setRelease(reldur * Fs);
}
selfmod::~selfmod() {}

void selfmod::init()
{
	rebuff.reset();

	clear();
	commit();
	update();
}

void selfmod::deactivate() {}

void selfmod::update(int samples)
{
	bool force = (samples == 0); //if no processing samples, we are initializing

	//depth --------------------------------------
	if(depth.isnew() || force){
		float vdepth = 0.0f;
		if (depth.get() <= 25.0f) {
			vdepth = boundlinmap(depth, 0.0f, 25.0f, 0.0f, 1.0f);
		}
		else {
			vdepth = boundlogmap(depth, 25.0f, 100.0f, 1.0f, 10.0f);
		}
		slewDepth.target(vdepth * 2.0f * constants.pi, force);
		depth.clearnew();
	}

	//rolloff --------------------------------------
	if (rolloff.isnew() || force) {
		float vrolloff = 0.0f;
		if (rolloff.get() < 25.0f) {
			vrolloff = boundlinmap(rolloff, 0.0f, 25.0f, 100.0f, 250.0f);
		}
		else {
			vrolloff = boundlogmap(rolloff, 25.0f, 100.0f, 250.0f, 5000.0f);
		}
		slewRolloff.target(vrolloff / (0.5f * Fs), force);
		rolloff.clearnew();
	}
	if (slewRolloff.slew(samples) || force) {
		sof::coefs clpf, chpf;
		filt::crossover2(clpf, chpf, slewRolloff); //NOT compensated, we want to cut deeply toward 0 to reduce 90deg path artifacts
		lpf.setCoefs(clpf, force);
		hpf0.setCoefs(chpf, force);
		hpf90.setCoefs(chpf, force);
	}

	//wet level --------------------------------------
	if (wet.isnew() || force) {
		float vgain = boundlinmap(wet, 0.0f, 100.0f, 0.0f, 2.0f);
		wetGain.setLevel(vgain, force);
		wet.clearnew();
	}

	//dry level --------------------------------------
	if (dry.isnew() || force) {
		float vgain = boundlinmap(dry, 0.0f, 100.0f, 0.0f, 2.0f);
		dryGain.setLevel(vgain, force);
		dry.clearnew();
	}

}

void selfmod::clear()
{
	filt0.clear();
	filt90.clear();

	filtdir0.clear();
	filtdir90.clear();
	filtapf.clear();
	filtctrl.clear();

	detect.clear();

	lpf.clear();
	hpf0.clear();
	hpf90.clear();
}

void selfmod::step(int len)
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

		//YE BE WARNED <-- while I find an algo that actually sounds good, I'm tearing this up and ignoring phase alignment!

		//LEFT OFF HERE <-- I think this sounds like shit and isn't worth finding a way to salvage
		//	-thought it might sound cool, but FM/PM is just too noisy on real signals
		//	-I think the sound I'm looking for could be generated from a filtered bit-crusher/aliaser
		//	-or just drop this as an idea and look at other effects for now, unrelated: like a bit-crusher

		//split the input signal into 0 and 90 deg paths
		filt90.stepBlock(inscr, scr1, curlen);
		filt0.stepBlock(inscr, inscr, curlen);

		//generate modulatee signals
		filtdir90.stepBlock(scr1, scr1, curlen);
		filtdir0.stepBlock(inscr, scr4, curlen);
		//filtapf.stepBlock(inscr, inscr, curlen); //phase align dry path

		//generate the control signal
		filtctrl.stepBlock(scr4, scr2, curlen); //includes APF cascade and bandlimit for alignment
		vabs(scr2.ptr(), scr3.ptr(), curlen);
		detect.stepBlock(scr3, scr3, curlen);
		vadd(scr3.ptr(), 1.0f / maxboost, scr3.ptr(), curlen);
		vdiv(scr2.ptr(), scr3.ptr(), scr2.ptr(), curlen);
		limit.stepBlock(scr2, scr3, curlen); 
		vmult(scr2.ptr(), scr3.ptr(), scr2.ptr(), curlen); //filtered, then compressed and saturated to unity amplitude

		//convert control signal to 0 and 90 deg gain vectors
		if (slewDepth.check()) { //slewing depth
			slewDepth.slewBlock(scr3, curlen);
			vmult(scr2.ptr(), scr3.ptr(), scr2.ptr(), curlen);
		}
		else { //fixed depth multiplier
			vmult(scr2.ptr(), slewDepth.get(), scr2.ptr(), curlen);
		}
		vcos(scr2.ptr(), scr3.ptr(), curlen);
		vsin(scr2.ptr(), scr2.ptr(), curlen);

		//apply rolloff processing -- NOTE that THIS WON'T ACTUALLY WORK and I need to split into two sections around the nonlinearity
		//lpf.stepBlock(inscr, pout, curlen);
		//hpf0.stepBlock(scr4, scr4, curlen);
		//hpf90.stepBlock(scr1, scr1, curlen);

		//apply gain vectors and mix to wet and dry outputs
		vmult(scr1.ptr(), scr2.ptr(), scr1.ptr(), curlen);
		vmultaccum(scr4.ptr(), scr3.ptr(), scr1.ptr(), curlen); //wet highs
		//vadd(scr1.ptr(), pout, scr1.ptr(), curlen); //wet plus dry low end
		//vadd(inscr.ptr(), pout, pout, curlen); //phase aligned fullband dry

		//dry/wet mix
		//dryGain.stepBlock(pout, pout, curlen);
		//wetGain.stepBlock(scr1, scr1, curlen);
		//vadd(scr1.ptr(), pout, pout, curlen);
		vcopy(scr1.ptr(), pout, curlen);
		
		if (safeOutput(pout, curlen)) {
			clear();
		}

		//step to next sub-block
		pin += curlen;
		pout += curlen;
		len -= curlen;
	}
}
