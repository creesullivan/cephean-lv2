
//--------------------------------------------------
// Defines the module's audio processing behavior
//--------------------------------------------------

#include "feedback.h"

feedback::feedback(float setFs) :
	plugin(PLUG_CONTROLS, PLUG_INPUTS, PLUG_OUTPUTS, setFs),
	rebuff(maxlen),
	inscr(maxlen),
	scscr(maxlen),
	fbscr(maxlen),

	Rlen((bufflen - chunklen)/DSR + 1),
	decay0((tstring * Fs) / (sqrtf((float)(maxper * minper)) * getNominalUSR(Fs))),
	sbeta0(1.0f - expf(-1.0f / decay0)),

	detect(bufflen * getNominalUSR(Fs), chunklen * getNominalUSR(Fs), DSR * getNominalUSR(Fs), maxlen),
	
	stateString((int)roundf(Fs * stringdur2), stringParam()), //rate changes to stringdur1 during sustain
	buffString(maxper * getNominalUSR(Fs) + 1),
	stmem(maxlen),
	fbGain((int)roundf(Fs* slewdur)),

	dryGain((int)roundf(Fs * slewdur)),

	R(Rlen),
	stscr(maxlen),
	scr1(maxlen),
	scr2(maxlen),
	scr3(maxlen)
{
	//set fixed correlation parameters
	detect.setThresh(thresh);
	detect.setMinPeriod(minper * getNominalUSR(Fs));

	//string model filters
	filtString.setCoefs(filt::lowpass1(fstring / (Fs / 2.0f),
		20.0f / (Fs / 2.0f), 20000.0f / (Fs / 2.0f)), true);
	DCreject.setCoefs(filt::highpass1(fDC / (Fs / 2.0f)), true);

	//feedback filter
	filtFeedback.setCoefs(filt::lowpass2(2000.0f / (Fs / 2.0f), 0.707f), true);
}
feedback::~feedback() {}

feedback::stringParam::stringParam(float setSustain, float setPeriod, float setConfidence) :
	sustain(setSustain), period(setPeriod), confidence(setConfidence)
{}
feedback::stringParam::~stringParam() {}
bool feedback::stringParam::operator ==(stringParam other) const {
	return (sustain == other.sustain) && (period == other.period) && (confidence == other.confidence);
}
bool feedback::stringParam::operator !=(stringParam other) const {
	return !(*this == other);
}

void feedback::init()
{
	rebuff.reset();

	clear();
	commit();
	update();
}

void feedback::deactivate() {}

void feedback::update(int samples)
{
	bool force = (samples == 0); //if no processing samples, we are initializing

	//amount --------------------------------------
	if (amount.isnew() || force) {
		amount.clearnew();

		float vamount = 0.0f;
		if (amount <= 10.0f) {
			vamount = boundlinmap(amount, 0.0f, 10.0f, 0.0f, 0.001f);
		}
		else {
			vamount = boundlogmap(amount, 10.0f, 100.0f, 0.001f, 1.0f);
		}
		fbGain.setLevel(vamount, force);
	}

	//dry --------------------------------------
	if (dry.isnew() || force) {
		dry.clearnew();

		float vdry = boundlinmap(dry, 0.0f, 100.0f, 0.0f, 2.0f);
		dryGain.setLevel(vdry, force);
	}
}

void feedback::clear()
{
	detect.clear();

	stateString.target(stringParam(0.0f, (float)(maxper * getNominalUSR(Fs)), 0.0f), true);
	buffString.clear();
	stmem.clear();
	filtString.clear();
	DCreject.clear();

	filtFeedback.clear();
}

void feedback::step(int len)
{
	//grab I/O pointers
	const float* pin = in;
	const float* psc = sc;
	const float* pfb = fb;
	float* pout = out;

	//commit control values at the top of the block
	commit();

	//run processing loop
	int curlen = 0;
	while (len > 0) {
		curlen = rebuff.next(len); //current sub-block length
		safeInput(pin, inscr, curlen);
		safeInput(psc, scscr, curlen);
		safeInput(pfb, fbscr, curlen);

		//update to slewed parameters
		update(curlen);

		//inject the feedback signal and step the last block's string buffer
		filtFeedback.stepBlock(fbscr, fbscr, curlen); //filter the feedback input
		vadd(stmem.get(maxlen), fbscr.ptr(), fbscr.ptr(), curlen); //add to string feedback
		buffString.put(fbscr, curlen); //step string buffer

		//detect fundamental period and sustain margin on sidechain input
		detect.put(scscr, curlen);
		detect.corr(R);
		float newper = 0.0f;
		float Rper = detect.firstpeak(R, newper, Rmax);
		Rper = boundlinmap(Rper, Rmin, Rmax, 0.0f, 1.0f);
		if (Rper >= stateString.next().confidence) { //update period and confidence
			stateString.target(stringParam(min(Rper, stateString.next().sustain),	//save worst case sustain
				newper * DSR,														//save most confident period
				Rper));																//and its confidence for checks
		}
		else {
			stateString.target(stringParam(min(Rper, stateString.next().sustain),	//save worst case sustain
				stateString.next().period,											//latch period
				stateString.next().confidence));									//and its confidence for checks
		}

		//prepare to run string model
		if (stateString.swap()) {
			stateString.target(stringParam(1.0f, stateString.current().period, 0.0f)); //init next state
		}
		float curdur = min(stateString.last().sustain, stateString.current().sustain);
		curdur = boundlogmap(curdur, 0.0f, 1.0f, stringdur2, stringdur1);
		stateString.setFade((int)roundf(Fs * curdur)); //fade across states faster when less sustained

		stateString.fadeBlock(scr1, curlen); //g
		vsub(1.0f, scr1.ptr(), scr2.ptr(), curlen); //1 - g
		vmult(scr1.ptr(), expf(-1.0f / (decay0 * stateString.current().sustain)), scr1.ptr(), curlen); //g * salpha1
		vmult(scr2.ptr(), expf(-1.0f / (decay0 * stateString.last().sustain)), scr2.ptr(), curlen); //(1 - g) * salpha2

		//get string feedback signal based on the 2-voice faded period
		buffString.getLinear(scr3, curlen, max(stateString.current().period - curlen, 0.0f));
		vmult(scr3.ptr(), scr1.ptr(), stscr.ptr(), curlen); //voice #1
		buffString.getLinear(scr3, curlen, max(stateString.last().period - curlen, 0.0f));
		vmultaccum(scr3.ptr(), scr2.ptr(), stscr.ptr(), curlen); //accumulate voice #2

		//apply string model output processing
		DCreject.stepBlock(stscr, fbscr, curlen); //low end reject
		fbGain.stepBlock(fbscr, fbscr, curlen);	//feedback gain

		//apply string reflection processing and save for next block
		filtString.stepBlock(stscr, stscr, curlen); //string reflection LPF
		stmem.put(stscr, curlen);

		//apply dry gain and output mix
		dryGain.stepBlock(inscr, inscr, curlen);
		vadd(inscr.ptr(), fbscr.ptr(), pout, curlen);
		
		if (safeOutput(pout, curlen)) {
			clear();
		}

		//step to next sub-block
		pin += curlen;
		psc += curlen;
		pfb += curlen;
		pout += curlen;
		len -= curlen;
	}
}
