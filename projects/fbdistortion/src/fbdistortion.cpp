
//--------------------------------------------------
// Defines the module's audio processing behavior
//--------------------------------------------------

#include "fbdistortion.h"

fbdistortion::fbdistortion(float setFs) :
	plugin(PLUG_CONTROLS, PLUG_INPUTS, PLUG_OUTPUTS, setFs),
	rebuff(maxlen),
	inscr(maxlen),
	scscr(maxlen),

	Rlen((bufflen - chunklen)/DSR + 1),
	decay0((tstring * Fs) / (sqrtf((float)(maxper * minper)) * permult * getNominalUSR(Fs))),
	sbeta0(1.0f - expf(-1.0f / decay0)),
	
	slewColor((int)roundf(Fs * slewdur)),

	detect(bufflen * getNominalUSR(Fs), chunklen * getNominalUSR(Fs), DSR * getNominalUSR(Fs), maxlen),
	
	stateString((int)roundf(Fs * stringdur2), stringParam()), //rate changes to stringdur1 during sustain
	buffString((int)ceilf(maxper * permult * getNominalUSR(Fs)) + 1),
	fbGain((int)roundf(Fs* slewdur)),

	envAttack(),
	envHold(0.0f),
	envRelease(),

	driveGain((int)roundf(slewdur * Fs)), //TODO figure out slew/fade and optimize these once it plays well
	compress((int)roundf(slewdur* Fs)),
	distort((int)roundf(slewdur* Fs)),
	makeupGain((int)roundf(slewdur* Fs)),
	postFilter((int)roundf(slewdur* Fs)),

	filtFeedback(),
	outGain((int)roundf(slewdur* Fs)),

	R(Rlen),
	stscr(maxlen),
	scr1(maxlen),
	scr2(maxlen),
	scr3(maxlen)
{
	//set global correlation parameters
	detect.setThresh(thresh);
	detect.setMinPeriod(minper * getNominalUSR(Fs));

	//string model filters
	filtString.setCoefs(filt::lowpass1(fstring / (Fs / 2.0f),
		20.0f / (Fs / 2.0f), 20000.0f / (Fs / 2.0f)), true);
	DCreject.setCoefs(filt::highpass1(fDC / (Fs / 2.0f)), true);

	//envelope detection
	envAttack.setAttack(attdur * Fs);
	envAttack.setRelease(attdur * Fs);
	envHold.setQuarterLength((int)roundf(holddur * Fs / 4.0f));
	envRelease.setAttack(0.0f);
	envRelease.setRelease(reldur * Fs);

	filtFeedback.setCoefs(filt::lowpass2(2000.0f / (Fs / 2.0f), 0.707f), true);

	//TESTING fastratio M-compositions--------------------------------
	/*
	const float rat = 0.5f; //target output ratio
	const float drive = 10.0f; //target drive gain
	const float x0 = 0.01f;	//x input value used to derive the gain rolloff
	const int M = 8; //stages
	const int N = 80;
	const float subrat = powf(rat, 1.0f / M);
	const float subdrive = powf(1.0f / x0, (rat - 1.0f) / (M * (subrat - 1.0f))) * x0;
	const float makeup = powf(drive, -rat);

	fvec X(N);
	fvec Y1(N);
	fvec Y2(N);
	fvec T1(N);
	fvec T2(N);
	fvec x(N);
	fvec y1(N);
	fvec y2(N);

	vincspace(X.ptr(), -1.0f * N, N / (N - 1.0f), N); //dB
	for (int n = 0; n < N; ++n) {
		x[n] = powf(10.0f, X[n] / 20.0f); //to linear
	}
	vmult(x.ptr(), drive, x.ptr(), N); //apply drive

	fastratio R1;
	R1.set(rat, 1.0f, true);
	R1.stepBlock(x, y1, N);
	vmult(x.ptr(), y1.ptr(), y1.ptr(), N);
	vmult(y1.ptr(), makeup, y1.ptr(), N); //apply makeup gain

	fastratio R2;
	R2.set(subrat, subdrive, true);
	for (int m = 0; m < M; ++m) {
		R2.stepBlock(x, y2, N); //M applications
		vmult(x.ptr(), y2.ptr(), x.ptr(), N);
	}
	vcopy(x.ptr(), y2.ptr(), N);
	vmult(y2.ptr(), makeup, y2.ptr(), N); //apply makeup gain

	for (int n = 0; n < N; ++n) {
		Y1[n] = 20.0f * log10f(y1[n]); //to dB
		Y2[n] = 20.0f * log10f(y2[n]);
	}

	vmult(X.ptr(), rat, T1.ptr(), N); //target curves
	vmult(X.ptr(), subrat, T2.ptr(), N);

	fvec E1(N);
	fvec E2(N);
	vsub(T1.ptr(), Y1.ptr(), E1.ptr(), N);
	vsub(T1.ptr(), Y2.ptr(), E2.ptr(), N);

	float bpdummy = 0.0f; //plot

	*/

	//TESTING fastratio 2-compositions--------------------------------
	/*
	const float rat = 0.5f; //target output ratio
	const float drive = 1.0f; //target drive gain
	const float x0 = 0.01f;	//x input value used to derive the gain rolloff
	const float amount = 0.2f; //ratio on 0 to 1 where 1 is all distortion, and 0 is all compression
	const int N = 80;

	const float rat1 = powf(rat, 1.0f - amount);
	const float rat2 = powf(rat, amount);
	const float subdrive = powf(x0, (1.0f - rat1) * (1.0f - rat2) / (2.0f - rat1 - rat2));
	const float makeup1 = powf(drive, -rat1);
	const float makeup2 = powf(drive, -rat2);
	const float makeup = powf(drive, -rat);

	fvec X(N);
	fvec Y(N);
	fvec Y1(N);
	fvec Y2(N);
	fvec Z(N);

	fvec T(N);
	fvec T1(N);
	fvec T2(N);

	fvec x(N);
	fvec y(N);
	fvec y1(N);
	fvec y2(N);
	fvec z(N);

	vincspace(X.ptr(), -1.0f * N, N / (N - 1.0f), N); //dB
	for (int n = 0; n < N; ++n) {
		x[n] = powf(10.0f, X[n] / 20.0f); //to linear
	}
	vmult(x.ptr(), drive, x.ptr(), N); //apply drive

	fastratio R; //reference
	R.set(rat, 1.0f, true);
	R.stepBlock(x, y, N);
	vmult(y.ptr(), makeup, y.ptr(), N); //apply makeup gain

	fastratio R1; //compressor stage
	R1.set(rat1, subdrive, true);
	R1.stepBlock(x, y1, N);

	fastratio R2; //distortion stage
	R2.set(rat2, subdrive, true);
	R2.stepBlock(x, y2, N);
	R2.stepBlock(y1, z, N); //nested application

	vmult(y1.ptr(), makeup1, y1.ptr(), N); //apply makeup gain
	vmult(y2.ptr(), makeup2, y2.ptr(), N);
	vmult(z.ptr(), makeup, z.ptr(), N);

	for (int n = 0; n < N; ++n) {
		Y[n] = 20.0f * log10f(y[n]);
		Y1[n] = 20.0f * log10f(y1[n]);
		Y2[n] = 20.0f * log10f(y2[n]);
		Z[n] = 20.0f * log10f(z[n]);
	}

	vmult(X.ptr(), rat, T.ptr(), N); //target curves
	vmult(X.ptr(), rat1, T1.ptr(), N);
	vmult(X.ptr(), rat2, T2.ptr(), N);

	fvec EY(N);
	fvec EZ(N);
	vsub(T.ptr(), Y.ptr(), EY.ptr(), N);
	vsub(T.ptr(), Z.ptr(), EZ.ptr(), N);

	float bpdummy = 0.0f; //plot
	*/
}
fbdistortion::~fbdistortion() {}

fbdistortion::stringParam::stringParam(float setSustain, float setPeriod, float setConfidence) :
	sustain(setSustain), period(setPeriod), confidence(setConfidence)
{}
fbdistortion::stringParam::~stringParam() {}
bool fbdistortion::stringParam::operator ==(stringParam other) const {
	return (sustain == other.sustain) && (period == other.period) && (confidence == other.confidence);
}
bool fbdistortion::stringParam::operator !=(stringParam other) const {
	return !(*this == other);
}

void fbdistortion::init()
{
	rebuff.reset();

	clear();
	commit();
	update();
}

void fbdistortion::deactivate() {}

void fbdistortion::update(int samples)
{
	bool force = (samples == 0); //if no processing samples, we are initializing

	//drive, shape, dirty, and color  --------------------------------------
	if (drive.isnew() || shape.isnew() || dirty.isnew() || color.isnew() || force) {
		drive.clearnew();
		shape.clearnew();
		dirty.clearnew();
		color.clearnew();

		float vdrive = 0.0f;
		if (drive <= 25.0f) {
			vdrive = boundlinmap(drive, 0.0f, 25.0f, 0.01f, 0.1f); //true 0 drive unsupported
		}
		else{
			vdrive = boundlogmap(drive, 25.0f, 100.0f, 0.1f, 10.0f);
		}

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

		float vdirty = 0.0f; 
		if (dirty <= 25.0f) {
			vdirty = boundlinmap(dirty, 0.0f, 25.0f, 0.0f, 0.1f);
		}
		else {
			vdirty = boundlogmap(dirty, 25.0f, 100.0f, 0.1f, 1.0f);
		}

		float vcfreq = 0.0f;
		float vcmug = 0.0f;
		if (color <= 25.0f) {
			vcfreq = boundlinmap(color, 0.0f, 25.0f, 3000.0f, 4000.0f);
			vcmug = boundlogmap(color, 0.0f, 25.0f, 1.0f, 4.0f);
		}
		else if (color <= 50.0f) {
			vcfreq = boundlinmap(color, 25.0f, 50.0f, 4000.0f, 6000.0f);
			vcmug = boundlinmap(color, 25.0f, 50.0f, 4.0f, 4.0f);
		}
		else if (color <= 75.0f) {
			vcfreq = boundlinmap(color, 50.0f, 75.0f, 6000.0f, 9000.0f);
			vcmug = boundlinmap(color, 50.0f, 75.0f, 4.0f, 4.0f);
		}
		else {
			vcfreq = boundlogmap(color, 75.0f, 100.0f, 8000.0f, 16000.0f);
			vcmug = boundlogmap(color, 75.0f, 100.0f, 4.0f, 16.0f);
		}

		const float x0 = 0.01f;	//x input value used to derive the gain rolloff
		const float rat1 = powf(vshape, 1.0f - vdirty);
		const float rat2 = powf(vshape, vdirty);
		const float subdrive = powf(x0, (1.0f - rat1) * (1.0f - rat2) / (2.0f - rat1 - rat2));
		const float makeup = powf(vdrive, -vshape);

		driveGain.setLevel(vdrive, force);
		compress.set(rat1, subdrive, force);
		distort.set(rat2, subdrive, vcfreq / (0.5f * Fs), force);
		postFilter.setCoefs(filt::peaking2(vcfreq / (0.5f * Fs), 0.707f, vcmug), force);
		makeupGain.setLevel(makeup, force);
	}

	//feedback --------------------------------------
	if (feedback.isnew() || force) {
		feedback.clearnew();

		float vfb = 0.0f;
		if (feedback <= 25.0f) {
			vfb = boundlinmap(feedback, 0.0f, 25.0f, 0.0f, 0.1f);
		}
		else {
			vfb = boundlogmap(feedback, 25.0f, 100.0f, 0.1f, 1.0f);
		}
		fbGain.setLevel(vfb, force);
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

void fbdistortion::clear()
{
	slewColor.converge();
	
	detect.clear();

	stateString.target(stringParam(0.0f, maxper * permult * getNominalUSR(Fs), 0.0f), true);
	buffString.clear();
	filtString.clear();
	DCreject.clear();

	envAttack.clear();
	envHold.clear();
	envRelease.clear();

	distort.clear();

	postFilter.clear();

	filtFeedback.clear();
}

void fbdistortion::step(int len)
{
	//grab I/O pointers
	const float* pin = in;
	const float* psc = sc;
	float* pout = out;

	//commit control values at the top of the block
	commit();

	//run processing loop
	int curlen = 0;
	while (len > 0) {
		curlen = rebuff.next(len); //current sub-block length
		safeInput(pin, inscr, curlen);
		safeInput(psc, scscr, curlen);

		//update to slewed parameters
		update(curlen);

		//detect fundamental period and sustain margin on sidechain input
		detect.put(scscr, curlen);
		detect.corr(R);
		float newper = 0.0f;
		float Rper = detect.firstpeak(R, newper, Rmax);
		Rper = boundlinmap(Rper, Rmin, Rmax, 0.0f, 1.0f);
		if (Rper >= stateString.next().confidence) { //update period and confidence
			stateString.target(stringParam(min(Rper, stateString.next().sustain),	//save worst case sustain
				permult * newper * DSR,												//save most confident period
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
		buffString.getLinear(scr3, curlen, stateString.current().period - curlen);
		vmult(scr3.ptr(), scr1.ptr(), stscr.ptr(), curlen); //voice #1
		buffString.getLinear(scr3, curlen, stateString.last().period - curlen);
		vmultaccum(scr3.ptr(), scr2.ptr(), stscr.ptr(), curlen); //accumulate voice #2

		//apply string reflection processing and save for later
		filtString.stepBlock(stscr, scr3, curlen); //string reflection LPF

		//apply string model feedback processing on the way into distortion
		DCreject.stepBlock(stscr, stscr, curlen); //low end reject
		fbGain.stepBlock(stscr, stscr, curlen);	//feedback gain
		
		//distortion preprocessing
		vadd(inscr.ptr(), stscr.ptr(), inscr.ptr(), curlen); //include feedback from string model on main input
		driveGain.stepBlock(inscr.ptr(), inscr.ptr(), curlen); //drive gain

		//envelope detection
		vabs(inscr.ptr(), pout, curlen);
		envAttack.stepBlock(pout, pout, curlen);
		envHold.stepBlock(pout, pout, curlen);
		envRelease.stepBlock(pout, pout, curlen);

		//form, process, and apply gain
		compress.stepGainBlock(pout, pout, curlen);
		vmult(inscr.ptr(), pout, pout, curlen); //apply compression gain
		distort.stepBlock(pout, scr2, curlen);
		makeupGain.stepBlock(scr2, pout, curlen);
		postFilter.stepBlock(pout, pout, curlen);

		//inject the distortion feedback signal and step the string buffer
		filtFeedback.stepBlock(pout, scr1, curlen); //filter the distorted feedback input
		vmultaccum(scr1.ptr(), sbeta0, scr3.ptr(), curlen); //include normalization gain
		buffString.put(scr3, curlen);

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
