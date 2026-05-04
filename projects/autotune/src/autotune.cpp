
//--------------------------------------------------
// Defines the module's audio processing behavior
//--------------------------------------------------

#include "autotune.h"

autotune::autotune(float setFs) :
	plugin(PLUG_CONTROLS, PLUG_INPUTS, PLUG_OUTPUTS, setFs),
	rebuff(maxlen),
	inscr(maxlen),

	buff(bufflen * getNominalUSR(Fs), chunklen * getNominalUSR(Fs), DSR * getNominalUSR(Fs), maxlen),
	shift(1.0f / maxShift, maxShift, blendlen * getNominalUSR(Fs), false, buff),

	R(Rlen),
	scr(maxlen)
{
	//set global correlation parameters
	buff.setThresh(thresh);
	buff.setMinPeriod(minper * getNominalUSR(Fs));

	//set shifter parameters
	shift.setBlendThresh(Tblend);
	shift.setSnapThresh(0.0f); //ignore transient snapping
	shift.setMinPeriod(minper * getNominalUSR(Fs));
}
autotune::~autotune() {}

void autotune::init()
{
	rebuff.reset();

	clear();
	commit();
	update();
}

void autotune::deactivate() {}

void autotune::update(int samples)
{
	bool force = (samples == 0); //if no processing samples, we are initializing

	//fundamental, mode, period, pitch1, pitch2 --------------------------------------
	if(fundamental.isnew() || mode.isnew() || tuning.isnew() || pitch1.isnew() || pitch2.isnew() || force){
		
		curscale.set((scale::fundamental)fundamental.get(), (scale::mode)mode.get(), tuning.get());
		maxPitchPeriod = bound((int)ceilf((Fs / DSR) / curscale.toHz(pitch1.get())), minper / DSR + 1, maxper / DSR - 1);
		minPitchPeriod = bound((int)floorf((Fs / DSR) / curscale.toHz(pitch2.get())), minper / DSR + 1, maxper / DSR - 1);
		minPitchPeriod = min(minPitchPeriod, maxPitchPeriod);
		
		fundamental.clearnew();
		mode.clearnew();
		tuning.clearnew();
		pitch1.clearnew();
		pitch2.clearnew();
	}

	//smooth ------------------------------------------------------------------------
	if (smooth.isnew() || force) {
		float vsmooth = smooth.get();
		if (vsmooth <= 50.0f) {
			vsmooth = boundlogmap(vsmooth, 0.0f, 50.0f, 0.001f, 0.04f);
		}
		else {
			vsmooth = boundlogmap(vsmooth, 50.0f, 100.0f, 0.04f, 0.25f);
		}

		salpha = expf(-1.0f / (vsmooth * Fs));
		sbeta = 1.0f - salpha;
		
		smooth.clearnew();
	}
	
}

void autotune::clear()
{
	buff.clear();
	shift.clear();
	curdsamp = 1.0f;
}

void autotune::step(int len)
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

		//manage autocorrelation sequence and pitch detection
		buff.put(inscr, curlen);
		buff.corr(R);
		
		float curper = 0.0f;
		float Rpeak = buff.firstpeak(R, curper, Tblend, minPitchPeriod, maxPitchPeriod);
		curper *= DSR;

		//autotune pitch shift computation logic
		float nextdsamp = curdsamp;
		if (Rpeak >= Ttrans) { //steady state -- do pitch target logic and smoothing
			float updsamp = 0.0f;
			float downdsamp = 0.0f;
			curscale.snap(Fs / curper, downdsamp, updsamp);
			downdsamp = curper / (Fs / downdsamp); //detected period / nearest valid period of a lower frequency = sample/sample downshift < 1
			updsamp = curper / (Fs / updsamp); //detected period / nearest valid period of a higher frequency = sample/sample upshift >= 1
			if (downdsamp * downdsamp >= (1.0f / updsamp)) { //close enough to force target the downshift
				nextdsamp = downdsamp;
			}
			else if (updsamp * updsamp <= (1.0f / downdsamp)) { //close enough to force target the upshift
				nextdsamp = updsamp;
			}
			else { //in the gray area between valid notes of the scale, so latch to whatever is nearer to our current state (hysteresis)
				if (fabsf(downdsamp - curdsamp) <= fabsf(updsamp - curdsamp)) {
					nextdsamp = downdsamp;
				}
				else {
					nextdsamp = updsamp;
				}
			}
		} //otherwise in a transient where we latch the shift (nextdsamp = curdsamp)

		//smooth the pitch shift sample step
		for (int i = 0; i < curlen; ++i) {
			curdsamp = salpha * curdsamp + sbeta * nextdsamp;
			scr[i] = curdsamp;
		}
		
		/*
		//test pitch oscillator
		float curper = 0.0f;
		float Rpeak = buff.firstpeak(R, curper, Tblend, minPitchPeriod, maxPitchPeriod);
		for (int i = 0; i < curlen; ++i) {
			scr[i] = 1.1f + 0.1f * sinf(2.0f * constants.pi * ph);
			ph += (1.0f / 192000.0f); //about one period per 4 seconds
		}
		*/
		/*
		//test fixed pitch
		int Rdel = 0;
		float Rpeak = buff.corrpeak(R, Rdel);
		float Rpeak = buff.firstpeak(R, curper, Tblend);
		vset(scr.ptr(), 1.3f, len);
		*/

		//apply pitch shifting
		// so this doesn't work perfectly, it's almost there, but still some rare glitches
		// I think it's good enough for now, let's use it for a while then come back if I don't like it
		shift.stepBlock(inscr, scr, pout, len, R, Rpeak);

		//do a dry/wet mix for DEBUGGING THE LATENCY
		//vadd(inscr.ptr(), pout, pout, len);
		
		if (safeOutput(pout, curlen)) {
			clear();
		}

		//step to next sub-block
		pin += curlen;
		pout += curlen;
		len -= curlen;
	}
}
