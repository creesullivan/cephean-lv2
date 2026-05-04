
//--------------------------------------------------
// Defines the module's audio processing behavior
//--------------------------------------------------

#include "linkwitz.h"

linkwitz::linkwitz(float setFs) :
	plugin(PLUG_CONTROLS, PLUG_INPUTS, PLUG_OUTPUTS, setFs),
	rebuff(maxlen),
	inscr1(maxlen),
	inscr2(maxlen),

	slewFreq((int)roundf(slewdur * Fs)),
	lpf(maxlen),
	hpf(maxlen)
{}
linkwitz::~linkwitz() {}

void linkwitz::init()
{
	rebuff.reset();

	clear();
	commit();
	update();
}

void linkwitz::deactivate() {}

void linkwitz::update(int samples)
{
	bool force = (samples == 0); //if no processing samples, we are initializing

	//frequency --------------------------------------
	if (frequency.isnew() || force) {
		float vfrequency = 0.0f;
		if (frequency.get() < 25.0f) {
			vfrequency = boundlinmap(frequency, 0.0f, 25.0f, 20.0f, 200.0f);
		}
		else if (frequency.get() < 75.0f) {
			vfrequency = boundlogmap(frequency, 25.0f, 75.0f, 200.0f, 5000.0f);
		}
		else {
			vfrequency = boundlogmap(frequency, 75.0f, 100.0f, 5000.0f, 20000.0f);
		}
		slewFreq.target(vfrequency / (0.5f * Fs), force);
		frequency.clearnew();
	}
	if (slewFreq.slew(samples) || force) {
		sof::coefs clpf, chpf;
		filt::crossover2(clpf, chpf, slewFreq, 20.0f / (0.5f * Fs), 20000.0f / (0.5f * Fs));
		lpf.setCoefs(clpf, force);
		hpf.setCoefs(chpf, force);
	}
}

void linkwitz::clear()
{
	lpf.clear();
	hpf.clear();
}

void linkwitz::step(int len)
{
	//grab I/O pointers
	const float* pinLow = inLow;
	const float* pinHigh = inHigh;
	float* pout = out;

	//commit control values at the top of the block
	commit();

	//run processing loop
	int curlen = 0;
	while (len > 0) {
		curlen = rebuff.next(len); //current sub-block length

		//update to slewed parameters
		update(curlen);
		
		safeInput(pinLow, inscr1, curlen);
		safeInput(pinHigh, inscr2, curlen);

		//apply the filter pair and sum to output
		lpf.stepBlock(inscr1, inscr1, curlen);
		hpf.stepBlock(inscr2, inscr2, curlen);
		vadd(inscr1.ptr(), inscr2.ptr(), pout, curlen);
		
		if (safeOutput(pout, curlen)) {
			clear();
		}

		//step to next sub-block
		pinLow += curlen;
		pinHigh += curlen;
		pout += curlen;
		len -= curlen;
	}
}
