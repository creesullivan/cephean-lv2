
//--------------------------------------------------
// Defines the module's audio processing behavior
//--------------------------------------------------

#include "reguitar.h"

reguitar::reguitar(float setFs) :
	plugin(PLUG_CONTROLS, PLUG_INPUTS, PLUG_OUTPUTS, setFs),
	rebuff(maxlen),
	inscr(maxlen),

	filt((int)roundf(slewdur * setFs))
{}
reguitar::~reguitar() {}

void reguitar::init()
{
	rebuff.reset();

	clear();
	commit();
	update();
}

void reguitar::deactivate() {}

void reguitar::update(int samples)
{
	bool force = (samples == 0); //if no processing samples, we are initializing

	//all parameters used to update the self-slewed filter --------------------------------------
	if(freq1.isnew() || res1.isnew() || freq2.isnew() || res2.isnew() || gain.isnew() || force){
		
		float f1 = freq1.get() / (Fs / 2.0f);
		float r1 = expf(constants.db2nats * res1.get());
		float f2 = freq2.get() / (Fs / 2.0f);
		float r2 = expf(constants.db2nats * res2.get());
		float g = expf(constants.db2nats * gain.get());

		//design input/output filters
		sof::coefs coef1 = filt::lowpass2(f1, r1);
		sof::coefs coef2 = filt::lowpass2(f2, r2);
		
		//generate and apply matching filter (a wonky shelf)
		g *= (float)((coef2.b0 + coef2.b1 + coef2.b2) / (coef1.b0 + coef1.b1 + coef1.b2)); //net gain including DC norm
		coef2.b0 = g;
		coef2.b1 = -g * coef1.na1;
		coef2.b2 = -g * coef1.na2;
		filt.setCoefs(coef2, force);

		freq1.clearnew();
		res1.clearnew();
		freq2.clearnew();
		res2.clearnew();
		gain.clearnew();
	}
}

void reguitar::clear() {}

void reguitar::step(int len)
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

		//apply the algorithm as a single second order filter
		filt.stepBlock(pin, pout, curlen);
		
		if (safeOutput(pout, curlen)) {
			clear();
		}

		//step to next sub-block
		pin += curlen;
		pout += curlen;
		len -= curlen;
	}
}
