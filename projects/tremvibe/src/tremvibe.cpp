
//--------------------------------------------------
// Defines the module's audio processing behavior
//--------------------------------------------------

#include "tremvibe.h"

tremvibe::tremvibe(float setFs) :
	plugin(PLUG_CONTROLS, PLUG_INPUTS, PLUG_OUTPUTS, setFs),
	rebuff(maxlen),

	maxdelN((int)roundf(maxdel * setFs)),
	vbuff(maxdelN + maxlen)
{
	//TESTING linwave and quadwave
	fvec x(1000);
	vincspace(x.ptr(), -2.0f * constants.pi, 3.0f * 2.0f * constants.pi * 0.001f, x.size());
	fvec ylin(x.size());
	fvec yquad(x.size());

	for (int i = 0; i < x.size(); ++i) {
		ylin[i] = triwave(x[i]);
		yquad[i] = quadwave(x[i]);
	}

	float bpdummy = 0.0f;
}
tremvibe::~tremvibe() {}

void tremvibe::init()
{
	rebuff.reset();

	vbuff.clear();
	ph = 0.0f;

	commit();
	update();
}

void tremvibe::deactivate() {}

void tremvibe::update(int samples)
{
	bool force = (samples == 0); //if no processing samples, we are initializing

	//tremelo --------------------------------------
	if (tremelo < 25.0f) {
		tdepth = boundlinmap(tremelo, 0.0f, 25.0f, -1.0f, -0.25f);
	}
	else if (tremelo < 50.0f) {
		tdepth = boundlinmap(tremelo, 25.0f, 50.0f, -0.25f, 0.0f);
	}
	else if (tremelo < 75.0f) {
		tdepth = boundlinmap(tremelo, 50.0f, 75.0f, 0.0f, 0.25f);
	}
	else {
		tdepth = boundlinmap(tremelo, 75.0f, 100.0f, 0.25f, 1.0f);
	}

	//vibrato --------------------------------------
	if (vibrato < 50.0f) {
		vdepth = boundlogmap(vibrato, 0.0f, 50.0f, powf(2.0f, 0.0f / 1200.0f), powf(2.0f, 40.0f / 1200.0f)) - 1.0f;
	}
	else if (vibrato < 75.0f) {
		vdepth = boundlogmap(vibrato, 50.0f, 75.0f, powf(2.0f, 40.0f / 1200.0f), powf(2.0f, 80.0f / 1200.0f)) - 1.0f;
	}
	else {
		vdepth = boundlogmap(vibrato, 75.0f, 100.0f, powf(2.0f, 80.0f / 1200.0f), powf(2.0f, 200.0f / 1200.0f)) - 1.0f;
	}

	//rate --------------------------------------
	if (rate < 50.0f) {
		dph = boundlogmap(rate, 0.0f, 50.0f, 0.7f, 4.5f) * 2.0f * constants.pi / Fs;
	}
	else {
		dph = boundlogmap(rate, 50.0f, 100.0f, 4.5f, 20.0f) * 2.0f * constants.pi / Fs;
	}

	//random --------------------------------------
	rmult = boundlinmap(random, 0.0f, 100.0f, 0.0f, 0.5f);

	if (force) {
		newperiod();
	}
}

void tremvibe::newperiod()
{
	curdph = dph * (1.0f + rmult * (randf() - 0.5f));
	gT = bound(tdepth * (1.0f + rmult * (randf() - 0.5f)), -1.0f, 1.0f);
	gV = bound((0.5f * constants.pi * vdepth / curdph) * (1.0f + rmult * (randf() - 0.5f)), 0.0f, (float)maxdelN);
}

void tremvibe::step(int len)
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

		//put new data to shift the delay buffer
		vbuff.put(pin, curlen);

		float temp1, temp2;
		const float* pbuff = vbuff.get(maxdelN + 1, curlen - 1);
		for (int i = 0; i < curlen; ++i) { //sample process
			temp1 = ph;
			ph += curdph;
			if (ph.get() < temp1) { //indicates overrun
				newperiod();
			}

			//apply vibrato
			//temp1 = gV * (0.5f - 0.5f * cosf(ph)); //current delay
			temp1 = gV * (0.5f - 0.5f * quadwave(ph.get() + 0.5f * constants.pi));
			temp1 = interp(pbuff++, maxdelN - temp1, maxdelN + 1);

			//apply tremelo
			//temp2 = 1.0f + gT * sinf(ph); //current gain
			temp2 = 1.0f + gT * triwave(ph);
			pout[i] = temp1 * temp2;
		}

		//step to next sub-block
		pin += curlen;
		pout += curlen;
		len -= curlen;
	}
}
