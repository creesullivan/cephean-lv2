
//--------------------------------------------------
// Defines the module's audio processing behavior
//--------------------------------------------------

#include "phaser.h"

phaser::phaser(float setFs) :
	plugin(PLUG_CONTROLS, PLUG_INPUTS, PLUG_OUTPUTS, setFs),
	rebuff(maxlen),
	slewedGain((int)roundf(slewdur * setFs))
{}
phaser::~phaser() {}

void phaser::init()
{
	rebuff.reset();

	commit();
	update();
}

void phaser::deactivate() {}

void phaser::update(int samples)
{
	bool force = (samples == 0); //if no processing samples, we are initializing

	//level --------------------------------------
	float vlev = boundlinmap(level, 0.0f, 100.0f, 0.0f, 2.0f);
	slewedGain.setLevel(vlev, force);
}

void phaser::step(int len)
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

		//apply the gain algorithm
		slewedGain.stepBlock(pin, pout, curlen);

		//step to next sub-block
		pin += curlen;
		pout += curlen;
		len -= curlen;
	}
}
