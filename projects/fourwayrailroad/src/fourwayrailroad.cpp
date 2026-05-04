
//--------------------------------------------------
// Defines the module's audio processing behavior
//--------------------------------------------------

#include "fourwayrailroad.h"

fourwayrailroad::fourwayrailroad(float setFs) :
	plugin(PLUG_CONTROLS, PLUG_INPUTS, PLUG_OUTPUTS, setFs),
	rebuff(maxlen),

	gain1((int)roundf(slewdur * Fs)),
	gain2((int)roundf(slewdur * Fs)),
	gain3((int)roundf(slewdur * Fs)),
	gain4((int)roundf(slewdur * Fs))
{}
fourwayrailroad::~fourwayrailroad() {}

void fourwayrailroad::init()
{
	rebuff.reset();

	clear();
	commit();
	update();
}

void fourwayrailroad::deactivate() {}

void fourwayrailroad::update(int samples)
{
	bool force = (samples == 0); //if no processing samples, we are initializing

	//position --------------------------------------
	if(position.isnew() || force){
		int select = (int)roundf(position.get());
		switch (select)
		{
		case 1:
			gain1.setLevel(1.0f, force);
			gain2.setLevel(0.0f, force);
			gain3.setLevel(0.0f, force);
			gain4.setLevel(0.0f, force);
			break;
		case 2:
			gain1.setLevel(0.0f, force);
			gain2.setLevel(1.0f, force);
			gain3.setLevel(0.0f, force);
			gain4.setLevel(0.0f, force);
			break;
		case 3:
			gain1.setLevel(0.0f, force);
			gain2.setLevel(0.0f, force);
			gain3.setLevel(1.0f, force);
			gain4.setLevel(0.0f, force);
			break;
		case 4:
			gain1.setLevel(0.0f, force);
			gain2.setLevel(0.0f, force);
			gain3.setLevel(0.0f, force);
			gain4.setLevel(1.0f, force);
			break;
		default: //something is wrong, mute
			gain1.setLevel(0.0f, force);
			gain2.setLevel(0.0f, force);
			gain3.setLevel(0.0f, force);
			gain4.setLevel(0.0f, force);
			break;
		}
		
		position.clearnew();
	}
}

void fourwayrailroad::clear() {}

void fourwayrailroad::step(int len)
{
	//grab I/O pointers
	const float* pin = in;
	float* pout1 = out1;
	float* pout2 = out2;
	float* pout3 = out3;
	float* pout4 = out4;

	//commit control values at the top of the block
	commit();

	//run processing loop
	int curlen = 0;
	while (len > 0) {
		curlen = rebuff.next(len); //current sub-block length

		//update to slewed parameters
		update(curlen);
		
		//for this very simple process, we don't bother wasting processing time
		//checking the inputs or outputs
		//safeInput(pin, inscr, curlen);

		//apply the slewed sub-gains
		gain1.stepBlock(pin, pout1, curlen);
		gain2.stepBlock(pin, pout2, curlen);
		gain3.stepBlock(pin, pout3, curlen);
		gain4.stepBlock(pin, pout4, curlen);
		
		//if (safeOutput(pout, curlen)) {
		//	clear();
		//}

		//step to next sub-block
		pin += curlen;
		pout1 += curlen;
		pout2 += curlen;
		pout3 += curlen;
		pout4 += curlen;
		len -= curlen;
	}
}
