
//--------------------------------------------------
// Defines the module's audio processing behavior
//--------------------------------------------------

#include "fourwayswitch.h"

fourwayswitch::fourwayswitch(float setFs) :
	plugin(PLUG_CONTROLS, PLUG_INPUTS, PLUG_OUTPUTS, setFs),
	rebuff(maxlen),
	inscr1(maxlen), inscr2(maxlen), inscr3(maxlen), inscr4(maxlen),

	gain1((int)roundf(slewdur * Fs)),
	gain2((int)roundf(slewdur * Fs)),
	gain3((int)roundf(slewdur * Fs)),
	gain4((int)roundf(slewdur * Fs))
{}
fourwayswitch::~fourwayswitch() {}

void fourwayswitch::init()
{
	rebuff.reset();

	clear();
	commit();
	update();
}

void fourwayswitch::deactivate() {}

void fourwayswitch::update(int samples)
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

void fourwayswitch::clear() {}

void fourwayswitch::step(int len)
{
	//grab I/O pointers
	const float* pin1 = in1;
	const float* pin2 = in2;
	const float* pin3 = in3;
	const float* pin4 = in4;
	float* pout = out;

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
		gain1.stepBlock(pin1, inscr1, curlen);
		gain2.stepBlock(pin2, inscr2, curlen);
		gain3.stepBlock(pin3, inscr3, curlen);
		gain4.stepBlock(pin4, inscr4, curlen);

		//sum to output
		vadd(inscr1.ptr(), inscr2.ptr(), pout, curlen);
		vadd(pout, inscr3.ptr(), pout, curlen);
		vadd(pout, inscr4.ptr(), pout, curlen);
		
		//if (safeOutput(pout, curlen)) {
		//	clear();
		//}

		//step to next sub-block
		pin1 += curlen;
		pin2 += curlen;
		pin3 += curlen;
		pin4 += curlen;
		pout += curlen;
		len -= curlen;
	}
}
