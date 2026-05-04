
//--------------------------------------------------
// Defines the module's audio processing behavior
//--------------------------------------------------

#include "doubleswitch.h"

doubleswitch::doubleswitch(float setFs) :
	plugin(PLUG_CONTROLS, PLUG_INPUTS, PLUG_OUTPUTS, setFs),
	rebuff(maxlen),
	inscr1(maxlen), inscr2(maxlen), inscr3(maxlen), inscr4(maxlen),

	gain1((int)roundf(slewdur * Fs)),
	gain2((int)roundf(slewdur * Fs)),
	gain3((int)roundf(slewdur * Fs)),
	gain4((int)roundf(slewdur * Fs))
{}
doubleswitch::~doubleswitch() {}

void doubleswitch::init()
{
	rebuff.reset();

	clear();
	commit();
	update();
}

void doubleswitch::deactivate() {}

void doubleswitch::update(int samples)
{
	bool force = (samples == 0); //if no processing samples, we are initializing

	//switch1 & switch2 --------------------------------------
	if(switch1.isnew() || switch2.isnew() || force){
		int select = (int)roundf(switch1.get()) + (int)roundf(2.0f * switch2.get()) + 1; //1, 2, 3, 4
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
		
		switch1.clearnew();
		switch2.clearnew();
	}
}

void doubleswitch::clear() {}

void doubleswitch::step(int len)
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
