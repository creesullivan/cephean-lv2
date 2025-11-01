
//--------------------------------------------------
// Defines the module's audio processing behavior
//--------------------------------------------------

#include "ensemble.h"

ensemble::ensemble(float setFs) :
	plugin(PLUG_CONTROLS, PLUG_INPUTS, PLUG_OUTPUTS, setFs),
	rebuff(maxlen),
	inscr(maxlen),
	delay(12)		//FOR POW2BUFFER, ~495 -> 423 TICKS/1000 SEC in whole step function
	//delay(4096)		//FOR CIRCBUFFER, ~438 TICKS/1000 SEC in whole step function
{
	/*
	//TESTING MY MAYBE VERY FAST CIRCULAR INDEXING IDEA
	fvec data(4);
	vincspace(data.ptr(), 0.0f, 1.0f, data.size());
	pow2circint ind(2, 0);

	int N = 1024;
	fvec output(N);

	for (int n = 0; n < N; ++n) {
		output[n] = data[ind++];
	}
	
	float bpdummy = 0.0f;
	//LEFT OFF HERE <-- hooray! this worked, next step is making a quick
	//time delay with a larger buffer.
	//
	// THEN I want to profile this vs a circbuff.getNN() and a flatbuff.getNN(),
	// probably gotta write a pow2buff that closely mirrors circbuff. Use 4096
	// as the representative duration.
	*/

}
ensemble::~ensemble() {}

void ensemble::init()
{
	rebuff.reset();
	
	clear();
	commit();
	update();
}

void ensemble::deactivate() {}

void ensemble::update(int samples)
{
	bool force = (samples == 0); //if no processing samples, we are initializing

	/*
	//level --------------------------------------
	if(level.isnew() || force){
		float vlev = boundlinmap(level, 0.0f, 100.0f, 0.0f, 2.0f);
		slewedGain.setLevel(vlev, force);
		
		level.clearnew();
	}
	*/
}

void ensemble::clear()
{
	delay.clear();
}

void ensemble::step(int len)
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

		//TESTING <-- run time delay test
		delay.put(inscr, curlen);
		delay.get(pout, curlen, 100);

		//copy to output
		//vcopy(inscr.ptr(), pout, curlen);
		
		if (safeOutput(pout, curlen)) {
			clear();
		}

		//step to next sub-block
		pin += curlen;
		pout += curlen;
		len -= curlen;
	}
}
