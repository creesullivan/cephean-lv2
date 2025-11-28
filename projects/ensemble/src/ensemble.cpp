
//--------------------------------------------------
// Defines the module's audio processing behavior
//--------------------------------------------------

#include "ensemble.h"

ensemble::ensemble(float setFs) :
	plugin(PLUG_CONTROLS, PLUG_INPUTS, PLUG_OUTPUTS, setFs),
	rebuff(maxlen),
	byp(maxlen, (int)(bypassdur * Fs)),
	inscr(maxlen),
	delscr(maxlen),
	vscr(maxlen),
	gscr(maxlen),
	gtotscr(maxlen),

	del0(roundf(sepDelay*Fs)),

	crossoverFreq((int)roundf(slewdur* Fs)),
	lpf(maxlen),
	hpf(maxlen),

	amountGain((int)roundf(slewdur * Fs)),
	outputGain((int)roundf(slewdur* Fs)),
	buff(maxlen, 13)
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

	//TESTING voice functions
	/*
	voice v;
	v.set(0.25f, 64, 3.0f, 128, 256);
	v.clear();

	int len = 4096;
	vec<uint16_t> del(len);
	fvec g(len);
	v.stepBlock(del, g, len);
	float bpdummy = 0.0f;
	*/

	//HARCODING PARAMETERS FOR NOW------------------
	/*
	Nvoice = 8;
	for (int i = 0; i < Nvoice; ++i) {
		v[i].set(expf(constants.two2nats * (0.25f / 12.0f)) - 1.0f,
			(uint16_t)roundf(0.02f * Fs),
			1.5f,
			(unsigned int)roundf(0.25f * Fs), (unsigned int)roundf(0.5f * Fs));
	}
	amountGain.target(1.0f, true);
	outputGain.setLevel(1.0f, true);
	*/
}
ensemble::~ensemble() {}

void ensemble::init()
{
	rebuff.reset();
	
	//clear();
	commit();
	update();
	clear(); //IN THIS CASE we clear after initial update to distribute voice
			 //starting points evenly given the initial parameters
}

void ensemble::deactivate() {}

void ensemble::update(int samples)
{
	bool force = (samples == 0); //if no processing samples, we are initializing

	//mode ---------------------------------------
	if (cmode.isnew() || force) {
		switch (cmode.get())
		{
		case 0: //bypass
			byp.set(true, force);
			break;
		case 1: //doubler
			byp.set(false, force);
			Nvoice = 1;
			break;
		case 2: //2 voices
			byp.set(false, force);
			Nvoice = 2;
			break;
		case 3: //4 voices
			byp.set(false, force);
			Nvoice = 4;
			break;
		case 4: //8 voices
			byp.set(false, force);
			Nvoice = 8;
			break;
		}
		cmode.clearnew();
	}

	//voice properties --------------------------------------
	if (ctime.isnew() || cpitch.isnew() || crate.isnew() || force) {
		float vtime = boundlogmap(ctime, 0.0f, 100.0f, 0.01f, 0.1f);
		float vpitch = 0.0f;
		if (cpitch.get() < 50.0f) {
			vpitch = boundlinmap(cpitch, 0.0f, 50.0f, 0.1f, 0.25f);
		}
		else {
			vpitch = boundlinmap(cpitch, 50.0f, 100.0f, 0.25f, 1.0f);
		}
		float vlevel = 3.0f;
		float vrate = 0.0f;
		if (crate.get() < 50.0f) {
			vrate = boundlinmap(crate, 0.0f, 50.0f, 1.0f, 0.3f);
		}
		else {
			vrate = boundlinmap(crate, 50.0f, 100.0f, 0.3f, 0.05f);
		}

		//assign to all 8 voices in case of mode switch
		vtime *= Fs;
		vpitch = expf(constants.two2nats * (vpitch / 12.0f)) - 1.0f;
		unsigned int vper1 = (unsigned int)(vrate * Fs);
		unsigned int vper2 = vper1 * 2;
		for (int i = 0; i < 8; ++i) {
			v[i].set(vpitch, vtime, vlevel, vper1, vper2);
		}

		ctime.clearnew();
		cpitch.clearnew();
		crate.clearnew();
	}

	//amount --------------------------------------
	if (camount.isnew() || force) {
		float vamount = boundlinmap(camount, 0.0f, 100.0f, 0.0f, 1.0f);
		amountGain.target(vamount, force);
		camount.clearnew();
	}

	//rolloff --------------------------------------
	if (crolloff.isnew() || force) {
		float vrolloff = 0.0f;
		if (crolloff.get() < 25.0f) {
			vrolloff = boundlinmap(crolloff, 0.0f, 25.0f, 20.0f, 100.0f);
		}
		else {
			vrolloff = boundlogmap(crolloff, 25.0f, 100.0f, 100.0f, 1000.0f);
		}
		crossoverFreq.target(vrolloff / (0.5f * Fs), force);
		crolloff.clearnew();
	}
	if (crossoverFreq.slew(samples) || force) {
		sof::coefs clpf, chpf;
		filt::crossover2(clpf, chpf, crossoverFreq, 20.0f / (0.5f * Fs), 20000.0f / (0.5f * Fs));
		lpf.setCoefs(clpf, force);
		hpf.setCoefs(chpf, force);
	}

	//level --------------------------------------
	if (clevel.isnew() || force) {
		float vlevel = boundlinmap(clevel, 0.0f, 100.0f, 0.0f, 2.0f);
		outputGain.setLevel(vlevel, force);
		clevel.clearnew();
	}

	if (byp.starting()) { //converge on activation
		amountGain.converge();
		outputGain.converge();
	}
}

void ensemble::clear()
{
	buff.clear();
	lpf.clear();
	hpf.clear();
	for (int i = 0; i < Nvoice; ++i) {
		v[i].clear();
	}
}

ensemble::voice::voice()
{
	set(0.0f, 64, 0.0f, 128, 256); //dummy values
	clear();
}
ensemble::voice::~voice() {}

//pitchMax - samples/sample, [-pitchMax, pitchMax]
//delayMax - samples, [1, delayMax]
//levelMax - dB >= 0, [-levelMax, levelMax]
//periodMin/Max - sample range of each period
void ensemble::voice::set(float pitchMax, float delayMax, float levelMax,
	unsigned int periodMin, unsigned int periodMax)
{
	pmax = pitchMax;
	nmax = delayMax;
	Gmax = levelMax;
	Lmin = (int)periodMin;
	Lmax = (int)periodMax;
}

//note that the g trajectory should be normalized downstream to ensure
//steady power across multiple voices, ie gnorm = 1/sqrt(sum_n{g_n})
void ensemble::voice::stepBlock(float* del, float* g, int len)
{
	float temp = 0.0f;
	float pcur = 0.0f;
	for (int i = 0; i < len; ++i) {
		if (l++ >= L) {
			n += 0.5f * (plast + pnext) * L; //accumulate from last period
			update();
		}

		temp = ((float)l) / L;
		g[i] = glast * (1.0f - temp) + gnext * temp;
		del[i] = n + (plast + 0.5f * (pnext - plast) * temp) * (temp * L);
	}
	
}

//clears current state to random settings within the current range
void ensemble::voice::clear()
{
	n = randf() * (nmax - nmin) + nmin;
	pnext = (2.0f * randf() - 1.0f) * pmax; //immediately assigned to plast by update()
	gnext = expf(constants.db2nats * (2.0f * randf() - 1.0f) * Gmax);

	//absolute bound on pnext given n
	if (pnext < 0.0f) {
		pnext = max(pnext, (nmin - n) * (2.0f / Lmax));
	}
	else {
		pnext = min(pnext, (nmax - n) * (2.0f / Lmax));
	}

	update();
}

void ensemble::voice::update()
{
	l = 0;
	L = (int)roundf(randf() * (Lmax - Lmin)) + Lmin;
	plast = pnext;
	if (n > nmax) { //catch up to changes in delay boundaries
		pnext = -pmax;
	}
	else if (n < nmin) {
		pnext = pmax;
	}
	else {
		pnext = (2.0f * randf() - 1.0f) * pmax;
	}
	glast = gnext;
	gnext = expf(constants.db2nats * (2.0f * randf() - 1.0f) * Gmax);

	//smart bounds on pnext and L given n, plast, and planned L
	if (pnext < 0.0f) {
		//nmin <= n + (L/2)*plast + ((L+Lmax)/2)*pnext)
		float pbound = (nmin - n - (L / 2.0f) * plast) * (2.0f / (L + Lmax));
		pnext = max(pnext, min(pbound, 0.0f));
	}
	else if (pnext > 0.0f) {
		//nmax >= n + (L/2)*plast + ((L+Lmax)/2)*pnext)
		float pbound = (nmax - n - (L / 2.0f) * plast) * (2.0f / (L + Lmax));
		pnext = min(pnext, max(pbound, 0.0f));
	}
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

		//fixed processing regardless of bypass status
		buff.put(inscr, curlen);
		vset(pout, 0.0f, curlen);

		if (byp.active()) {
			//apply each voice, accumulating into the wet output
			vset(gtotscr.ptr(), 0.0f, curlen);
			for (int i = 0; i < Nvoice; ++i) {
				v[i].stepBlock(delscr, gscr, curlen); //get voice properties
				vadd(delscr.ptr(), del0, delscr.ptr(), curlen); //add fixed voice delay
				buff.getLinear(vscr, curlen, delscr); //retrieve repitched voice

				vmultaccum(vscr.ptr(), gscr.ptr(), pout, curlen); //accumulate voices into wet
				vmultaccum(gscr.ptr(), gscr.ptr(), gtotscr.ptr(), curlen); //prep to normalize power
			}

			amountGain.slewBlock(gscr.ptr(), curlen); //apply slewed amount gain to wet path
			vmult(gscr.ptr(), pout, pout, curlen);

			vmult(gscr.ptr(), gscr.ptr(), gscr.ptr(), curlen); //apply amount gain power to total
			vmult(gscr.ptr(), gtotscr.ptr(), gtotscr.ptr(), curlen);
			vadd(gtotscr.ptr(), 1.0f, gtotscr.ptr(), curlen); //include dry path power
			vsqrt(gtotscr.ptr(), gtotscr.ptr(), curlen); //form output normalization gain

			vadd(inscr.ptr(), pout, pout, curlen); //mix dry and wet in amount ratio
			vdiv(pout, gtotscr.ptr(), pout, curlen); //apply normalization

			lpf.stepBlock(inscr, inscr, curlen); //apply bass crossover
			hpf.stepBlock(pout, pout, curlen);
			vadd(inscr.ptr(), pout, pout, curlen); //mix

			outputGain.stepBlock(pout, pout, curlen); //output gain
		}

		//manage raw vs processed selection & fade to output
		byp.stepBlock(inscr, pout, pout, curlen);

		if (safeOutput(pout, curlen)) {
			clear();
		}

		//step to next sub-block
		pin += curlen;
		pout += curlen;
		len -= curlen;
	}
}
