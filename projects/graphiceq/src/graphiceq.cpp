
//--------------------------------------------------
// Defines the module's audio processing behavior
//--------------------------------------------------

#include "graphiceq.h"

graphiceq::graphiceq(float setFs) :
	plugin(PLUG_CONTROLS, PLUG_INPUTS, PLUG_OUTPUTS, setFs),
	rebuff(maxlen),
	casc((int)roundf(slewdur * Fs))
{
	//collect band controls into a pointer
	bands[0] = &band1;
	bands[1] = &band2;
	bands[2] = &band3;
	bands[3] = &band4;
	bands[4] = &band5;
	bands[5] = &band6;
	bands[6] = &band7;
	bands[7] = &band8;
	bands[8] = &band9;

	//TESTING THE NEW FILTERS
	/*
	sof::coefs c = filt::lowshelf2(0.1f, Qs, 0.5f);

	int npts = 4096;
	int K = getKForNpts(npts);
	fvec h(npts);
	fvec H(K);
	fft F(npts);

	impz(c, h, npts);
	F.rmag(h, H);
	vmagdB(H.ptr(), H.ptr(), K);

	float bpdummy = 0.0f;
	*/

	/*
	//test the gain matrix solver
	fvec gtarg(Nband);
	vset(gtarg.ptr(), 1.0f, Nband);
	gtarg[0] = 10.0f; //dummy target
	fvec g(Nband);
	gainSolve.resolve(gtarg, g);
	fvec gresult(Nband);
	gainSolve.check(Heq, g, gresult);
	*/
	//TEST a simple 4x4 solver
	/*
	solver S(4);
	fmat EQ(4, 4);
	vset(EQ.flatptr(), 0.0f, 16);
	for (int i = 0; i < 16; ++i) {
		EQ.flatptr()[i] = randf(); //almost surely positive definite
	}
	fvec SOL(4);
	SOL[0] = 1.0f; SOL[1] = 2.0f; SOL[2] = 3.0f; SOL[3] = 4.0f;
	fvec x(4);
	S.prepare(EQ);
	S.resolve(SOL, x);
	fvec y(4);
	S.check(EQ, x, y);
	*/
	/*
	//TEST the LR bank response with the gains
	int niter = 128;
	int npts = niter * maxlen;
	int K = getKForNpts(npts);
	fvec h0(npts);
	fmat h(Nband, npts);
	fvec H0(K);
	fmat H(Nband, K);
	fft F(npts);

	vset(h0.ptr(), 0.0f, npts);
	h0[0] = 1.0f; //impulse
	int l = 0;
	for (int i = 0; i < niter; ++i) {
		bank.analyzeBlock(h0.ptr() + l, sbscr.ptr(), maxlen);
		for (int j = 0; j < Nband; ++j) {
			vcopy(sbscr[j], h[j] + l, maxlen);
		}
		bank.synthesizeBlock(sbscr.ptr(), h0.ptr() + l, maxlen);
		l += maxlen;
	}

	F.rmag(h0.ptr(), H0.ptr(), npts);
	vmagdB(H0.ptr(), H0.ptr(), K);
	for (int j = 0; j < Nband; ++j) {
		F.rmag(h[j], H[j], npts);
		vadd(H[j], 0.00001f, H[j], K);
		vmagdB(H[j], H[j], K);
	}

	float bpdummy = 0.0f;
	*/
}
graphiceq::~graphiceq() {}

void graphiceq::init()
{
	rebuff.reset();

	casc.clear();

	commit();
	update();
}

void graphiceq::deactivate() {}

void graphiceq::update(int samples)
{
	bool force = (samples == 0); //if no processing samples, we are initializing

	for (int n = 0; n < Nband; ++n) {
		if (bands[n]->isnew() || force) {

			//map parameter to linear gain value
			float val = bands[n]->get();
			if (val <= 25.0f) {
				val = boundlogmap(val, 0.0f, 25.0f, lev0, lev25);
			}
			else if (val <= 75.0f) {
				val = boundlogmap(val, 25.0f, 75.0f, lev25, lev75);
			}
			else {
				val = boundlogmap(val, 75.0f, 100.0f, lev75, lev100);
			}

			//design the new filter coefficients
			if (n == 0) { //lowshelf filter
				casc[0].setCoefs(filt::lowshelf2(fx[0] / (0.5f * Fs), Qs, val), force);
			}
			else if (n == (Nband-1)) { //highshelf filter
				casc[Nband - 1].setCoefs(filt::highshelf2(fx[Nband - 2] / (0.5f * Fs), Qs, val), force);
			}
			else { //peaking filters
				casc[n].setCoefs(filt::peaking2(f0[n] / (0.5f * Fs), Qp, val), force);
			}

			bands[n]->clearnew();
		}
	}
}

void graphiceq::step(int len)
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

		//apply filter cascade
		casc.stepBlock(pin, pout, curlen);

		//step to next sub-block
		pin += curlen;
		pout += curlen;
		len -= curlen;
	}
}
