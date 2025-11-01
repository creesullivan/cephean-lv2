
#pragma once

//--------------------------------------------------
// Defines first, second, and higher order filters
//--------------------------------------------------

#include "cephean-lv2.h"

namespace cephean
{

//==================================================

//First order filter class, defining the fof::coefs structure
//and permitting mono, stereo, or multichannel filtering with
//fixed or time varying slewed filters.
class fof : public monoalg, public stereoalg, public multialg
{
	template<unsigned int N> friend class lrcascade;

public:
	struct coefs
	{
		double na1;
		double b0;
		double b1;
	};

	fof(int slew = 1, int maxNch = 1, bool blockTransposed = false);
	~fof();

	void reallocate(int slew = 1, int maxNch = 1, bool blockTransposed = false);

	void setCoefs(coefs newCoefs, bool converge = false);
	coefs getCoefs() const;

	void clear();

	float step(float x) override;
	void stepBlock(const float* x, float* y, int len) override;

	void step(float xL, float xR, float& yL, float& yR) override;
	void stepBlock(const float* xL, const float* xR, float* yL, float* yR, int len) override;

	void step(const float* x, float* y, int chan) override;
	void stepBlock(const float* const* x, float* const* y, int chan, int len) override;

private:
	slewedvec<double> sc; //coefficient slew
	coefs c; //last set target coefficients

	double memL, memR; //mono/stereo filter state
	dvec vmem; //multichannel filter state
};

namespace filt
{

// First order filter definitions
//====================================================

//Inverts the phase of the passed in coefficients structure
fof::coefs invert(const fof::coefs& coef);

//Returns the filter inverse of a stable minimum phase coefs structure
fof::coefs inverse(const fof::coefs& coef);

//Returns an allpass filter that phase matches a double application of
// //the provided filter with one stage's numerator coefficients flipped
fof::coefs matchphase(const fof::coefs& coef);

//Flips the numerator coefficients of the provided filter, switching
//min phase <-> max phase, usually for phase aligned double filtering
fof::coefs mp2mp(const fof::coefs& coef);

//====================================================

//First order unity filter, returns default
fof::coefs unity1();

//First order trivial gain filter with gain multiplier g
fof::coefs gain1(float g);

//====================================================

//First order lowpass with rolloff frequency f
fof::coefs lowpass1(float f);

//First order compensated lowpass with rolloff frequency f and
//bounds minf, maxf for stability
fof::coefs lowpass1(float f, float minf, float maxf);

//First order highpass with rolloff frequency f
fof::coefs highpass1(float f);

//First order compensated highpass with rolloff frequency f and
//bounds minf, maxf for stability
fof::coefs highpass1(float f, float minf, float maxf);

//====================================================

//First order low shelf with crossover frequency f and skirt gain multiplier g
fof::coefs lowshelf1(float f, float g);

//First order compensated low shelf with crossover f, skirt gain multiplier g,
//and bounds minf, maxf for stability
fof::coefs lowshelf1(float f, float g, float minf, float maxf);

//First order high shelf with crossover frequency f and skirt gain multiplier g
fof::coefs highshelf1(float f, float g);

//First order compensated high shelf with crossover f, skirt gain multiplier g,
//and bounds minf, maxf for stability
fof::coefs highshelf1(float f, float g, float minf, float maxf);

//====================================================

//First order Pade delay approximation for del in samples
fof::coefs pade1(float del);

//First order allpass with crossover frequency f, phase >f is 180 deg
fof::coefs allpass1(float f);

//Simple one-pole smoother applied as a first order filter
fof::coefs smooth1(float N);

}

//==================================================

//Second order filter class, defining the sof::coefs structure
//and permitting mono, stereo, or multichannel filtering with
//fixed or time varying slewed filters.
class sof : public monoalg, public stereoalg, public multialg
{
	template<unsigned int N> friend class sofcasc;
	template<unsigned int N> friend class lrcascade;

public:
	struct coefs
	{
		double na1;
		double na2;
		double b0;
		double b1;
		double b2;
	};

	sof(int slew = 1, int maxNch = 1, bool blockTransposed = false);
	~sof();

	void reallocate(int slew = 1, int maxNch = 1, bool blockTransposed = false);

	void setCoefs(coefs newCoefs, bool converge = false);
	coefs getCoefs() const;

	void clear();

	float step(float x) override;
	void stepBlock(const float* x, float* y, int len) override;

	void step(float xL, float xR, float& yL, float& yR) override;
	void stepBlock(const float* xL, const float* xR, float* yL, float* yR, int len) override;

	void step(const float* x, float* y, int chan) override;
	void stepBlock(const float* const* x, float* const* y, int chan, int len) override;

private:
	slewedvec<double> sc; //coefficient slew
	coefs c; //last set target coefficients

	double mem1L, mem1R, mem2L, mem2R; //mono/stereo filter state
	dvec vmem1, vmem2; //multichannel filter state

};

namespace filt
{

	// Second order filter definitions
	//====================================================

	//Inverts the phase of the passed in coefficients structure
	sof::coefs invert(const sof::coefs& coef);

	//Returns the filter inverse of a stable minimum phase coefs structure
	sof::coefs inverse(const sof::coefs& coef);

	//Merge two first order filters into a single second order one
	sof::coefs merge(const fof::coefs& coef1, const fof::coefs& coef2);

	//Applies bilinear warping so that content at frequency f1 appears at
	//frequency f2 in the output coefficient set
	sof::coefs warp(const sof::coefs& coef, float f1, float f2);

	//Returns an allpass filter that phase matches a double application of
	//the provided filter with one stage's numerator coefficients flipped
	sof::coefs matchphase(const sof::coefs& coef);

	//Flips the numerator coefficients of the provided filter, switching
	//min phase <-> max phase, usually for phase aligned double filtering
	sof::coefs mp2mp(const sof::coefs& coef);

	//====================================================

	//Second order unity filter, returns default
	sof::coefs unity2();

	//Second order trivial gain filter with gain multiplier g
	sof::coefs gain2(float g);

	//2-pole (second order) smoother with time constant of roughly L
	//samples and resonance of r >= 1. 
	sof::coefs smoother2(float L, float r = 1.0f);

	//==================================================

	//Second order lowpass with rolloff frequency f and resonance
	//multiplier r (r = Q)
	sof::coefs lowpass2(float f, float r);

	//Second order lowpass with rolloff frequency f, resonance
	//multiplier r (r = Q), and skirt depth multiplier d
	sof::coefs lowpass2(float f, float r, float d);

	//Second order compensated lowpass with rolloff frequency f, resonance
	//multiplier r, bounds minf, maxf for stability, and optional skirt depth
	//multiplier d
	sof::coefs lowpass2(float f, float r, float minf, float maxf, float d = 0.0f);

	//Second order highpass with rolloff frequency f and resonance
	//multiplier r (r = Q)
	sof::coefs highpass2(float f, float r);

	//Second order highpass with rolloff frequency f, resonance
	//multiplier r (r = Q), and skirt depth multiplier d
	sof::coefs highpass2(float f, float r, float d);

	//Second order compensated highpass with rolloff frequency f, resonance
	//multiplier r, bounds minf, maxf for stability, and optional skirt depth
	//multiplier d
	sof::coefs highpass2(float f, float r, float minf, float maxf, float d = 0.0f);

	//Second order lowpass/highpass pair with shared frequency and r of 0.5,
	//forming a Linkwitz-Riley crossover that sums coherently to flat.
	void crossover2(sof::coefs& lpf, sof::coefs& hpf, float f);

	//Second order compensated lowpass/highpass pair with shared frequency
	//and bounds and r of 0.5, forming a Linkwitz-Riley crossover that sums
	//coherently to flat, and that approaches true single-band operation as
	//the frequency approaches the bounds.
	void crossover2(sof::coefs& lpf, sof::coefs& hpf, float f, float minf, float maxf);

	//====================================================

	//Second order bandpass with center frequency f and octave bandwidth oct
	sof::coefs bandpass2(float f, float oct);

	//TODO -- COMPENSATED BANDPASS

	//Second order resonator with center frequency f and resonance multiplier r
	sof::coefs resonator2(float f, float r);

	//TODO -- COMPENSATED RESONATOR

	//====================================================

	//Second order low shelf filter with crossover frequency f, quality factor
	//Q, and low frequency linear gain g
	sof::coefs lowshelf2(float f, float Q, float g);

	//TODO -- COMPENSATED LOWSHELF

	//Second order high shelf filter with crossover frequency f, quality factor
	//Q, and high frequency linear gain g
	sof::coefs highshelf2(float f, float Q, float g);

	//TODO -- COMPENSATED HIGHSHELF

	//====================================================

	//Second order boost/cut peaking filter with center frequency f, quality
	//factor Q, and linear peak gain g
	sof::coefs peaking2(float f, float Q, float g);

	//TODO -- COMPENSATED... PEAKING...?

	//====================================================

	//TODO -- PHASE

}


//==================================================

//Cascade of N second order filter sections
template<unsigned int N> class sofcasc : public monoalg, public stereoalg, public multialg
{
public:
	typedef sof::coefs coefs[N];

	sofcasc(int slew = 1, int maxNch = 1, bool blockTransposed = false) :
		multialg(maxNch, blockTransposed)
	{
		for (int n = 0; n < N; ++n) {
			p[n].reallocate(slew, maxNch, blockTransposed);
			ct[n] = p[n].sc.get();
		}
	}
	~sofcasc() {}

	const sof& operator[](int ind) const
	{
		return p[ind];
	}
	sof& operator[](int ind)
	{
		return p[ind];
	}

	void reallocate(int slew = 1, int maxNch = 1, bool blockTransposed = false)
	{
		multialg::reallocate(maxNch, blockTransposed);
		for (int n = 0; n < N; ++n) {
			p[n].reallocate(slew, maxNch, blockTransposed);
			ct[n] = p[n].sc.get();
		}
	}

	void setCoefs(const coefs& newCoefs, bool converge = false)
	{
		for (int n = 0; n < N; ++n) {
			p[n].setCoefs(newCoefs[n], converge);
		}
	}
	void getCoefs(coefs& putCoefs) const
	{
		for (int n = 0; n < N; ++n) {
			putCoefs[n] = p[n].getCoefs();
		}
	}

	void clear()
	{
		for (int n = 0; n < N; ++n) {
			p[n].clear();
		}
	}

	float step(float x) override
	{
		double xtemp[N];
		double ytemp = x;
		for (int n = 0; n < N; ++n) {
			p[n].sc.slew(); //slew first
		}
		for (int n = 0; n < N; ++n) {
			xtemp[n] = ytemp + ct[n][0] * p[n].mem1L + ct[n][1] * p[n].mem2L;
			ytemp = ct[n][2] * xtemp[n] + ct[n][3] * p[n].mem1L + ct[n][4] * p[n].mem2L;
		}
		for (int n = 0; n < N; ++n) {
			p[n].mem2L = p[n].mem1L;
			p[n].mem1L = xtemp[n];
		}
		return (float)ytemp;
	}

	void stepBlock(const float* x, float* y, int len) override
	{
		bool slewing = false; //determine if we need to slew
		for (int n = 0; n < N; ++n) {
			if (p[n].sc.check()) {
				slewing = true;
			}
		}
		if (slewing) { //slew every sample
			monoalg::stepBlock(x, y, len);
		}
		else { //faster block-based application
			double xtemp[N];
			double ytemp = 0.0;
			for (int i = 0; i < len; ++i) {
				ytemp = x[i];
				for (int n = 0; n < N; ++n) {
					xtemp[n] = ytemp + ct[n][0] * p[n].mem1L + ct[n][1] * p[n].mem2L;
					ytemp = ct[n][2] * xtemp[n] + ct[n][3] * p[n].mem1L + ct[n][4] * p[n].mem2L;
				}
				for (int n = 0; n < N; ++n) {
					p[n].mem2L = p[n].mem1L;
					p[n].mem1L = xtemp[n];
				}
				y[i] = (float)ytemp;
			}
		}
	}

	void step(float xL, float xR, float& yL, float& yR) override
	{
		double xtempL[N];
		double xtempR[N];
		double ytempL = xL;
		double ytempR = xR;
		for (int n = 0; n < N; ++n) {
			p[n].sc.slew(); //slew first
		}
		for (int n = 0; n < N; ++n) {
			xtempL[n] = ytempL + ct[n][0] * p[n].mem1L + ct[n][1] * p[n].mem2L;
			xtempR[n] = ytempR + ct[n][0] * p[n].mem1R + ct[n][1] * p[n].mem2R;
			ytempL = ct[n][2] * xtempL[n] + ct[n][3] * p[n].mem1L + ct[n][4] * p[n].mem2L;
			ytempR = ct[n][2] * xtempR[n] + ct[n][3] * p[n].mem1R + ct[n][4] * p[n].mem2R;
		}
		for (int n = 0; n < N; ++n) {
			p[n].mem2L = p[n].mem1L;
			p[n].mem2R = p[n].mem1R;
			p[n].mem1L = xtempL[n];
			p[n].mem1R = xtempR[n];
		}
		yL = (float)ytempL;
		yR = (float)ytempR;
	}
	void stepBlock(const float* xL, const float* xR, float* yL, float* yR, int len) override
	{
		bool slewing = false; //determine if we need to slew
		for (int n = 0; n < N; ++n) {
			if (p[n].sc.check()) {
				slewing = true;
			}
		}
		if (slewing) { //slew every sample
			stereoalg::stepBlock(xL, xR, yL, yR, len);
		}
		else { //faster block-based application
			double xtempL[N];
			double xtempR[N];
			double ytempL = 0.0;
			double ytempR = 0.0;
			for (int i = 0; i < len; ++i) {
				ytempL = xL[i];
				ytempR = xR[i];
				for (int n = 0; n < N; ++n) {
					xtempL[n] = ytempL + ct[n][0] * p[n].mem1L + ct[n][1] * p[n].mem2L;
					xtempR[n] = ytempR + ct[n][0] * p[n].mem1R + ct[n][1] * p[n].mem2R;
					ytempL = ct[n][2] * xtempL[n] + ct[n][3] * p[n].mem1L + ct[n][4] * p[n].mem2L;
					ytempR = ct[n][2] * xtempR[n] + ct[n][3] * p[n].mem1R + ct[n][4] * p[n].mem2R;
				}
				for (int n = 0; n < N; ++n) {
					p[n].mem2L = p[n].mem1L;
					p[n].mem2R = p[n].mem1R;
					p[n].mem1L = xtempL[n];
					p[n].mem1R = xtempR[n];
				}
				yL[i] = (float)ytempL;
				yR[i] = (float)ytempR;
			}
		}
	}

	void step(const float* x, float* y, int chan) override
	{
		double xtemp[N];
		double ytemp = 0.0;
		for (int n = 0; n < N; ++n) {
			p[n].sc.slew(); //slew first
		}
		for (int m = 0; m < chan; ++m) {
			ytemp = x[m];
			for (int n = 0; n < N; ++n) {
				xtemp[n] = ytemp + ct[n][0] * p[n].vmem1[m] + ct[n][1] * p[n].vmem2[m];
				ytemp = ct[n][2] * xtemp[n] + ct[n][3] * p[n].vmem1[m] + ct[n][4] * p[n].vmem2[m];
			}
			for (int n = 0; n < N; ++n) {
				p[n].vmem2[m] = p[n].vmem1[m];
				p[n].vmem1[m] = xtemp[n];
			}
			y[m] = (float)ytemp;
		}
	}
	void stepBlock(const float* const* x, float* const* y, int chan, int len) override
	{
		bool slewing = false; //determine if we need to slew
		for (int n = 0; n < N; ++n) {
			if (p[n].sc.check()) {
				slewing = true;
			}
		}
		if (slewing) { //slew every sample
			multialg::stepBlock(x, y, chan, len);
		}
		else { //faster block-based application
			double xtemp[N];
			double ytemp = 0.0;
			if (transposed) { //len by chan
				for (int i = 0; i < len; ++i) {
					for (int m = 0; m < chan; ++m) {
						ytemp = x[i][m];
						for (int n = 0; n < N; ++n) {
							xtemp[n] = ytemp + ct[n][0] * p[n].vmem1[m] + ct[n][1] * p[n].vmem2[m];
							ytemp = ct[n][2] * xtemp[n] + ct[n][3] * p[n].vmem1[m] + ct[n][4] * p[n].vmem2[m];
						}
						for (int n = 0; n < N; ++n) {
							p[n].vmem2[m] = p[n].vmem1[m];
							p[n].vmem1[m] = xtemp[n];
						}
						y[i][m] = (float)ytemp;
					}
				}
			}
			else { //chan by len
				for (int m = 0; m < chan; ++m) {
					for (int i = 0; i < len; ++i) {
						ytemp = x[m][i];
						for (int n = 0; n < N; ++n) {
							xtemp[n] = ytemp + ct[n][0] * p[n].vmem1[m] + ct[n][1] * p[n].vmem2[m];
							ytemp = ct[n][2] * xtemp[n] + ct[n][3] * p[n].vmem1[m] + ct[n][4] * p[n].vmem2[m];
						}
						for (int n = 0; n < N; ++n) {
							p[n].vmem2[m] = p[n].vmem1[m];
							p[n].vmem1[m] = xtemp[n];
						}
						y[m][i] = (float)ytemp;
					}
				}
			}
		}
	}

private:
	sof p[N]; //list of sof stages
	const double* ct[N]; //preallocated list of pointers to slewed coefs
};

namespace filt
{

// Second order filter cascade definitions
//====================================================

//Applies bilinear warping so that content at frequency f1 appears at
//frequency f2 in the output coefficient set. Warping is applied in place
//destructively.
template<unsigned int N> const typename sofcasc<N>::coefs& warp(typename sofcasc<N>::coefs& coef, float f1, float f2)
{
	double A = tan(0.5 * constants.pi * f1) / tan(0.5 * constants.pi * f2);
	A = (1.0 - A) / (1.0 + A);
	double AA = A * A;
	double A2 = 2.0 * A;
	double AAp1 = AA + 1.0;
	double a0 = 0.0;

	sof::coefs ctemp;

	for (int n = 0; n < N; ++n) {
		ctemp = coef[n];
		a0 = 1.0 - A * ctemp.na1 - AA * ctemp.na2;

		coef[n].b0 = (ctemp.b0 + A * ctemp.b1 + AA * ctemp.b2) / a0;
		coef[n].b1 = (A2 * ctemp.b0 + AAp1 * ctemp.b1 + A2 * ctemp.b2) / a0;
		coef[n].b2 = (AA * ctemp.b0 + A * ctemp.b1 + ctemp.b2) / a0;
		coef[n].na1 = (-A2 + AAp1 * ctemp.na1 + A2 * ctemp.na2) / a0;
		coef[n].na2 = (-AA + A * ctemp.na1 + ctemp.na2) / a0;
	}
	return coef;
}

//====================================================

//Fourth order lowpass with rolloff frequency f, and resonance
//multiplier r (r = Q), and spread factor s where 0 is tightest
//and 1 is widest (s of 1 is a double filter)
const sofcasc<2>::coefs& lowpass4(sofcasc<2>::coefs& coef, float f, float r, float s = 0.0f);

//Fourth order lowpass with rolloff frequency f, resonance
//multiplier r (r = Q), spread factor s on [0, 1], and skirt
//depth multiplier d
const sofcasc<2>::coefs& lowpass4(sofcasc<2>::coefs& coef, float f, float r, float s, float d);

//Fourth order compensated lowpass with rolloff frequency f, resonance
//multiplier r, spread factor s on [0, 1], bounds minf, maxf for stability,
//and optional skirt depth multiplier d
const sofcasc<2>::coefs& lowpass4(sofcasc<2>::coefs& coef, float f, float r, float s,
	float minf, float maxf, float d = 0.0f);
	
//====================================================

//Fourth order bandpass with center frequency f and octave bandwidth oct derived
//from a frequency matched pair of lowpass/highpass filters with unity gain at
//the center frequency
const sofcasc<2>::coefs& bandpass4(sofcasc<2>::coefs& coef, float f, float oct);

//====================================================

//8th order (4 section) antialiasing lowpass filter with a bandwidth
//of 1/DSR the Nyquist rate
const sofcasc<4>::coefs& antialias8(sofcasc<4>::coefs& coef, int DSR);

}

//==================================================

//Linkwitz-Riley cascade filterbank providing perfect equalization
//via analysis/synthesis processing that splits a signal into several
//frequency bands, then recombines them again with no phase issues
template<unsigned int N> class lrcascade
{
public:
	lrcascade(int slew = 1, bool blockTransposed = false)
	{
		reallocate(slew, blockTransposed);
	}
	~lrcascade() {}

	void reallocate(int slew = 1, bool blockTransposed = false)
	{
		transposed = blockTransposed;
		for (int m = 0; m < N-1; ++m) {
			lpf[m].reallocate(slew, 1, false);
			clpf[m] = lpf[m].sc.get();

			hpf[m].reallocate(slew, 1, false);
			chpf[m] = hpf[m].sc.get();

			apf[m].reallocate(slew, 1, false);
			capf[m] = apf[m].sc.get();
		}
	}

	//Updates the normalized crossover frequencies, for an N-band bank,
	//there should be N-1 crossover frequencies
	void setCrossoverFrequencies(const float* f, bool converge = false)
	{
		for (int m = 0; m < N-1; ++m) {
			lpf[m].setCoefs(filt::lowpass2(f[m], 0.5f), converge);
			hpf[m].setCoefs(filt::invert(filt::highpass2(f[m], 0.5f)), converge);
			apf[m].setCoefs(filt::allpass1(f[m]), converge);
		}
	}

	sof::coefs getLPFCoefs(int n) const
	{
		return lpf[n].getCoefs();
	}
	sof::coefs getHPFCoefs(int n) const
	{
		return hpf[n].getCoefs();
	}
	sof::coefs getAPFCoefs(int n) const
	{
		return apf[n].getCoefs();
	}

	void clear()
	{
		for (int m = 0; m < N-1; ++m) {
			lpf[m].clear();
			hpf[m].clear(); //memory actually unused
			apf[m].clear();
		}
	}

	void analyze(float x, float* X)
	{
		double xtemp[N-1];
		double ztemp = 0.0;
		double ytemp = x;
		for (int m = 0; m < N-1; ++m) {
			lpf[m].sc.slew(); //slew analysis filters first
			hpf[m].sc.slew();
		}
		for (int m = 0; m < N-1; ++m) {
			xtemp[m] = ytemp + clpf[m][0] * lpf[m].mem1L + clpf[m][1] * lpf[m].mem2L; //shared denom & memory
			ztemp = clpf[m][2] * xtemp[m] + clpf[m][3] * lpf[m].mem1L + clpf[m][4] * lpf[m].mem2L; //LPF
			ytemp = chpf[m][2] * xtemp[m] + chpf[m][3] * lpf[m].mem1L + chpf[m][4] * lpf[m].mem2L; //HPF
			X[m] = (float)ztemp; //lowpass output
		}
		X[N-1] = (float)ytemp; //final highpass output
		for (int m = 0; m < N-1; ++m) { //only use LPF memory for efficiency
			lpf[m].mem2L = lpf[m].mem1L;
			lpf[m].mem1L = xtemp[m];
		}
	}
	void analyzeBlock(const float* x, float* const* X, int len)
	{
		bool slewing = lpf[0].sc.check(); //all filter slews are coupled
		if (slewing) { //slew every sample
			if (transposed) {
				for (int i = 0; i < len; ++i) {
					analyze(x[i], X[i]);
				}
			}
			else {
				for (int i = 0; i < len; ++i) {
					analyze(x[i], scr);
					for (int n = 0; n < N; ++n) {
						X[n][i] = scr[n];
					}
				}
			}
		}
		else { //faster block-based application
			double xtemp[N-1];
			double ztemp = 0.0;
			double ytemp = 0.0;
			if (transposed) {
				for (int i = 0; i < len; ++i) {
					ytemp = (double)x[i];
					for (int m = 0; m < N-1; ++m) {
						xtemp[m] = ytemp + clpf[m][0] * lpf[m].mem1L + clpf[m][1] * lpf[m].mem2L; //shared denom & memory
						ztemp = clpf[m][2] * xtemp[m] + clpf[m][3] * lpf[m].mem1L + clpf[m][4] * lpf[m].mem2L; //LPF
						ytemp = chpf[m][2] * xtemp[m] + chpf[m][3] * lpf[m].mem1L + chpf[m][4] * lpf[m].mem2L; //HPF
						X[i][m] = (float)ztemp; //lowpass output
					}
					X[i][N - 1] = (float)ytemp; //final highpass output
					for (int m = 0; m < N-1; ++m) { //only use LPF memory for efficiency
						lpf[m].mem2L = lpf[m].mem1L;
						lpf[m].mem1L = xtemp[m];
					}
				}
			}
			else {
				for (int i = 0; i < len; ++i) {
					ytemp = (double)x[i];
					for (int m = 0; m < N-1; ++m) {
						xtemp[m] = ytemp + clpf[m][0] * lpf[m].mem1L + clpf[m][1] * lpf[m].mem2L; //shared denom & memory
						ztemp = clpf[m][2] * xtemp[m] + clpf[m][3] * lpf[m].mem1L + clpf[m][4] * lpf[m].mem2L; //LPF
						ytemp = chpf[m][2] * xtemp[m] + chpf[m][3] * lpf[m].mem1L + chpf[m][4] * lpf[m].mem2L; //HPF
						X[m][i] = (float)ztemp; //lowpass output
					}
					X[N - 1][i] = (float)ytemp; //final highpass output
					for (int m = 0; m < N-1; ++m) { //only use LPF memory for efficiency
						lpf[m].mem2L = lpf[m].mem1L;
						lpf[m].mem1L = xtemp[m];
					}
				}
			}
		}
	}

	float synthesize(const float* Y)
	{
		double xtemp[N-1];
		double ytemp = (double)Y[0];
		for (int m = 0; m < N-1; ++m) {
			apf[m].sc.slew(); //slew synthesis filters first
		}
		for (int m = 1; m < N-1; ++m) {
			xtemp[m] = ytemp + capf[m][0] * apf[m].memL;
			ytemp = capf[m][1] * xtemp[m] + capf[m][2] * apf[m].memL;
			ytemp += Y[m]; //accumulate with phase
		}
		ytemp += Y[N - 1]; //last highpassed output is already phase aligned
		for (int m = 0; m < N-1; ++m) {
			apf[m].memL = xtemp[m];
		}
		return (float)ytemp;
	}
	void synthesizeBlock(const float*const* Y, float* y, int len)
	{
		bool slewing = apf[0].sc.check(); //all filter slews are coupled
		if (slewing) { //slew every sample
			if (transposed) {
				for (int i = 0; i < len; ++i) {
					y[i] = synthesize(Y[i]);
				}
			}
			else {
				for (int i = 0; i < len; ++i) {
					for (int n = 0; n < N; ++n) {
						scr[n] = Y[i][n];
					}
					y[i] = synthesize(scr);
				}
			}
		}
		else { //faster block-based application
			double xtemp[N-1];
			double ytemp = 0.0;
			if (transposed) {
				for (int i = 0; i < len; ++i) {
					ytemp = (double)Y[i][0];
					for (int m = 1; m < N-1; ++m) {
						xtemp[m] = ytemp + capf[m][0] * apf[m].memL;
						ytemp = capf[m][1] * xtemp[m] + capf[m][2] * apf[m].memL;
						ytemp += Y[i][m]; //accumulate with phase
					}
					ytemp += Y[i][N - 1]; //last highpassed output is already phase aligned
					for (int m = 0; m < N-1; ++m) {
						apf[m].memL = xtemp[m];
					}
					y[i] = (float)ytemp;
				}
			}
			else {
				for (int i = 0; i < len; ++i) {
					ytemp = (double)Y[0][i];
					for (int m = 1; m < N-1; ++m) {
						xtemp[m] = ytemp + capf[m][0] * apf[m].memL;
						ytemp = capf[m][1] * xtemp[m] + capf[m][2] * apf[m].memL;
						ytemp += Y[m][i]; //accumulate with phase
					}
					ytemp += Y[N - 1][i]; //last highpassed output is already phase aligned
						for (int m = 0; m < N-1; ++m) {
							apf[m].memL = xtemp[m];
						}
					y[i] = (float)ytemp;
				}
			}
		}
	}

private:
	sof lpf[N - 1]; //list of LPF stages
	sof hpf[N - 1]; //list of HPF stages
	fof apf[N - 1]; //list of APF stages

	const double* clpf[N - 1]; //preallocated list of pointers to slewed LPF coefs
	const double* chpf[N - 1]; //preallocated list of pointers to slewed HPF coefs
	const double* capf[N - 1]; //preallocated list of pointers to slewed APF coefs

	bool transposed = false;
	float scr[N]; //multiband transposition scratch
};

}