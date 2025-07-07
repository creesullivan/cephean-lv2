
#pragma once

//--------------------------------------------------
// Defines spectral classes and functions
//--------------------------------------------------

#include "cephean-lv2.h"
#include "cephean-filter.h"

#include "fftw3.h"

namespace cephean
{
	
//=====================================================
	
inline int getKForNpts(int npts)
{
	return npts/2 + 1;
}

//log-spaced frequency points to represent K linear spectral points
//with at least Noct octave resolution, 1/48 by default
int getKlog(int K, float Noct=1.0f / 48.0f);

//log-spaced frequency points to represent K linear spectral points
//over a range f1 to f2 in normalized frequency and at least Noct
//octave resolution, 1/48 by default
int getKlog(int K, float f1, float f2, float Noct = 1.0f / 48.0f);

//Remaps spectra X -> Y with linear interpolation at Ky log spaced values
//between f1 and f2, defaulting to original X range excluding DC. Get Ky
//from getKlog() to fix a certain octave resolution.
void lin2log(const float* X, float* Y, int Kx, int Ky, float f1=0.0f, float f2=1.0f);

//=====================================================

//Class that performs an FFT of fixed size using the fftw3
//library, both forwards and backwards
class fft
{
public:
	fft(int setNpts);
	~fft();

	int getNpts() const;

	int getK() const;

	//Applies a real->cplx FFT from npts -> K with automatic padding
	void rfwd(const float* td, complex<float>* fd, int len=-1);
	
	//Windowed version of rfwd
	void rfwd(const float* td, complex<float>* fd, const float* wdw, int len=-1);

	//Applies a cplx->real IFFT from K -> npts with automatic padding
	void rinv(const complex<float>* fd, float* td, int len=-1);
	
	//Windowed version of rinv
	void rinv(const complex<float>* fd, float* td, const float* wdw, int len=-1);

	//Applies a real->real magnitude FFT from npts->len with automatic padding
	void rmag(const float* td, float* fdmag, int len=-1);
	
	//Windowed version of rmag
	void rmag(const float* td, float* fdmag, const float* wdw, int len=-1);

	//2x time aliased windowed version of rmag (len == 2*npts)
	void rmagA2(const float* td, float* fdmag, const float* wdw);

	//Applies a real->real linear phase IFFT from K -> npts with optional
	//circular shifting (by default latency of npts/2)
	void rlpinv(const float* fd, float* td, int lat=-1);
	
	//Windowed version of rlpinv
	void rlpinv(const float* fd, float* td, const float* wdw, int lat=-1);

	//Linear phase transform of a spectral magnitude response over K samples
	//assuming a latency of lat samples, npts/2 by default. It is usually faster
	//to do a lp inverse, as the lp part can be applied as a time shift.
	void linphase(const float* Hmag, complex<float>* H, int lat=-1);

	//Minimum phase transform of a spectral magnitude response over K samples.
	//Best practices are to use spectral oversampling -- aliasing may cause
	//errors in the result for critically sampled magnitude spectra.
	void minphase(const float* Hmag, complex<float>* H);

	//Minimum phase transform of a spectral magnitude response into a half-
	//filter assuming we are operating at a critically sampled rate
	//void minphaseA2(const float* Hmag, complex<float>* H);
	
	//Applies a real->real minimum phase IFFT from K->len with automatic padding
	void rmpinv(const float* fd, float* td, int len = -1);

	//Applies a real->real 2x aliased minimum phase IFFT from K->len with
	//automatic padding
	//void rmpinvA2(const float* fd, float* td, int len = -1);

private:
	const int npts; //number of points of this transform
	const int K;
	
	fvec tdscr;
	cfvec fdscr;
	fftwf_plan fftp_fwd, fftp_inv;

	//fvec hkernA2; //hilbert kernel for 2x aliased operation
};


//=====================================================

void impz(const fof::coefs& c, float* h, int len, int shift = 0);

void impz(const sof::coefs& c, float* h, int len, int shift = 0);

template<unsigned int N> void impz(const typename sofcasc<N>::coefs& c, float* h, int len, int shift = 0)
{
	if (len > 0) {
		sofcasc<N> filt;
		filt.setCoefs(c, true);
		vset(h, 0.0f, len);
		if ((shift < len) && (shift >= 0)) {
			h[shift] = 1.0f;
			filt.stepBlock(h, h, len);
		}
		else if (shift < 0) {
			int rem = -shift;
			vset(h, 0.0f, len);
			h[0] = 1.f;
			while (rem > 0) {
				int curlen = min(rem, len);
				filt.stepBlock(h, h, rem);
				vset(h, 0.0f, rem);
				rem -= curlen;
			}
			filt.stepBlock(h, h, len);
		}
	}
}


//=====================================================

//periodic linear phase Hann window with npts/2 sample latency
void hann(float* w, int npts);

//minimum phase Hann window approximation from a critically damped
//second order lowpass filter's truncated impulse response
void mphann(float* w, int npts, int shift=0);


//=====================================================

//Applies frequency domain 2x Hann smoothing as a linear interpolation
//kernel in the frequency domain. X == Y is NOT supported! X and Y are
//assumed to be real magnitude spectra of a half-length FFT.
void hannsmooth(const float* X, float* Y, int K);

//Fast spectral magnitude/power smoothing object with several
//available algorithms for different tasks.
class specsmooth
{
public:
	specsmooth(int setK, float setoct = (1.0f / 3.0f), float f0 = 0.0f);
	~specsmooth();

	//This invalidates any previously set weighting
	void setResolution(float newoct, float newf0 = 0.0f);

	//This invalidates any previously set weighting
	void setResolution(float newoct, float newfmin, float newfmax);

	//This invalidates any previously set weighting
	void setAsymResolution(float newoctup, float newoctdown,
		float newf0up = 0.0f, float newf0down = -1.0f);

	//sets the weighting vector used by wtsmooth()
	void setWeighting(const float* wt, int len);

	//sets the weighting vector to all 1s in the passband
	void setUnityWeighting();

	//linear symmetric smoothing with out-of-band decay
	void smooth(const float* X, float* Y, int len) const;

	//fast weighted moment smoothing with previously set weights in passband
	void wtsmooth(const float* X, float* Y, int len) const;

	//nonlinear masking-based symmetric smoothing
	void maxsmooth(const float* X, float* Y, int len) const;

private:
	fvec alphau, betau, alphad, betad; //vectors of subband smoothing feedback coefficients

	fvec w; //weight multiplier
	fvec wden; //precomputed smoothed weight denominator
};

//Complex-valued symmetric multichannel smoother
class cvsmo
{
public:
	cvsmo(int Nch, float smoothingInSamples = 0.0f);
	~cvsmo();

	void setSmooth(float samples);

	void clear(complex<float> val = 0.0f);
	void step(const complex<float>* x, complex<float>* y, int chan);

private:
	float alpha, beta; //fixed smoothing time constants
	cfvec env; //state memory
};

//=====================================================

//FFT enhanced FIR filter operating in the strict frequency domain
//via OLS with built-in sample-accurate slew
class fastfir
{
public:
	//FFT object must be have npts >= taps + bsize - 1
	fastfir(int taps, fft& Fref, int slew = 1);
	~fastfir();

	//sets the FIR filter via time domain coefficients
	void set(const float* hset, int hlen, bool converge = false);

	//sets the FIR filter via complex frequency domain coefficients
	void set(const complex<float>* Hset, int Kset, bool converge = false);

	//returns FIR filter complex spectral memory for direct modification
	complex<float>* set();

	void converge();

	//applies the filter to time domain buffer xbuff of length xlen,
	//generating y of length ylen. xlen >= taps + ylen - 1
	void apply(const float* xbuff, float* y, int xlen, int ylen);

private:
	fft& F; //FFT reference
	const int npts;
	const int K;

	fadedvec<complex<float>> H; //complex spectral FIR crossfade

	cfvec X; //complex spectral scratch for input data
	cfvec Y; //complex spectral scratch for output data
	fvec tscratch; //time domain scratch
};

}