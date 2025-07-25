
//--------------------------------------------------
// Defines spectral classes and functions
//--------------------------------------------------

#include "cephean-fft.h"

namespace cephean
{

//==================================================

//log-spaced frequency points to represent K linear spectral points
//with at least Noct octave resolution
int getKlog(int K, float Noct)
{
	return (int)ceilf(log2f(max(K, 1) - 1.0f) / Noct) + 1;
}

//log-spaced frequency points to represent K linear spectral points
//over a range f1 to f2 in normalized frequency and at least Noct
//octave resolution, 1/48 by default
int getKlog(int K, float f1, float f2, float Noct)
{
	f1 = bound(f1 * (K - 1), 1.0f, K - 1.0f);
	f2 = bound(f2 * (K - 1), f1, K - 1.0f);
	return (int)ceilf(log2f(f2 / f1) / Noct) + 1;
}

//Remaps spectra X -> Y with linear interpolation at Ky log spaced values
//between f1 and f2, defaulting to original X range excluding DC. Get Ky
//from getKlog() to fix a certain octave resolution.
void lin2log(const float* X, float* Y, int Kx, int Ky, float f1, float f2)
{
	f1 = bound(f1 * (Kx - 1), 1.0f, Kx - 1.0f); //convert to indexes
	f2 = bound(f2 * (Kx - 1), f1, Kx - 1.0f);

	const float logf1 = logf(f1);
	const float logf2 = logf(f2);
	float kxf = 0.0f;
	int kx1 = 0;
	int kx2 = 0;
	for (int ky = 0; ky < Ky; ++ky) {
		//get fractional index
		kxf = bound(expf(linmap((float)ky, 0.0f, Ky - 1.0f, logf1, logf2)), f1, f2);

		kx1 = (int)floorf(kxf); //linear interpolation
		kx2 = (int)ceilf(kxf);
		kxf = kxf - (float)kx1;
		Y[ky] = X[kx1] * (1.0f - kxf) + X[kx2] * kxf;
	}
}

//==================================================

fft::fft(int setNpts) : npts(setNpts), K(setNpts/2+1),
	tdscr(setNpts), fdscr(setNpts/2+1)//, hkernA2(setNpts)
{
	fftp_fwd = fftwf_plan_dft_r2c_1d(npts,
		tdscr.ptr(),
		reinterpret_cast<fftwf_complex*>(fdscr.ptr()),
		FFTW_ESTIMATE);
	fftp_inv = fftwf_plan_dft_c2r_1d(npts,
		reinterpret_cast<fftwf_complex*>(fdscr.ptr()),
		tdscr.ptr(),
		FFTW_ESTIMATE);
	/*
	hkernA2[0] = 0.0f; //embeds sqrt and norm operation in cepstral space
	for (int i = 1; i < npts; ++i) {
		hkernA2[i] = (0.5f / npts) * cosf((constants.pi * i) / npts);
	}
	*/
}

fft::~fft()
{
	fftwf_destroy_plan(fftp_fwd);
	fftwf_destroy_plan(fftp_inv);
}

//Returns the npts of the transform
int fft::getNpts() const
{
	return npts;
}

//Returns npts/2+1, the number of unique spectral indices resulting from
//a real forward transform
int fft::getK() const
{
	return K;
}

//Applies a real->cplx FFT from npts -> K with automatic padding
void fft::rfwd(const float* td, complex<float>* fd, int len)
{
	len = (len < 0 ? npts : len);
	vpadcopy(td, tdscr.ptr(), 0.0f, len, npts);
	fftwf_execute(fftp_fwd); //run FFT core
	vcopy(fdscr.ptr(), fd, K);
}

//Windowed version of rfwd
void fft::rfwd(const float* td, complex<float>* fd, const float* wdw, int len)
{
	len = (len < 0 ? npts : len);
	vpadcopy(td, tdscr.ptr(), 0.0f, len, npts);
	vmult(tdscr.ptr(), wdw, tdscr.ptr(), len); //window
	fftwf_execute(fftp_fwd); //run FFT core
	vcopy(fdscr.ptr(), fd, K);
}

//Applies a cplx->real IFFT from K -> npts with automatic padding
void fft::rinv(const std::complex<float>* fd, float* td, int len)
{
	len = (len < 0 ? npts : len);
	vcopy(fd, fdscr.ptr(), K);
	fftwf_execute(fftp_inv); //run IFFT core
	vmult(tdscr.ptr(), 1.0f / npts, tdscr.ptr(), npts); //normalize
	vpadcopy(tdscr.ptr(), td, 0.0f, npts, len);
}

//Windowed version of rinv
void fft::rinv(const complex<float>* fd, float* td, const float* wdw, int len)
{
	len = (len < 0 ? npts : len);
	vcopy(fd, fdscr.ptr(), K);
	fftwf_execute(fftp_inv); //run IFFT core
	vmult(tdscr.ptr(), 1.0f / npts, tdscr.ptr(), npts); //normalize
	vpadcopy(tdscr.ptr(), td, 0.0f, npts, len);
	vmult(td, wdw, td, len); //window
}

//Applies a real->real magnitude FFT from npts->K with automatic padding
void fft::rmag(const float* td, float* fdmag, int len)
{
	len = (len < 0 ? npts : len);
	vpadcopy(td, tdscr.ptr(), 0.0f, len, npts);
	fftwf_execute(fftp_fwd); //run FFT core
	vabs(fdscr.ptr(), fdmag, K); //take magnitude on the way out
}

//Windowed version of rmag
void fft::rmag(const float* td, float* fdmag, const float* wdw, int len)
{
	len = (len < 0 ? npts : len);
	vpadcopy(td, tdscr.ptr(), 0.0f, len, npts);
	vmult(tdscr.ptr(), wdw, tdscr.ptr(), len); //window
	fftwf_execute(fftp_fwd); //run FFT core
	vabs(fdscr.ptr(), fdmag, K); //take magnitude on the way out
}

//2x time aliased windowed version of rmag (len == 2*npts)
void fft::rmagA2(const float* td, float* fdmag, const float* wdw)
{
	vmult(td, wdw, tdscr.ptr(), npts);
	vmultaccum(td + npts, wdw + npts, tdscr.ptr(), npts);
	fftwf_execute(fftp_fwd); //run FFT core
	vabs(fdscr.ptr(), fdmag, K); //take magnitude on the way out
}

//Applies a real->real linear phase IFFT from K -> npts with optional
//circular shifting (by default latency of npts/2)
void fft::rlpinv(const float* fd, float* td, int lat)
{
	lat = (lat < 0 ? npts / 2 : lat);
	vcopy(fd, fdscr.ptr(), K); //copy data in as a 0 phase spectrum
	fftwf_execute(fftp_inv); //run IFFT core
	vmult(tdscr.ptr(), 1.0f / npts, tdscr.ptr(), npts); //normalize
	vcirccopy(tdscr.ptr(), td, npts, lat, npts); //circ shift on the way out
}

//Windowed version of rlpinv
void fft::rlpinv(const float* fd, float* td, const float* wdw, int lat)
{
	lat = (lat < 0 ? npts / 2 : lat);
	vcopy(fd, fdscr.ptr(), K); //copy data in as a 0 phase spectrum
	fftwf_execute(fftp_inv); //run IFFT core
	vmult(tdscr.ptr(), 1.0f / npts, tdscr.ptr(), npts); //normalize
	vcirccopy(tdscr.ptr(), td, npts, lat, npts); //circ shift on the way out
	vmult(td, wdw, td, npts); //window
}

//Linear phase transform of a spectral magnitude response over K samples
//assuming a latency of lat samples, npts/2 by default. It is usually faster
//to do a lp inverse, as the lp part can be applied as a time shift.
void fft::linphase(const float* Hmag, complex<float>* H, int lat)
{
	lat = (lat < 0 ? npts / 2 : lat);
	vincspace(tdscr.ptr(), 0.0f, -(2.0f * constants.pi * lat) / npts, K); //LP angle
	for (int i = 0; i < K; ++i) {
		H[i] = expi(tdscr[i]); //convert angle to complex spectra
		H[i] *= Hmag[i]; //embed magnitude
	}
}

//Minimum phase transform of a spectral magnitude response over K samples.
//Best practices are to use spectral oversampling -- aliasing may cause
//errors in the result for critically sampled magnitude spectra.
void fft::minphase(const float* Hmag, complex<float>* H)
{
	for (int i = 0; i < K; ++i) {
		fdscr[i] = logf(Hmag[i]);
	}
	fftwf_execute(fftp_inv); //cast to cepstrum

	//apply Hilbert kernel
	tdscr[0] = 0.0f; //clear DC
	vmult(tdscr.ptr() + 1, 1.0f / npts, tdscr.ptr() + 1, K - 2); //norm front-half
	tdscr[K - 1] = 0.0f; //clear NQ
	vmult(tdscr.ptr() + K, -1.0f / npts, tdscr.ptr() + K, K - 2); //invert norm back-half

	fftwf_execute(fftp_fwd); //cast to MP angle
	for (int i = 0; i < K; ++i) {
		H[i] = expi(fdscr[i].imag()); //convert angle to complex spectra
		H[i] *= Hmag[i]; //embed magnitude
	}
}
/*
//Minimum phase transform of a spectral magnitude response into a half-
//filter assuming we are operating at a critically sampled rate
void fft::minphaseA2(const float* Hmag, complex<float>* H)
{
	for (int i = 0; i < K; ++i) {
		fdscr[i] = logf(Hmag[i]);
	}
	fftwf_execute(fftp_inv); //cast to cepstrum

	vmult(tdscr.ptr(), hkernA2.ptr(), tdscr.ptr(), npts); //apply hilbert kernel

	fftwf_execute(fftp_fwd); //cast to MP angle
	for (int i = 0; i < K; ++i) {
		H[i] = expi(fdscr[i].imag()); //form net complex spectra
		H[i] *= sqrtf(Hmag[i]);
	}
}
*/
//Applies a real->real minimum phase IFFT from K->len with automatic padding
void fft::rmpinv(const float* fd, float* td, int len)
{
	len = (len < 0 ? npts : len);
	minphase(fd, fdscr.ptr()); //yeah this is safe
	fftwf_execute(fftp_inv); //run IFFT core
	vmult(tdscr.ptr(), 1.0f / npts, tdscr.ptr(), npts); //normalize
	vpadcopy(tdscr.ptr(), td, 0.0f, npts, len);
}
/*
//Applies a real->real 2x aliased minimum phase IFFT from K->len with
//automatic padding
void fft::rmpinvA2(const float* fd, float* td, int len)
{
	len = (len < 0 ? npts : len);
	minphaseA2(fd, fdscr.ptr()); //yeah this is safe
	fftwf_execute(fftp_inv); //run IFFT core
	vmult(tdscr.ptr(), 1.0f / npts, tdscr.ptr(), npts); //normalize
	vpadcopy(tdscr.ptr(), td, 0.0f, npts, len);
}
*/

//=====================================================

void impz(const fof::coefs& c, float* h, int len, int shift)
{
	if (len > 0) {
		fof filt;
		filt.setCoefs(c, true);
		vset(h, 0.0f, len);
		if ((shift < len) && (shift >= 0)) {
			h[shift] = 1.0f;
			filt.stepBlock(h, h, len);
		}
		else if(shift < 0){
			int rem = -shift;
			vset(h, 0.0f, len);
			h[0] = 1.f;
			while (rem > 0) {
				int curlen = min(rem, len);
				filt.stepBlock(h, h, rem);
				vset(h, 0.0f, len);
				rem -= curlen;
			}
			filt.stepBlock(h, h, len);
		}
	}
}

void impz(const sof::coefs& c, float* h, int len, int shift)
{
	if (len > 0) {
		sof filt;
		filt.setCoefs(c, true);
		vset(h, 0.0f, len);
		if((shift < len) && (shift >= 0)){
			h[shift] = 1.0f;
			filt.stepBlock(h, h, len);
		}
		else if(shift < 0) {
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

void freqz(sof::coefs c, const float* f, complex<float>* H, int K)
{
	complex<double> E;
	complex<double> EE;
	for (int k = 0; k < K; ++k) {
		E = expi(-constants.pi * f[k]);
		EE = E * E;
		H[k] = (c.b0 + c.b1 * E + c.b2 * EE);
		H[k] /= (1.0 - c.na1 * E - c.na2 * EE);
	}
}


//=====================================================

//periodic linear phase Hann window with npts/2 sample latency
void hann(float* w, int npts)
{
	const float wmult = 2.0f * constants.pi / npts;
	for (int i = 0; i < npts; ++i) {
		w[i] = 0.5f - 0.5f * cosf(wmult * i);
	}
}

//minimum phase Hann window approximation from a critically damped
//second order lowpass filter's truncated impulse response
void mphann(float* w, int npts, int shift)
{
	fvec scratch(npts);
	impz(filt::lowpass2(2.0f / npts, 0.5f), scratch.ptr(), npts, -shift);
	vrevcopy(scratch.ptr(), w, npts);
	float wnorm = mean(w, npts);
	vmult(w, 0.5f / wnorm, w, npts); //normalize DC to 0.5 like a Hann window
}

//=====================================================

//Applies frequency domain 2x Hann smoothing as a linear interpolation
//kernel in the frequency domain. X == Y is NOT supported! X and Y are
//assumed to be real magnitude spectra of a half-length FFT.
void hannsmooth(const float* X, float* Y, int K)
{
	vadd(X, X + 2, Y + 1, K - 2);
	vmult(Y + 1, 0.5f, Y + 1, K - 2);
	Y[0] = X[1];
	Y[K - 1] = X[K - 2];
	vadd(X, Y, Y, K);
	vmult(Y, 0.5f, Y, K);
}

specsmooth::specsmooth(int setK, float setoct, float setf0) :
	alphau(setK), betau(setK), alphad(setK), betad(setK), w(setK), wden(setK)
{
	setResolution(setoct, setf0); //init
	setUnityWeighting();
}
specsmooth::~specsmooth() {}

//This invalidates any previously set weighting
void specsmooth::setResolution(float newoct, float newf0)
{
	setAsymResolution(newoct, newoct, newf0, newf0);
}

//This invalidates any previously set weighting
void specsmooth::setResolution(float newoct, float newfmin, float newfmax)
{
	const int K = alphau.size();
	const float df = 1.0f / (K - 1.0f); //norm freq per spectral bin
	const float octmult = powf(2.0f, newoct / 2.0f) - 1.0f; //half octave mult
	newfmin /= 2.0f; //half linear bandwidth
	newfmax /= 2.0f; //half linear bandwidth

	float bw;
	for (int k = 0; k < K; ++k) {
		bw = bound(k * df * octmult, newfmin, newfmax); //norm freq half-bandwidth of this bin
		alphau[k] = expf(-df / bw); //feedback constants
		betau[k] = 1.0f - alphau[k];
	}
	vcopy(alphau.ptr(), alphad.ptr(), K);
	vcopy(betau.ptr(), betad.ptr(), K); //symmetric
}

//This invalidates any previously set weighting
void specsmooth::setAsymResolution(float newoctup, float newoctdown,
	float newf0up, float newf0down)
{
	const int K = alphau.size();
	const float df = 1.0f / (K - 1.0f); //norm freq per spectral bin
	const float octmultup = powf(2.0f, newoctup / 2.0f) - 1.0f; //half octave mult
	const float octmultdown = powf(2.0f, newoctdown / 2.0f) - 1.0f; //half octave mult
	newf0up /= 2.0f; //half linear bandwidth
	newf0down = (newf0down < 0.0f) ? newf0up : (newf0down / 2.0f);

	float bw;
	for (int k = 0; k < K; ++k) {
		bw = k * df * octmultup + newf0up; //norm freq half-bandwidth of this bin
		alphau[k] = expf(-df / bw); //feedback constants
		betau[k] = 1.0f - alphau[k];
	}
	for (int k = 0; k < K; ++k) {
		bw = k * df * octmultdown + newf0down; //norm freq half-bandwidth of this bin
		alphad[k] = expf(-df / bw); //feedback constants
		betad[k] = 1.0f - alphad[k];
	}
}

//sets the weighting vector to all 1s in the passband
void specsmooth::setUnityWeighting()
{
	int K = alphau.size();
	vset(w.ptr(), 1.0f, w.size());
	setWeighting(w.ptr(), K);
}

//sets the weighting vector used by wtsmooth() with no wt argument
void specsmooth::setWeighting(const float* wt, int len)
{
	const int K = alphau.size();
	vpadcopy(wt, w.ptr(), 0.0f, len, K);

	smooth(w.ptr(), wden.ptr(), K);

	//stabilize
	vadd(wden.ptr(), constants.eps, wden.ptr(), K);
}

//linear symmetric smoothing with out-of-band decay
void specsmooth::smooth(const float* X, float* Y, int len) const
{
	UNUSED(len); //len better equal K
	const int K = alphau.size();

	//upwards pass
	float mem = 0.0f;
	for (int k = 0; k < K; ++k) {
		mem = alphau[k] * mem + betau[k] * X[k];
		Y[k] = mem;
	}

	//downwards pass
	mem *= (betad[K - 1] * alphau[K - 1]) / (1 - alphau[K - 1] * alphad[K - 1]);
	for (int k = (K - 1); k >= 0; --k) {
		mem = alphad[k] * mem + betad[k] * Y[k];
		Y[k] = mem;
	}
}

//fast weighted moment smoothing with previously set weights in passband
void specsmooth::wtsmooth(const float* X, float* Y, int len) const
{
	//apply weights
	vmult(X, w.ptr(), Y, len);

	//apply smoothing
	smooth(Y, Y, len);

	//apply fixed weighting
	vdiv(Y, wden.ptr(), Y, len);
}

//nonlinear masking-based symmetric smoothing
void specsmooth::maxsmooth(const float* X, float* Y, int len) const
{
	UNUSED(len); //len better equal K
	const int K = alphau.size();

	//upwards pass
	float mem = 0.0f;
	for (int k = 0; k < K; ++k) {
		mem = max(X[k], alphau[k] * mem);
		Y[k] = mem;
	}

	//downwards pass
	mem *= alphau[K - 1];
	for (int k = (K - 1); k >= 0; --k) {
		mem = max(Y[k], alphad[k] * mem);
		Y[k] = mem;
	}
}

cvsmo::cvsmo(int Nch, float smoothingInSamples) : env(Nch)
{
	setSmooth(smoothingInSamples); //init
	clear();
}
cvsmo::~cvsmo() {}

void cvsmo::setSmooth(float samples)
{
	alpha = (samples > 0.0f) ? expf(-1.0f / samples) : 0.0f;
	beta = 1.0f - alpha;
}

void cvsmo::clear(complex<float> val)
{
	vset(env.ptr(), val, env.size());
}
void cvsmo::step(const complex<float>* x, complex<float>* y, int chan)
{
	for (int i = 0; i < chan; ++i) {
		env[i] = env[i] * alpha + beta * x[i];
	}
	vcopy(env.ptr(), y, chan);
}

//=====================================================

fastfir::fastfir(int taps, fft& Fref, int slew) :
	F(Fref), npts(Fref.getNpts()), K(Fref.getK()),
	H(Fref.getK(), slew), X(Fref.getK()), Y(Fref.getK()),
	tscratch(Fref.getNpts())
{
	vset(H.target(), 0.0f, K);
	H.converge(); //init to null filter
}
fastfir::~fastfir() {}

//sets the FIR filter via time domain coefficients
void fastfir::set(const float* hset, int hlen, bool converge)
{
	F.rfwd(hset, H.target(), hlen);
	if (converge) {
		H.converge();
	}
}

//sets the FIR filter via complex frequency domain coefficients
void fastfir::set(const complex<float>* Hset, int Kset, bool converge)
{
	UNUSED(Kset); //K better == Kset
	vcopy(Hset, H.target(), K); 
	if (converge) {
		H.converge();
	}
}

//returns FIR filter complex spectral memory for direct modification
complex<float>* fastfir::set()
{
	return H.target();
}

void fastfir::converge()
{
	H.converge();
}

//applies the filter to time domain buffer xbuff of length xlen,
//generating y of length ylen. xlen >= taps + ylen - 1
void fastfir::apply(const float* xbuff, float* y, int xlen, int ylen)
{
	H.swap(); //permit changes to the data at the top of the block

	F.rfwd(xbuff, X.ptr(), xlen); //fast conv number 1
	vmult(X.ptr(), H.current(), Y.ptr(), K);
	F.rinv(Y.ptr(), tscratch.ptr(), npts);
	vcopy(tscratch.ptr() + (xlen - ylen), y, ylen);

	if (H.check()) { //if cross-fading, run another fast conv and fade
		vmult(X.ptr(), H.last(), Y.ptr(), K);
		F.rinv(Y.ptr(), tscratch.ptr(), npts);

		float* pscr = tscratch.ptr() + (xlen - ylen);
		for (int i = 0; i < ylen; ++i) {
			H.fade();
			y[i] *= H.gain();
			y[i] += pscr[i] * (1.0f - H.gain());
		}
	}
}

} //namespace cephean

