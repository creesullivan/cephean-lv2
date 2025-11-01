
#pragma once

//--------------------------------------------------
// Defines DSP classes and functions
//--------------------------------------------------

#include "cephean-lv2.h"
#include "cephean-filter.h"
#include "cephean-fft.h"

#include <cstdlib> //random numbers

namespace cephean
{

//==================================================

//Returns a random float value between 0 and 1 using the
//seed set by srand().
inline float randf()
{
	return ((float)rand()) / RAND_MAX;
}

//==================================================

//Periodic triangle wave from -1 to +1 as a function of phase x
//where triwave(0) = triwave(2*pi) = 0.0f ala a sine wave.
float triwave(float x);

//Periodic quadratic wave from -1 to +1 as a function of phase x
//where quadwave(0) = quadwave(2*pi) = 0.0f ala a sine wave.
float quadwave(float x);

//==================================================

//Circularly addressing integer on some range between 0 and a given 
//wrapping length. Helpful for circular buffers of data where increment
//and decrement operations are guaranteed to be withing +- one wrapping
//length to avoid a modulo operation.
class fastcircint
{
public:
	enum setlenmode
	{
		ZERO = 0,
		WRAP = 1,
		BOUND = 2
	};

	fastcircint(int len = 1, int vset = 0);
	fastcircint(const fastcircint& other);
	~fastcircint();

	void set(int vset);
	void setlength(int len, setlenmode mode = setlenmode::WRAP);
	void operator=(int vset);

	int length() const;
	int get() const;
	operator int() const;

	void operator+=(int i);
	void operator-=(int i);
	void operator+=(unsigned int i);
	void operator-=(unsigned int i);
	int operator++();
	int operator++(int);
	int operator--();
	int operator--(int);

	fastcircint operator+(int i) const;
	fastcircint operator-(int i) const;

private:
	int val = 0;
	int wrap = 1;
};

//Vector of fastcircint with a helpful constructor that automatically
//sets all owned objects to the same wrapping length
class fcivec : public vec<fastcircint>
{
public:
	fcivec(int len, int wrap=1);
	~fcivec();

	void setwrap(int wrap, fastcircint::setlenmode mode = fastcircint::setlenmode::WRAP);
};

//Circularly addressing integer on some range between 0 and a given 
//wrapping length that is equal to 2^r with r <= 16. This is similar
//to fastcircint, but has unlimited safety up to 16-bit overflow and
//should run faster due to high efficiency bitwise operations.
class pow2circint
{
public:
	enum setlenmode
	{
		ZERO = 0,
		WRAP = 1,
		BOUND = 2
	};

	pow2circint(unsigned int rset = 16, uint16_t vset = 0);
	pow2circint(const pow2circint& other);
	~pow2circint();

	void set(uint16_t vset);
	void setlength(unsigned int rset, setlenmode mode = setlenmode::WRAP);
	void operator=(uint16_t vset);

	uint16_t wrap(uint16_t v) const;
	void wrap(const uint16_t* x, uint16_t* y, int N) const;

	unsigned int length() const;
	uint16_t get() const;
	operator uint16_t() const;

	void operator+=(int i);
	void operator-=(int i);
	uint16_t operator++();
	uint16_t operator++(int);
	uint16_t operator--();
	uint16_t operator--(int);

	pow2circint operator+(int i) const;
	pow2circint operator-(int i) const;

private:
	uint16_t val = 0;
	unsigned int len = 65536;
	unsigned int r = 16;
	uint16_t wrapmask = 0xFFFF;
};

//Circularly wrapping float on some range between 0 and a given 
//wrapping length. Helpful for phase state and circular buffers of
//data with fractional indexing. This object is "fast" because it
//assumes all increment/decrement operations will wrap the value by
//at most one period.
class fastcircfloat
{
public:
	enum setwrapmode
	{
		ZERO = 0,
		WRAP = 1,
		BOUND = 2
	};

	fastcircfloat(float wrapset = 1.0f, float valset = 0.0f);
	fastcircfloat(const fastcircfloat& other);
	~fastcircfloat();

	void set(float newval);
	void setwrap(float newwrap, setwrapmode mode = setwrapmode::WRAP);
	void operator=(float newval);

	float getwrap() const;
	float get() const;
	operator float() const;

	//inc must be positive! use safeinc() for mixed sign
	void operator+=(float inc);

	//dec must be positive! use safeinc() for mixed sign
	void operator-=(float dec);

	//incdec may be positive (increment) or negative (decrement)
	void safeinc(float incdec);

	fastcircfloat operator+(float inc) const;
	fastcircfloat operator-(float inc) const;

private:
	float val = 0.0f;
	float wrap = 1.0f;
};

//Vector of fastcircfloat with a helpful constructor that automatically
//sets all owned objects to the same wrapping point
class fcfvec : public vec<fastcircfloat>
{
public:
	fcfvec(int len, float wrap = 1.0f);
	~fcfvec();

	void setwrap(float wrap, fastcircfloat::setwrapmode mode = fastcircfloat::setwrapmode::WRAP);
};

//==================================================

//Returns proclen + datalen - 1, the number of samples needed to apply
//proclen-length processing convolution style with datalen output samples
inline int getDataLengthFor(int proclen, int datalen)
{
	return proclen + datalen - 1;
}

//Mono data buffer with simple put/get functions and flat,
//non-circular storage.
class flatbuffer : private vec<float>
{
public:
	flatbuffer(int len=1);
	~flatbuffer();

	int size() const;
	void reset(int newlen);

	void clear(float val = 0.0f);

	void put(const float* x, int len);

	float get(int del = 0, int hostlen = 1, int hostn = 0) const;
	void get(float* y, int len, int del = 0) const;
	float* get(int len, int del = 0);
	const float* get(int len, int del = 0) const;
};

//Mono data buffer with put/set/get functions and circular
//storage, stepping the buffer on every put. 
class circbuffer : private vec<float>
{
public:
	circbuffer(int len=1);
	~circbuffer();

	int size() const;
	void reset(int newlen);

	void clear(float val = 0.0f);

	void put(float x);
	void put(const float* x, int len);

	void set(float x, int del = 0, int hostlen = 1, int hostn = 0);
	void set(const float* x, int len, int del = 0);

	float get(int del = 0, int hostlen = 1, int hostn = 0) const;
	void get(float* y, int len, int del = 0) const;

	float getNN(float del, int hostlen = 1, int hostn = 0) const;
	//interpolated from a starting delay, adjusting ddel samp/samp for pitch shifting
	void getNN(float* y, int len, float del0 = 0, float ddel = 0.0f, int hostlen = -1, int hostn = 0) const;
	void getNN(float* y, const float* del, int len) const;

	float getLinear(float del, int hostlen = 1, int hostn = 0) const;
	void getLinear(float* y, int len, float del, int hostlen = -1, int hostn = 0) const;

private:
	fastcircint ind; //points to most recent sample
};

//Mono data buffer with put/set/get functions and circular
//storage with power-of-two length for speed, stepping the
//buffer on every put. 
class pow2buffer : private vec<float>
{
public:
	pow2buffer(unsigned int rset = 16);
	~pow2buffer();

	int size() const;
	void reset(unsigned int newr);

	void clear(float val = 0.0f);

	void put(float x);
	void put(const float* x, int len);

	void set(float x, int del = 0, int hostlen = 1, int hostn = 0);
	void set(const float* x, int len, int del = 0);

	float get(int del = 0, int hostlen = 1, int hostn = 0) const;
	void get(float* y, int len, int del = 0) const;

	//interpolation functions TODO

private:
	pow2circint ind; //points to most recent sample
};

//Mono flat data buffer with an allpass + sample delay between
//each node
class warpbuffer : private vec<float>
{
public:
	warpbuffer(int slew = 1, int len = 1);
	~warpbuffer();

	int size() const;
	void reset(int newlen);

	//Returns alpha for a per node delay in samples, including the extra
	//single sample of delay. For stability, samples strictly >1.
	static float getAlphaForDelay(float samples);
	void setAlpha(float newAlpha, bool convergeInstantly=false);

	void clear();

	//Quick and easy all-in-one step function, a new data sample comes in
	//at the front or a single node, the buffer is shifted once, and data
	//is read from the back or another single node
	float step(float x, int xnode = 0, int ynode = -1);

private:
	slewed<float> alpha; //alpha constant of the warping
	vec<float> xnew; //input injection memory
	vec<float> bmem; //numerator memory vector
};

//==================================================

//Simple mono slewed gain that can be applied by sample or block
class gain : public monoalg
{
public:
	gain(int slew = 1);
	~gain();

	void setLevel(float mult, bool converge = false);

	//empty clear function
	float step(float x) override;
	void stepBlock(const float* x, float* y, int len) override;

private:
	slewed<float> g; //gain multiplier
};


//==================================================

//Asymmetric attack/release mono envelope detector
class attrel : public monoalg
{
public:
	attrel(float attackInSamples = 0.0f, float releaseInSamples = 0.0f);
	virtual ~attrel();

	void setAttack(float samples);
	void setRelease(float samples);

	void clear(float val = 0.0f);
	virtual float step(float x) override;
	//default stepBlock function

private:
	float aalpha, abeta, ralpha, rbeta; //attack and release time constants
	float env; //state memory
};

//Symmetric multichannel smoother
class vsmo : public multialg
{
public:
	vsmo(int Nch, bool blockTransposed = false, float smoothingInSamples = 0.0f);
	virtual ~vsmo();

	void setSmooth(float samples);

	void clear(float val = 0.0f);
	virtual void step(const float* x, float* y, int chan) override;
	//default stepBlock function

private:
	float alpha, beta; //fixed smoothing time constants
	fvec env; //state memory
};

//Asymmetric attack/release multichannel envelope detector with instant release
class vatt : public multialg
{
public:
	vatt(int Nch, bool blockTransposed = false, float attackInSamples = 0.0f);
	virtual ~vatt();

	void setAttack(float samples);

	void clear(float val = 0.0f);
	virtual void step(const float* x, float* y, int chan) override;
	//default stepBlock function

private:
	float aalpha, abeta; //fixed attack time constants
	fvec env; //state memory
};

//Asymmetric attack/release multichannel envelope detector with instant attack
class vrel : public multialg
{
public:
	vrel(int Nch, bool blockTransposed = false, float releaseInSamples = 0.0f);
	virtual ~vrel();

	void setRelease(float samples);

	void clear(float val = 0.0f);
	virtual void step(const float* x, float* y, int chan) override;
	//default stepBlock function

private:
	float ralpha, rbeta; //fixed release time constants
	fvec env; //state memory
};

//==================================================

//Leaky moving average filter closely approximating a sliding
//FIR rectangular window filter
class lma : public monoalg
{
public:
	lma(int slew = 1, double leakInSamples = 48000.0, int maxLengthInSamples = 4096, int maxBlockSize = 128);
	~lma();

	void setLength(int samples, bool converge = false);

	void clear();
	virtual float step(float x) override;
	virtual void stepBlock(const float* x, float* y, int len) override;

private:
	const double L = 48000.0;

	double alpha, gain, betaN; //y[n] = alpha*y[n-1] + x[n]; z[n] = gain*(y[n] - betaN*y[n+N]);
	double mem = 0.0; //feedback memory
	circbuffer buff; //output buffer
	slewed<int> N; //slewed window length
	fvec scr; //scratch block

	void update();
};

//4th order hold algorithm approximating a moving max filter,
//defined by its quarter-length and clear value.
class hold4 : public monoalg
{
public:
	hold4(float clearval = 0.0f);
	virtual ~hold4();

	//Set the quarter-length of the hold process, the full hold = 4*Nq
	void setQuarterLength(int Nq);

	void clear();
	virtual float step(float x) override;
	//default stepBlock function

private:
	const float x0 = 0.0f; //clear value

	float mem[4]; //hold value memory
	fastcircint n; //hold counter from 0 to N
};

//4th order hold algorithm like hold4, but applied to a
//multichannel vectorized block of data.
class vhold4 : public multialg
{
public:
	vhold4(int Nch, bool blockTransposed = false, float clearval = 0.0f);
	virtual ~vhold4();

	//Set the quarter-length of the hold process, the full hold = 4*Nq
	void setQuarterLength(int Nq);

	void clear();
	virtual void step(const float* x, float* y, int chan) override;
	//default stepBlock function

private:
	const float x0 = 0.0f; //clear value
	const int M = 1; //fixed number of channels

	fvec mem; //flat hold value memory arranged in blocks of length N
	fcivec n; //hold counters from 0 to N
};

//==================================================

//Fixed saturating waveshaper with unity gain at (0, 0)
//and a unity saturation point. Apply pre-gain to drive
//the saturator, and post-gain to adjust the clip point.
//Returns the gain to multiply by the input, not the final
//output. This permits gain manipulations.
class saturator : public monoalg
{
public:
	saturator();
	float step(float x) override;
	using monoalg::stepBlock; //default
};

//Fast ratio waveshaping nonlinearity with built-in slew.
//This is valid over a range from roughly -40 dB to 0 dB
//when undriven, below which it linearizes and above which
//it saturates. Drive gain is available to reposition the
//valid range, adjusted so that unity always maps to unity
//under ideal assumptions (approx errors may result).
class fastratio : public monoalg, public stereoalg
{
public:
	//applied y = x*(a1/(1 + b1*|x|) + a2/(1 + b2*|x|))
	struct coefs {
		float a1 = 0.0f;
		float a2 = 0.0f;
		float b1 = 0.0f;
		float b2 = 0.0f;
	};

	//stationary points of undriven SOLE
	static const float x0; //0.01f; 
	static const float x1; //0.1f;
	static const float x2; //0.5f;
	//x3 = 1.0f;

	static coefs design(float ratio, float drive);

	fastratio(int slew = 1);

	void set(float newRatio, float newDrive, bool converge = false);

	float step(float x) override;
	using monoalg::stepBlock; //default

	void step(float xL, float xR, float& yL, float& yR) override;
	using stereoalg::stepBlock; //default

	float stepGain(float x);
	void stepGainBlock(const float* x, float* g, int len);

	//TO DO <-- block optimizations! at least for mono version

private:
	fadedvec<float> prop; //ratio and drive
	coefs coefA, coefB;

	void update();
};

//Fast ratio waveshaper applied in 8 stages with first order
//lowpass filters in between to antialias simultaneously
class simshaper
{
public:
	simshaper(int slew = 1, int maxlen = 128);

	void set(float newRatio, float newDrive, float newRolloff, bool converge = false);
	void clear();

	//x == y is UNSAFE
	void stepBlock(const float* x, float* y, int len);

private:
	const int M = 8; //stage count

	fadedvec<float> prop; //ratio, drive, and rolloff

	fvec yscr, gscr; //scratch memory for crossfade output and gain

	//A group for fade
	float driveA = 1.0f; //input drive gain
	float makeupA = 1.0f; //output makeup gain
	fof::coefs lpfA; //shared LPF coefs
	fastratio::coefs wsA; //shared waveshaper coefs
	float memA[8]; //memory slots for DF-II LPFs

	//B group for fade
	float driveB = 1.0f; //input drive gain
	float makeupB = 1.0f; //output makeup gain
	fof::coefs lpfB; //shared LPF coefs
	fastratio::coefs wsB; //shared waveshaper coefs
	float memB[8]; //memory slots for DF-II LPFs

	void update();
};

//==================================================

//Linear interpolation at fractional sample index ind, forming
// y = x[ind]. Extrapolation here is safe.
inline float interp(const float* x, float ind, int xlen)
{
	if (ind <= 0.0f) {
		return x[0];
	}
	else if (ind >= (xlen - 1.0f)) {
		return x[xlen - 1];
	}
	else {
		float gain2 = floorf(ind);
		int ind1 = (int)gain2;
		gain2 = ind - gain2;
		return (1.0f - gain2) * x[ind1] + gain2 * x[ind1 + 1];
	}
}

//Linear interpolation at fractional sample indices ind, forming
// y = x[ind]. Extrapolation is unsafe! Pre-bound ind to [0, xlen-1].
//
// x == y is unsafe, but ind == y is safe.
void interp(const float* x, const float* ind, float* y, int ylen);

//Linear interpolation engine that slews and saves a mapping from
// x -> y and reapplies it to save computation.
class interpolator
{
public:
	interpolator(int ylen, int slew = 1);
	~interpolator();

	//Sets the index map so that y = x[index], len == ylen and the indices
	//should be bounded by 0 and xlen - 1 for safety
	void set(const float* newIndex, int len, bool converge = false);

	//Returns a pointer to the index map for direct modification. Converge
	//manually if desired with the converge() function.
	float* set();

	void slewpause(bool paused);
	void converge();

	//Applies the previously set index.  x == y is UNSAFE
	void apply(const float* x, float* y, int xlen, int ylen);

private:
	void update();

	const int N;

	slewedvec<float> ind;

	fvec g1, g2;
	ivec i1, i2;
};


//==================================================

//Upsamples an input signal with aliasing and full power compensation.
//Typically, you will want an antialiasing filter to go along with this.
//
//The step functions require an input data length and a requested output
//data length, often rate * xlen. For applications upsampling a previously
//downsampled signal, sometimes block sizes may change, so manually specifying
//the output size is helpful.
class upsampler
{
public:
	upsampler(int rate = 2);
	~upsampler();

	void setRate(int newRate);
	int getRate() const;

	void clear();

	void step(float x, float* y, int ylen);
	void stepBlock(const float* x, float* y, int xlen, int ylen);

private:
	int R = 2; //upsampling rate multiplier
	int r = 0; //current upsampling index < R
};

//Downsamples an input signal with aliasing and no power compensation.
//Typically, you will want an antialiasing filter to go along with this.
//
//The step functions require an input data length and return the output
//data length, sometimes 0. The output data length can never be more than
//ceil(len/rate).
class downsampler
{
public:
	downsampler(int rate = 2);
	~downsampler();

	void setRate(int newRate);
	int getRate() const;

	void clear();

	int step(float x, float* y);
	int stepBlock(const float* x, float* y, int xlen);

private:
	int R = 2; //downsampling rate multiplier
	int r = 0; //current downsampling index < R
};

//==================================================

//Simple mono waveguide reverb algorithm with a slewed duration
//control and direct access to 16 individual waveguide delays
//defining the timbre's fine structure.
class waveverb
{
public:
	waveverb(int slew = 1, int maxbsize = 128, int nominalRate = 1);
	virtual ~waveverb();

	//Set the -20 dB decay point in samples
	virtual void setDecay(float samples, bool converge = false);

	//Retrieves the 16-length waveguide delay pointer for direct adjustment
	int* getDelayPointer();

	//Set the seed value on [0, 15], defining the waveverb's timbre
	void setSeed(unsigned int newSeed);

	void clear();
	//no step(), aggressively block optimized
	void stepBlock(const float* x, float* y, int len);

protected:
	const int M = 16; //order (fixed to 16)
	const int rate = 1; //nominal sampling rate (1, 2, or 4 re 48k)
	const int maxdel0 = 4800; //maximum tap delay of any buffer @ usr1
	const int maxlen = 128; //maximum block size
	const int N = 4224; //shared buffer length

	float decay = 0.0f; //-20 dB decay duration control in samples

	ivec delay; //M length waveguide delay vector in samples
	slewedvec<float> alpha; //alpha parameter vector = 0.1^(-delay/decay)

	fvec dryscr; //dry scratch vector
	circbuffer buffer[16]; //N by M circular buffer
	fvec xscr[16]; //maxlen by M input scratch vector
	fvec sscr[4]; //maxlen by 4 summation scratch vector
	fvec fscr; //maxlen feedback scratch vector

	static const int del0[16][16]; //tables of constants
};

//Waveverb algorithm with only 9 voices instead of 16, improving
//efficiency significantly for cases where short or low quality
//reverb is acceptable
class miniverb
{
public:
	miniverb(int slew = 1, int maxbsize=128, int nominalRate = 1);
	~miniverb();

	//Set the -20 dB decay point in samples
	void setDecay(float samples, bool converge = false);

	//Retrieves the 9-length waveguide delay pointer for direct adjustment
	int* getDelayPointer();

	//Set the seed value on [0, 15], defining the miniverb's timbre
	void setSeed(unsigned int newSeed);
	
	void clear();
	//no step(), aggressively block optimized
	void stepBlock(const float* x, float* y, int len);

private:
	const int M = 9; //order (fixed to 9)
	const int rate = 1; //nominal sampling rate (1, 2, or 4 re 48k)
	const int maxdel0 = 4800; //maximum tap delay of any buffer @ usr1
	const int maxlen = 128; //maximum block size
	const int N = 4224; //shared buffer length

	float decay = 0.0f; //-20 dB decay duration control in samples

	ivec delay; //M length waveguide delay vector in samples
	slewedvec<float> alpha; //alpha parameter vector = 0.1^(-delay/decay)

	fvec dryscr; //dry scratch vector
	circbuffer buffer[9]; //N by M circular buffer
	fvec xscr[9]; //maxlen by M input scratch vector
	fvec sscr[3]; //maxlen by 3 summation scratch vector
	fvec fscr; //maxlen feedback scratch vector

	static const int del0[16][9]; //tables of constants
};

//Specialized miniverb with a dry/wet mix that enforces unity total
//power responses within 1 dB accuracy. This only permits seed-based
//adjustment, not direct-set delays.
class diffuser
{
public:
	diffuser(int slew = 1, int maxbsize = 128, int nominalRate = 1);
	~diffuser();

	//Set the -20 dB decay point in samples
	void setDecay(float samples, bool converge = false);

	//Set the seed value on [0, 15], defining the miniverb's timbre
	void setSeed(unsigned int newSeed);

	void clear();
	//no step(), aggressively block optimized
	void stepBlock(const float* x, float* y, int len);

private:
	const float drymix = 0.1f;
	const float r = 1.7f;
	static const float g[16];

	const int rate = 1; //nominal sampling rate (1, 2, or 4 re 48k)

	int seed = 0; //index on [0, 15] for power norm calculation
	float tdecay = 0.0f; //decay time in seconds for power norm calculation
	float getNormGain() const;

	fvec scratch; //for wet/dry mix

	miniverb core;
	gain level;
};

//==================================================

//Core phaser processing
template<unsigned int order> class phasercore
{
public:
	phasercore(int slew = 1, int alphaslew = 1) :
		salpha(expf(-1.0f / ((float)alphaslew))), sbeta(1.0f - salpha),
		fb(slew), depth(slew)
	{
		update();
		clear();
	}
	~phasercore() {}

	//Returns alpha for a per node delay in samples, for stability,
	//samples must be strictly >0
	static float getAlphaForDelay(float samples)
	{
		return (samples - 1.0f) / (samples + 1.0f);
	}
	void setDepth(float newDepth, float newFeedback, bool convergeImmediately = false)
	{
		fb.target(newFeedback, convergeImmediately);
		depth.target(newDepth, convergeImmediately);
		if (convergeImmediately) {
			update();
		}
	}
	void clear()
	{
		vset(mem, 0.0f, order);
		fbmem = 0.0f;
	}
	float step(float x, float alpha)
	{
		//slew first
		bool shouldUpdate = depth.slew();
		shouldUpdate = fb.slew() || shouldUpdate;
		if (shouldUpdate) {
			update();
		}
		amem = salpha * amem + sbeta * alpha;

		const float acur = amem;
		float xtemp[order];
		float ytemp = x + fb * fbmem; //feedback
		for (int n = 0; n < order; ++n) {
			xtemp[n] = ytemp + acur * mem[n];
			ytemp = -acur * xtemp[n] + mem[n];
		}
		for (int n = 0; n < order; ++n) {
			mem[n] = xtemp[n];
		}
		fbmem = ytemp;
		return (gdry * x - depth * fbmem);
	}

	void stepBlock(const float* x, float alpha, float* y, int len)
	{
		bool shouldUpdate = false;
		float xtemp[order];
		float ytemp = 0.0f;

		for (int i = 0; i < len; ++i) {
			shouldUpdate = depth.slew(); //slew first
			shouldUpdate = fb.slew() || shouldUpdate;
			if (shouldUpdate) {
				update();
			}
			amem = salpha * amem + sbeta * alpha;
			const float acur = amem;
			const float* memcur = mem;

			ytemp = x[i] + fb * fbmem; //feedback
			for (int n = 0; n < order; ++n) {
				xtemp[n] = ytemp + acur * memcur[n];
				ytemp = -acur * xtemp[n] + memcur[n];
			}
			for (int n = 0; n < order; ++n) {
				mem[n] = xtemp[n];
			}
			fbmem = ytemp;
			y[i] = gdry * x[i] - depth * fbmem;
		}
	}

private:
	const float salpha = 0.0f; //alpha smoothing constant
	const float sbeta = 1.0f;

	slewed<float> fb; //feedback gain
	slewed<float> depth; //depth gain
	float gdry = 0.0f; //dry gain
	float amem = 0.0f; //alpha smoothing memory

	float mem[order]; //filter DF-II memory
	float fbmem = 0.0f; //feedback memory

	void update() { //updates gdry from fb and depth
		gdry = (1.0f + (depth / (1.0f - fb)));
	}
};

//Core flanger processing
class flangercore
{
public:
	flangercore(int slew = 1, int delslew = 1, int maxdel = 128, int maxbsize = 128);
	~flangercore();

	void setDepth(float newDepth, float newFeedback, bool convergeImmediately = false);

	void clear();
	float step(float x, float del);
	void stepBlock(const float* x, float del, float* y, int len);

private:
	const float salpha = 0.0f; //delay smoothing constant
	const float sbeta = 1.0f;

	slewed<float> fb; //feedback gain
	slewed<float> depth; //depth gain
	float gdry = 0.0f; //dry gain for DC normalization
	float delmem = 0.0f; //delay smoothing memory

	circbuffer buff; //core memory buffer
	float fbmem = 0.0f; //feedback memory

	void update(); //updates gdry from fb and depth
};

/*
//Core ensemble chorus/doubler voice processor, passed a shared
//circbuffer for memory savings in multivoice systems
class doubler
{
public:
	doubler(const circbuffer& buffref, int maxbsize = 128);
	~doubler();

	void setMaxDelay(int newMaxDelay);
	void setMaxStep(float newMaxSampleStep);
	void setMaxGain(float newMaxGainMult);
	void setPeriodRange(int mewMinPeriod, int newMaxPeriod);

	void clear(int setDelay = 0, float setShift = 1.0f, float setGainMult = 1.0f);
	void stepBlock(float* y, int len);

private:
	const circbuffer& buff; //core shared memory buffer

	fastcircint l; //current index in the buffer
	
};
*/

//==================================================

//Mono downsampled flat buffer that lets you compute a
//normalized autocorrelation sequence. This is useful as
//the core for pitch-shifting and pitch-detecting processes.
class corrbuffer
{
public:
	friend class downshift;
	friend class upshift;

	corrbuffer(int len = 2048, int chunk = 448, int dsr = 8, int maxbsize = 128);
	~corrbuffer();

	int size() const;
	int dsize() const;
	int Rsize() const;

	//Power threshold for division stability
	void setThresh(float powerThreshold);

	//Correlation samples less than this are ignored by peak finding
	void setMinPeriod(int samples);

	void clear();

	void put(const float* x, int len);

	//Computes the downsampled normalized autocorrelation sequence and copies to R
	void corr(float* R);

	//Searches R over its valid range for the largest positive value
	float corrpeak(const float* R, int& delay, int mindelay = -1, int maxdelay = -1) const;

	//Searches R over its valid range for its lowest-delay local maximum at least
	//allow away from the global peak, optionally constraining the search range
	float firstpeak(const float* R, float& delay, float allow=0.95f, int mindelay = -1, int maxdelay=-1) const;

private:
	const int Nu = 2048;
	const int Nd = 256;
	const int Md = 56;
	const int DSR = 8;
	const int maxlen = 128;

	const int Rlen = 201;
	const int K = 127;

	float T = 1e-8f; //threshold parameter for division safety
	int minperd = 12; //first samples in R we ignore for finding the peak thresh

	fft F; //internal FFT
	downsampler D; //internal downsampler for fractional block management
	sofcasc<4> aa; //antialiasing lowpass filter
	//no phase matching, it should be good enough without it
	lma ave; //leaky moving average for chunk power sequence

	flatbuffer buffd; //buffer at the downsampled rate
	flatbuffer bpowd; //current downsampled buffer chunk power sequence

	cfvec fscr1, fscr2; //complex spectral scratch vectors
	fvec tscr; //time domain scratch vector
};

//Manages an internal circbuffer of data and performs a fixed
//pitch downshift by some sample step multiplier <1. Initialize
//with a reference to a corrbuffer object used for correlation
//detection. The corrbuffer must be stepped manually before
//the stepBlock() function so that its results can be shared
//among several pitch shift objects.
class downshift
{
public:
	downshift(float sampleStep, int blendLength, const corrbuffer& Robj);
	~downshift();

	//Sets the absolute correlation threshold that flags a transient and snaps to front
	void setSnapThresh(float thresh);
	//Sets the relative correlation threshold that flags a period to blend across
	void setBlendThresh(float thresh);

	void clear();

	void stepBlock(const float* x, float* y, int len, const float* R, float Rmax, int Rmaxdel);

private:
	const int Nu = 2048;
	const int DSR = 8;
	const int maxlen = 128;
	const int Rlen = 201;

	const float ddel = 0.0f; //sample/sample delay step >=0
	const int blendlen = 256;

	const int maxdeld; //maxdel / dsr
	const int maxdel; //end-of-buffer delay that triggers a force snap to front

	float Tsnap = 0.5f;
	float Tblend = 0.9f;

	circbuffer buff; //data buffer
	sofcasc<2> lpf; //4th order BW post-lowpass for artifact reduction

	faded<bool> blend; //blending crossfade helper
	float del1 = 0.0f; //fractional sample delay for primary blending voice
	float del2 = 0.0f; //fractional sample delay for secondary blending voice
	bool intransient = false; //transient state latch

	fvec scr; //scratch voice memory
};

//Manages an internal circbuffer of data and performs a fixed
//pitch upshift by some sample step multiplier >1. Initialize
//with a reference to a corrbuffer object used for correlation
//detection. The corrbuffer must be stepped manually before
//the stepBlock() function so that its results can be shared
//among several pitch shift objects.
class upshift
{
public:
	upshift(float sampleStep, int blendLength, const corrbuffer& Robj);
	~upshift();

	//Sets the absolute correlation threshold that flags a transient and snaps to front
	void setSnapThresh(float thresh);
	//Sets the relative correlation threshold that flags a period to blend across
	void setBlendThresh(float thresh);
	//Minimum sample shift corresponding roughly to highest pitch period
	void setMinPeriod(int samples);

	void clear();

	void stepBlock(const float* x, float* y, int len, const float* R, float Rmax, int Rmaxdel);

private:
	const int Nu = 2048;
	const int DSR = 8;
	const int maxlen = 128;
	const int Rlen = 201;

	const float ddel = 0.0f; //sample/sample delay step <=0
	const int blendlen = 256;

	const int mindeld; //mindel / dsr
	const int mindel; //snap to front delay and earliest delay that can be blended to

	float Tsnap = 0.5f;
	float Tblend = 0.9f;
	int minperd = 12; //minper / dsr
	int minper = 96; //minimum sample shift by a non-snapping blend

	circbuffer buff; //data buffer
	sofcasc<2> lpf; //4th order BW pre-lowpass for artifact reduction

	faded<bool> blend; //blending crossfade helper
	float del1 = 0.0f; //fractional sample delay for primary blending voice
	float del2 = 0.0f; //fractional sample delay for secondary blending voice

	fvec scr; //scratch voice memory
};

}