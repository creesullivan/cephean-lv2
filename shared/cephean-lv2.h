
#pragma once

//--------------------------------------------------
// Defines classes and functions helpful for
// creating LV2 plugins and other very high level
// general purpose stuff
//--------------------------------------------------

#include <lv2/core/lv2.h>
#include <math.h>
#include <complex>
#include <cmath>
#include <stdint.h>
#include <stdlib.h>

#include "config/cephean-config.h"

#define UNUSED(x) (x)

namespace cephean
{

//==================================================

//Scoped denormal number disabler utility
class disableDenormals
{
public:
	disableDenormals();
	~disableDenormals();
private:
#if ISWINDOWS
	int MXCSR;
#endif
#if ISARM
	fenv_t fenv;
#endif
};

template<class ptype> using complex = std::complex<ptype>;

using std::abs;
using std::norm;
using std::arg;
using std::conj;

template<class ptype> inline complex<ptype> expi(ptype ph)
{
	return complex<ptype>(cos(ph), sin(ph));
}

//==================================================

//Incredibly simple RA-II vector of a type with a fixed length
//set on construction or destructively via reset
template<class ptype> class vec
{
public:
	vec(int len = 0) : N(len)
	{
		if (N > 0) {
			p = new ptype[N];
		}
	}
	vec(const ptype* data, int len) : vec(len)
	{
		for (int i = 0; i < N; ++i) {
			p[i] = data[i]; //deep copy
		}
	}
	vec(const vec& other) : vec(other.p, other.N) {}

	virtual ~vec()
	{
		if (N > 0) {
			delete[] p;
		}
	}

	void reset(int newlen)
	{
		if (N != newlen){
			if (N > 0) {
				delete[] p;
			}
			N = newlen;
			if (N > 0) {
				p = new ptype[N];
			}
		}
	}

	void operator=(const vec& other)
	{
		reset(other.N);
		for (int i = 0; i < N; ++i) {
			p[i] = other.p[i]; //deep copy
		}
	}

	int size() const { return N; }
	ptype* ptr() { return p; }
	const ptype* ptr() const { return p; }
	operator ptype* () { return p; }
	operator const ptype* () const { return p; }

protected:
	int N = 0;
	ptype* p = nullptr;
};
typedef vec<float> fvec;
typedef vec<int> ivec;
typedef vec<double> dvec;
typedef vec<complex<float>> cfvec;
typedef vec<complex<double>> cdvec;

//==================================================

//Incredibly simple RA-II 2D matrix, formed as a flat vector
//with built-in cast to pointers-of-pointers
template<class ptype> class mat
{
public:
	mat(int len1 = 0, int len2 = 0) : N(len1), M(len2), L(len1*len2)
	{
		if (L > 0){
			p = new ptype[L];
			pmat = new ptype * [N];
			for (int i = 0; i < N; ++i) {
				pmat[i] = (p + i * M);
			}
		}
	}
	mat(const ptype* data, int len1, int len2, bool transpose=false) :
		mat(transpose ? len2 : len1, transpose ? len1 : len2)
	{
		if (transpose) {
			for (int i = 0; i < N; ++i) {
				for (int j = 0; j < M; ++j) {
					pmat[i][j] = data[N*j + i]; //deep transposed copy
				}
			}
		}
		else {
			for (int l = 0; l < L; ++l) {
				p[l] = data[l]; //deep copy
			}
		}
	}
	mat(const ptype*const* data, int len1, int len2, bool transpose=false) :
		mat(transpose ? len2 : len1, transpose ? len1 : len2)
	{
		if (transpose) {
			for (int i = 0; i < N; ++i) {
				for (int j = 0; j < M; ++j) {
					pmat[i][j] = data[j][i]; //deep transposed copy
				}
			}
		}
		else {
			for (int i = 0; i < N; ++i) {
				for (int j = 0; j < M; ++j) {
					pmat[i][j] = data[i][j]; //deep copy
				}
			}
		}
	}
	mat(const mat& other, bool transpose=false) : mat(other.p, other.N, other.M, transpose) {}

	virtual ~mat()
	{
		if (L > 0) {
			delete[] p;
			delete[] pmat;
		}
	}

	void reset(int newlen1, int newlen2)
	{
		if ((N != newlen1) || (M != newlen2)) {
			if (L > 0) {
				delete[] p;
				delete[] pmat;
			}
			N = newlen1;
			M = newlen2;
			L = N * M;
			p = new ptype[L];
			pmat = new ptype * [N];
			for (int i = 0; i < N; ++i) {
				pmat[i] = (p + i * M);
			}
		}
	}

	void operator=(const mat& other)
	{
		reset(other.N, other.M);
		for (int l = 0; l < L; ++l) {
			p[l] = other.p[l]; //deep copy
		}
	}

	int size() const { return L; }
	void size(int& getN, int& getM) { getN = N; getM = M; }
	int cols() const { return N; }
	int rows() const { return M; }

	ptype* flatptr() { return p; }
	const ptype* flatptr() const { return p; }
	//no default cast to a flat pointer for clarity

	ptype* const* ptr() { return pmat; }
	const ptype* const* ptr() const { return pmat; }
	operator ptype* const* () { return pmat; }
	operator const ptype* const* () const { return pmat; }
	
protected:
	int L = 0;
	int N = 0;
	int M = 0;
	ptype* p = nullptr;
	ptype** pmat = nullptr;
};
typedef mat<float> fmat;
typedef mat<int> imat;
typedef mat<double> dmat;
typedef mat<complex<float>> cfmat;
typedef mat<complex<double>> cdmat;

//==================================================

//Core abstract class that defines the functions eventually
//exported via the LV2 descriptor
class plugin
{
public:
	plugin(int numControls, int numInputs, int numOutputs, 
		float setFs = 48000.0f,
		const char* bundle_path = nullptr,
		const LV2_Feature* const* features = nullptr);

	virtual ~plugin();

	void connect_port(uint32_t port, void* data);

	void commit();

	virtual void init();

	virtual void step(int len) = 0;

	virtual void deactivate();

	class control
	{
	public:
		control(int portIndex);
		virtual ~control();
		int getIndex() const;

		virtual void connect(void* data) = 0;
		virtual void commit() = 0;
		bool isnew() const;
		void clearnew();

	protected:
		const int index = 0;
		bool vnew = false;
	};

	template<class ptype> class tcontrol : public control
	{
	public:
		tcontrol(plugin* host, int portIndex) : control(portIndex)
		{
			host->vctrl[host->ctrl_count++] = this;
		}
		virtual ~tcontrol() {}

		void connect(void* data) override
		{
			p = static_cast<const ptype*>(data);
		}
		void commit() override
		{
			if (p != nullptr) {
				if (v != *p) {
					vnew = true;
					v = *p;
				}
			}
		}

		operator const ptype() const
		{
			return v;
		}
		ptype get() const
		{
			return v;
		}

	private:
		const ptype* p = nullptr;
		ptype v = 0;
	};
	typedef tcontrol<float> fcontrol;
	typedef tcontrol<int> icontrol;

	class input
	{
	public:
		input(plugin* host, int portIndex);
		~input();
		int getIndex() const;

		void connect(void* data);
		operator const float* () const;
		const float* get() const;

	private:
		const int index = 0;
		const float* p = nullptr;
	};

	class output
	{
	public:
		output(plugin* host, int portIndex);
		~output();
		int getIndex() const;

		void connect(void* data);
		operator float* ();
		float* get();

	private:
		const int index = 0;
		float* p = nullptr;
	};

protected:
	float Fs = 0.0f;

private:
	int ctrl_count = 0;
	int in_count = 0;
	int out_count = 0;

	vec<control*> vctrl;
	vec<input*> vin;
	vec<output*> vout;
};

//==================================================

//Returns 1, 2, or 4 to represet 1x, 2x, and 4x nominal
//sampling rates. Rates < 48k yield 1, and > 192 yield 4,
//otherwise the nearest rate is used:
//	1 - 48k
//	2 - 96k
//	4 - 192k
int getNominalUSR(float Fs);


//==================================================

struct Constants
{
	const float pi = 3.1415926535897932384626433832795f;
	const float e = 2.7182818284590452353602874713527f;
	const float eps = 1e-37f; //theoretical eps is around 1e-38
	const float p707 = 0.70710678118654752440084436210485f;

	const float db2nats = 0.11512925464970228420089957273422f;
	const float common2nats = 2.3025850929940456840179914546844f;
	const float two2nats = 0.69314718055994530941723212145818f;
};
//Collection of mathematical constants
extern const Constants constants;

//==================================================

template<class ptype> inline ptype max(ptype x, ptype y)
{
	return (x > y) ? x : y;
}
template<class ptype> inline ptype min(ptype x, ptype y)
{
	return (x < y) ? x : y;
}
template<class ptype> inline ptype bound(ptype x, ptype val1, ptype val2)
{
	if (val1 <= val2) {
		return (x > val1) ? ((x < val2) ? x : val2) : val1;
	}
	else {
		return (x > val2) ? ((x < val1) ? x : val1) : val2;
	}
}

template<class ptype> void vmax(const ptype* x, ptype y, ptype* z, int len)
{
	for (int i = 0; i < len; ++i) {
		z[i] = max(x[i], y);
	}
}
template<class ptype> void vmax(const ptype* x, const ptype* y, ptype* z, int len)
{
	for (int i = 0; i < len; ++i) {
		z[i] = max(x[i], y[i]);
	}
}

template<class ptype> void vmin(const ptype* x, ptype y, ptype* z, int len)
{
	for (int i = 0; i < len; ++i) {
		z[i] = min(x[i], y);
	}
}
template<class ptype> void vmin(const ptype* x, const ptype* y, ptype* z, int len)
{
	for (int i = 0; i < len; ++i) {
		z[i] = min(x[i], y[i]);
	}
}

template<class ptype> void vbound(const ptype* x, ptype val1, ptype val2, ptype* z, int len)
{
	for (int i = 0; i < len; ++i) {
		z[i] = bound(x[i], val1, val2);
	}
}

template<class ptype> ptype vfindmax(const ptype* x, int len)
{
	ptype ret = x[0];
	for (int i = 1; i < len; ++i) {
		ret = max(ret, x[i]);
	}
	return ret;
}
template<class ptype> ptype vfindmax(const ptype* x, int len, int& index)
{
	ptype ret = x[0];
	index = 0;
	for (int i = 1; i < len; ++i) {
		if (x[i] > ret) {
			ret = x[i];
			index = i;
		}
	}
	return ret;
}
template<class ptype> ptype vfindmin(const ptype* x, int len)
{
	ptype ret = x[0];
	for (int i = 1; i < len; ++i) {
		ret = min(ret, x[i]);
	}
	return ret;
}
template<class ptype> ptype vfindmin(const ptype* x, int len, int& index)
{
	ptype ret = x[0];
	index = 0;
	for (int i = 1; i < len; ++i) {
		if (x[i] < ret) {
			ret = x[i];
			index = i;
		}
	}
	return ret;
}
//Notes: y == x is safe, y must be of length len even if N < len, this is a slow sort
template<class ptype> void vsortmin(const ptype* x, ptype* y, int len, int N=-1)
{
	if (x != y) {
		for (int i = 0; i < len; ++i) {
			y[i] = x[i];
		}
	}
	if (N < 0) {
		N = len;
	}
	ptype curmin;
	int curind = 0;
	for (int n = 0; n < N; ++n) {
		curmin = vfindmin(y, len, curind);
		y[curind] = y[0]; //swap
		*y++ = curmin; //step
		--len;
	} //returns with the first N values of y holding the ascending minimal values
}
//Notes: y == x is safe, y and ind must be of length len even if N < len, this is a slow sort
template<class ptype> void vsortmin(const ptype* x, ptype* y, int* ind, int len, int N=-1)
{
	if (x != y) {
		for (int i = 0; i < len; ++i) {
			y[i] = x[i];
		}
	}
	if (N < 0) {
		N = len;
	}
	for (int i = 0; i < len; ++i) {
		ind[i] = i;
	}
	ptype curmin;
	int curind = 0;
	int nextind = 0;
	for (int n = 0; n < N; ++n) {
		curmin = vfindmin(y, len, curind);
		nextind = ind[curind];
		y[curind] = y[0]; //swap
		ind[curind] = ind[0];
		*y++ = curmin; //step
		*ind++ = nextind;
		--len;
	} //returns with the first N values of y/ind holding ascending minimal values and indices
}

//Wraps the input vector of 16-bit integer indices into y given the power-of-2
//wrapping length len = 2^r.
void vwrap(const uint16_t* x, uint16_t* y, unsigned int r);

//Converts a sorted index vector ind of length len into a place vector
//of the same length. Ex. [2, 0, 1] -> [1, 2, 0]. ind == place is UNSAFE
void unsort(const int* ind, int* place, int len);

inline float linmap(float x, float xmin, float xmax, float ymin, float ymax)
{
	return (x - xmin) * (ymax - ymin) / (xmax - xmin) + ymin;
}

float boundlinmap(float x, float xmin, float xmax, float ymin, float ymax);

float boundlogmap(float x, float xmin, float xmax, float ymin, float ymax);

//==================================================

template<class ptype> inline ptype sign(ptype a) {
	if (a < (ptype)0) {
		return (ptype)-1;
	}
	else if (a > (ptype)0) {
		return (ptype)1;
	}
	else {
		return (ptype)0;
	}
}

inline float fix(float a) {
	if (a >= 0.0f) {
		return floorf(a);
	}
	else {
		return ceilf(a);
	}
}
inline double fix(double a) {
	if (a >= 0.0) {
		return floor(a);
	}
	else {
		return ceil(a);
	}
}

//==================================================

template<class xtype, class ytype> void vset(xtype* x, ytype val, int len)
{
	for (int i = 0; i < len; ++i) {
		x[i] = val;
	}
}

template<class xtype, class ytype> void vincspace(xtype* x, ytype val0, ytype dval, int len)
{
	for (int i = 0; i < len; ++i) {
		x[i] = val0 + i * dval;
	}
}

template<class xtype, class ytype> void vcopy(const xtype* x, ytype* y, int len)
{
	for (int i = 0; i < len; ++i) {
		y[i] = (const ytype)x[i];
	}
}

//x == y is UNSAFE!
template<class xtype, class ytype> void vrevcopy(const xtype* x, ytype* y, int len)
{
	int j = len;
	for (int i = 0; i < len; ++i) {
		y[--j] = (const ytype)x[i];
	}
}

template<class xtype, class ytype> void vpadcopy(const xtype* x, ytype* y, ytype padval, int xlen, int ylen)
{
	xlen = min(xlen, ylen);
	for (int i = 0; i < xlen; ++i) {
		y[i] = (const ytype)x[i];
	}
	for (int i = xlen; i < ylen; ++i) {
		y[i] = padval;
	}
}

template<class xtype, class ytype> void vcirccopy(const xtype* x, ytype* y, int xlen, int xfirst, int ylen)
{
	xfirst %= xlen; //map the shift to a value between 0 and xlen
	if (xfirst < 0) {
		xfirst += xlen;
	}

	int curlen = 0;
	while (ylen > 0) {
		curlen = min(ylen, xlen - xfirst);
		vcopy(x + xfirst, y, curlen);
		ylen -= curlen;
		y += curlen;
		xfirst = 0;
	}
}

//==================================================

template<class ptype> ptype sum(const ptype* x, int len)
{
	ptype ret = (ptype)0;
	for (int i = 0; i < len; ++i) {
		ret += x[i];
	}
	return ret;
}
template<class ptype> ptype mean(const ptype* x, int len) {
	return sum(x, len) / len;
}
template<class ptype> ptype ssq(const ptype* x, int len) {
	ptype ret = (ptype)0;
	for (int i = 0; i < len; ++i) {
		ret += x[i] * x[i];
	}
	return ret;
}
float rms(const float* x, int len);
double rms(const double* x, int len);
float rssq(const float* x, int len);
double rssq(const double* x, int len);

template<class ptype> void vabs(const ptype* x, ptype* y, int len)
{
	for (int i = 0; i < len; ++i) {
		y[i] = abs(x[i]);
	}
}

template<class ptype> void vabs(const complex<ptype>* x, ptype* y, int len)
{
	for (int i = 0; i < len; ++i) {
		y[i] = abs(x[i]);
	}
}

template<class ptype> void vpow(const ptype* x, ptype* y, int len)
{
	for (int i = 0; i < len; ++i) {
		y[i] = x[i] * x[i];
	}
}

template<class ptype> void vpow(const complex<ptype>* x, ptype* y, int len)
{
	for (int i = 0; i < len; ++i) {
		y[i] = norm(x[i]);
	}
}

template<class ptype> void vconj(const complex<ptype>* x, complex<ptype>* y, int len)
{
	for (int i = 0; i < len; ++i) {
		y[i] = conj(x[i]);
	}
}

template<class ptype> ptype vdot(const ptype* x, const ptype* y, int len)
{
	ptype ret = 0.0f;
	for (int i = 0; i < len; ++i) {
		ret += x[i] * y[i];
	}
	return ret;
}

//==================================================

template<class ptype> void vmult(const ptype* x, ptype scalar, ptype* y, int len)
{
	if (x == y) {
		for (int i = 0; i < len; ++i) {
			y[i] *= scalar;
		}
	}
	else {
		for (int i = 0; i < len; ++i) {
			y[i] = x[i] * scalar;
		}
	}
}
template<class ptype> void vmult(const ptype* x1, const ptype* x2, ptype* y, int len)
{
	if (x1 == y) {
		if (x2 == y) {
			for (int i = 0; i < len; ++i) {
				y[i] *= y[i];
			}
		}
		else {
			for (int i = 0; i < len; ++i) {
				y[i] *= x2[i];
			}
		}
	}
	else if (x2 == y) {
		for (int i = 0; i < len; ++i) {
			y[i] *= x1[i];
		}
	}
	else {
		for (int i = 0; i < len; ++i) {
			y[i] = x1[i] * x2[i];
		}
	}
}
template<class ptype> void vadd(const ptype* x, ptype scalar, ptype* y, int len)
{
	if (x == y) {
		for (int i = 0; i < len; ++i) {
			y[i] += scalar;
		}
	}
	else {
		for (int i = 0; i < len; ++i) {
			y[i] = x[i] + scalar;
		}
	}
}
template<class ptype> void vadd(const ptype* x1, const ptype* x2, ptype* y, int len)
{
	if (x1 == y) {
		if (x2 == y) {
			for (int i = 0; i < len; ++i) {
				y[i] *= (ptype)2; //special case
			}
		}
		else {
			for (int i = 0; i < len; ++i) {
				y[i] += x2[i];
			}
		}
	}
	else if (x2 == y) {
		for (int i = 0; i < len; ++i) {
			y[i] += x1[i];
		}
	}
	else {
		for (int i = 0; i < len; ++i) {
			y[i] = x1[i] + x2[i];
		}
	}
}
template<class ptype> void vdiv(const ptype* x, ptype scalar, ptype* y, int len)
{
	if (x == y) {
		for (int i = 0; i < len; ++i) {
			y[i] /= scalar;
		}
	}
	else {
		for (int i = 0; i < len; ++i) {
			y[i] = x[i] / scalar;
		}
	}
}
template<class ptype> void vdiv(ptype scalar, const ptype* x, ptype* y, int len)
{
	if (x == y) {
		for (int i = 0; i < len; ++i) {
			y[i] = scalar / y[i];
		}
	}
	else {
		for (int i = 0; i < len; ++i) {
			y[i] = scalar / x[i];
		}
	}
}
template<class ptype> void vdiv(const ptype* x1, const ptype* x2, ptype* y, int len)
{
	if (x1 == y) {
		if (x2 == y) {
			for (int i = 0; i < len; ++i) {
				y[i] = (ptype)1; //special case
			}
		}
		else {
			for (int i = 0; i < len; ++i) {
				y[i] /= x2[i];
			}
		}
	}
	else if (x2 == y) {
		for (int i = 0; i < len; ++i) {
			y[i] = x1[i] / y[i];
		}
	}
	else {
		for (int i = 0; i < len; ++i) {
			y[i] = x1[i] / x2[i];
		}
	}
}
template<class ptype> void vsub(const ptype* x, ptype scalar, ptype* y, int len)
{
	if (x == y) {
		for (int i = 0; i < len; ++i) {
			y[i] -= scalar;
		}
	}
	else {
		for (int i = 0; i < len; ++i) {
			y[i] = x[i] - scalar;
		}
	}
}
template<class ptype> void vsub(ptype scalar, const ptype* x, ptype* y, int len)
{
	if (x == y) {
		for (int i = 0; i < len; ++i) {
			y[i] = scalar - y[i];
		}
	}
	else {
		for (int i = 0; i < len; ++i) {
			y[i] = scalar - x[i];
		}
	}
}
template<class ptype> void vsub(const ptype* x1, const ptype* x2, ptype* y, int len)
{
	if (x1 == y) {
		if (x2 == y) {
			for (int i = 0; i < len; ++i) {
				y[i] = (ptype)0; //special case
			}
		}
		else {
			for (int i = 0; i < len; ++i) {
				y[i] -= x2[i];
			}
		}
	}
	else if (x2 == y) {
		for (int i = 0; i < len; ++i) {
			y[i] = x1[i] - y[i];
		}
	}
	else {
		for (int i = 0; i < len; ++i) {
			y[i] = x1[i] - x2[i];
		}
	}
}

template<class ptype> void vmultaccum(const ptype* x, ptype scalar, ptype* y, int len)
{
	if (x == y) {
		scalar += 1.0f;
		for (int i = 0; i < len; ++i) {
			y[i] *= scalar;
		}
	}
	else {
		for (int i = 0; i < len; ++i) {
			y[i] += x[i] * scalar;
		}
	}
}
template<class ptype> void vmultaccum(const ptype* x1, const ptype* x2, ptype* y, int len)
{
	for (int i = 0; i < len; ++i) {
		y[i] += x1[i] * x2[i];
	}
}

//==================================================

void vmagdB(const float* x, float* y, int len);
void vpowdB(const float* x, float* y, int len);

void vlog(const float* x, float* y, int len);
void vexp(const float* x, float* y, int len);
void vsqrt(const float* x, float* y, int len);

void vexpi(const float* ph, complex<float>* y, int len);

//==================================================

//Solves for the fractional x value of a parabola that intersects
//all three points (-1, ym), (0, y0), (1, yp) when y0 is a local maxima
inline float parabsolve(float ym, float y0, float yp)
{
	return 0.5f * (ym - yp) / (ym + yp - 2.0f * (y0 + constants.eps));
}

//==================================================

//Replaces NaN values with 0 and saturates values >40 dB amplitude
//The input is const, so you MUST copy inputs to scratch blocks
void safeInput(const float* pin, float* x, int len);

//Replaces NaN values with 0 and saturates values >40 dB amplitude,
//also returns true if any NaN or value >80 dB amplitude is found
//so that the module can clear internal state
bool safeOutput(float* pout, int len);

//==================================================

//Linear SOE solver class: EQ*x = SOL
class solver
{
public:
	solver(int setN, double setTolerance = 0.0);
	~solver();

	//Returns x such that EQ*x = SOL using LU decomposition with pivoting
	//performed immediately in full
	bool solve(const float* const* EQ, const float* SOL, float* x);

	//Prepares internal state to quickly solve problems of the form EQ*x = SOL
	//with the resolve() function
	bool prepare(const float* const* EQ);

	//Returns x such that EQ*x = SOL using whatever the last set EQ was by
	//either solve() or prepare()
	void resolve(const float* SOL, float* x);

	//Evaluates EQ*x = SOL, returning SOL to check correctness
	void check(const float* const* EQ, const float* x, float* SOL) const;

private:
	const int N = 1;
	const double tol = 0.0;
	double** A = nullptr;
	double* b = nullptr;
	double* v = nullptr;
	int* P = nullptr;
};

//==================================================

//Object that breaks up incoming audio blocks into chunks <= a given
//maxlen while enforcing that no rebuffered chunk crosses a maxlen
//buffer boundary. This is helpful for guaranteeing performance of 
//some block-based algorithms.
class rebuffer
{
public:
	rebuffer(int maxlen);
	~rebuffer();

	//Clears internal sample count to 0. Typically called on init or
	//some other whole plugin reinitialization.
	void reset();

	//Returns the next application block size when remlen samples
	//remain in the I/O buffers. Typical usage:
	//	//pointer p holds data, len holds block length
	// while(len > 0){
	//	curlen = rebuffer.next(len);
	//	<DO PROCESSING ON p OVER curlen SAMPLES>
	//	p += curlen; //step the data pointer forward
	//  len -= curlen; //reduce the remaining samples
	// }
	int next(int remlen);

private:
	const int N = 0;
	int n = 0;
};


//==================================================

//Pure virtual class that defines common functions for mono
//DSP algorithms. Mostly helpful for automatically generating
//block processing functions from sample processing ones.
class monoalg
{
public:
	monoalg();
	virtual ~monoalg();

	virtual float step(float x) = 0;
	virtual void stepBlock(const float* x, float* y, int len);
};

//Pure virtual class that defines common functions for stereo
//DSP algorithms. Mostly helpful for automatically generating
//block processing functions from sample processing ones.
class stereoalg
{
public:
	stereoalg();
	virtual ~stereoalg();

	virtual void step(float xL, float xR, float& yL, float& yR) = 0;
	virtual void stepBlock(const float* xL, const float* xR,
		float* yL, float* yR, int len);
};

//Pure virtual class that defines common functions for multichannel
//DSP algorithms. Mostly helpful for automatically generating
//block processing functions from sample processing ones.
class multialg
{
public:
	//if blockTransposed = true, stepBlock calls use data blocks of size len by chan
	//instead of chan by length -- this saves some copy overhead
	multialg(int maxNch, bool blockTransposed = false);
	virtual ~multialg();

	virtual void reallocate(int maxNch, bool blockTransposed = false);

	virtual void step(const float* x, float* y, int chan) = 0;
	virtual void stepBlock(const float*const* x, float*const* y, int chan, int len);

protected:
	bool transposed = false;
	fvec xscr, yscr;
};


//==================================================

struct range
{
public:
	friend class multirange;

	enum class shape{
		LIN = 0,
		LOG
	};
	range(float setxmin = 0.0f, float setxmax = 1.0f,
		float setymin = 0.0f, float setymax = 1.0f,
		shape setshape = shape::LIN);
	~range();

	void set(float setxmin, float setxmax,
		float setymin, float setymax,
		shape setshape = shape::LIN);

	float map(float x) const;
	float unmap(float y) const;

private:
	float xmin = 0.0f;
	float xmax = 1.0f;
	float ymin = 0.0f;
	float ymax = 1.0f;
	shape rshape = shape::LIN;
};


//==================================================

class multirange
{
public:
	typedef range::shape shape;

	multirange(int stops);
	~multirange();

	//the shape applies from this x to the next one, unused at the list end
	void setStop(int index, float x, float y, shape setshape = shape::LIN);

	float map(float x) const;
	//no unmap, implement if necessary... but that's hard

private:
	fvec xstop;
	fvec ystop;
	vec<range> r;
};

//==================================================

//Precomputes values of a slow math function like log(), cos(), and exp()
//and stores them in a lookup table for quick access in cases where
//performance is more important than accuracy.
template<class xtype, class ytype> class lutop
{
public:
	lutop(xtype setxmin, xtype setxmax, int N) :
		xmin(setxmin), xmax(setxmax), xdelta((setxmax - setxmin) / N),
		lut(N)
	{}
	virtual ~lutop() {}

	inline xtype bnd(xtype x) const
	{
		return bound(x, xmin, xmax);
	}

	xtype* bnd(const xtype* x, xtype* y, int len) const
	{
		vbound(x, xmin, xmax, y, len);
		return y;
	}

	//Call bnd(x) first if it is unknown whether it will be within range
	virtual ytype operator()(xtype x) const
	{
		int ind = (int)std::round((x - xmin) / xdelta);
		return lut[ind];
	}

	virtual void operator()(const xtype* x, ytype* y, int len) const
	{
		int ind = 0;
		for (int i = 0; i < len; ++i) {
			ind = (int)std::round((x[i] - xmin) / xdelta);
			y[i] = lut[ind];
		}
	}

protected:
	const xtype xmin;
	const xtype xmax;
	const xtype xdelta;
	vec<ytype> lut;
};

class loglut : public lutop<float, float>
{
public:
	loglut(float setxmin, float setxmax, int N) : lutop(setxmin, setxmax, N)
	{
		float x = xmin;
		for (int i = 0; i < N; ++i) {
			lut[i] = logf(x);
			x += xdelta;
		}
	}
	~loglut() {}
};

//==================================================

//Smoothed scalar value object that uses a first order all-pole
//filter to update its value on every sample.
template<class ptype> class smoothed
{
public:
	smoothed(int samples = 1, ptype initval = (ptype)0)
	{
		setSmooth(samples);
		target(initval, true);
	}
	~smoothed() {}

	void setSmooth(int samples)
	{
		salpha = expf(-1.0f / samples);
		sbeta = 1.0f - salpha;
	}

	void target(ptype set, bool convergeInstantly = false)
	{
		next = (double)set;
		if (convergeInstantly) {
			converge();
		}
	}

	void smooth(int len = 1)
	{
		mem = salpha * mem + sbeta * next;
		val = (ptype)mem;
	}
	//smooths a block sample by sample, populating a value trajectory
	void smoothBlock(ptype* v, int len)
	{
		for (int n = 0; n < len; ++n) {
			mem = salpha * mem + sbeta * next;
			v[n] = (ptype)mem;
		}
		val = v[len - 1];
	}

	operator ptype() const
	{
		return val;
	}
	ptype get() const
	{
		return val;
	}

	void converge()
	{
		mem = next;
		val = (ptype)mem;
	}

private:
	double salpha = 0.0;
	double sbeta = 1.0;
	double mem = 0.0;
	double next = 0.0;

	ptype val = (ptype)0;
};

//Slewed scalar value object that allows retargetting, checking,
//updating, and retrieval.
template<class ptype> class slewed
{
public:
	slewed(int samples = 1, ptype initval = (ptype)0) :
		L(max(samples, 1)), delta(1.0 / (double)max(samples, 1)), norm(1.0),
		val(initval), last(initval), next(initval)
	{}
	~slewed() {}

	void setSlew(int samples)
	{
		L = max(samples, 1);
		delta = 1.0 / L;
	}
	const int getSlew() const
	{
		return L;
	}

	void target(ptype set, bool convergeInstantly = false)
	{
		if (val == set) {
			norm = 1.0;
			last = val;
			next = set;
		}
		else if (next != set) {
			norm = 0.0;
			last = val;
			next = set;
		} //otherwise, if next == set leave the trajectory alone!

		if (convergeInstantly) {
			converge();
		}
	}
	void pause()
	{
		active = false;
	}
	void unpause()
	{
		active = true;
	}
	bool check() const
	{
		return (norm < 1.0) && active ; //true if currently slewing
	}
	bool slew(int len = 1)
	{
		if (check()) {
			norm = min(norm + len * delta, 1.0);
			val = (ptype)(norm * next + (1.0 - norm) * last);
			return true; //new value
		}
		else {
			return false; //no new value
		}
	}
	//slews a block sample by sample, populating a value trajectory
	void slewBlock(ptype* v, int len)
	{
		double temp = 0.0;
		for (int i = 0; i < len; ++i) {
			temp = min(norm + i * delta, 1.0);
			v[i] = (ptype)(temp * next + (1.0 - temp) * last);
		}
		norm = temp;
		val = v[len - 1];
	}

	operator ptype() const
	{
		return val;
	}
	ptype get() const
	{
		return val;
	}

	bool converge()
	{
		norm = 1.0;
		if (val != next) {
			val = next;
			return true; //new value
		}
		else {
			return false; //no new value
		}
	}

private:
	int L = 1; //slew in samples
	double delta = 1.0; //norm step per sample
	double norm = 1.0; //current normalized value on 0 to 1
	bool active = true; //whether slewing is active or paused

	ptype val = (ptype)0; //output value
	ptype last = (ptype)0; //previous value
	ptype next = (ptype)0; //next value
};

//Slewed vector value object that allows retargetting, checking,
//updating, and retrieval over a vector of coupled parameters
template<class ptype> class slewedvec
{
public:
	slewedvec(int len, int samples = 1) : N(len),
		L(max(samples, 1)), delta(1.0 / (double)max(samples, 1)), norm(1.0)
	{
		if (N > 0) {
			pval = new ptype[N];
			plast = new ptype[N];
			pnext = new ptype[N];
		}
	}
	~slewedvec()
	{
		if (N > 0) {
			delete[] pval;
			delete[] plast;
			delete[] pnext;
		}
	}

	int size()
	{
		return N;
	}

	void setSlew(int samples)
	{
		L = max(samples, 1);
		delta = 1.0 / L;
	}

	ptype* target()
	{
		for (int n = 0; n < N; ++n) {
			plast[n] = pval[n];
		}
		norm = 0.0;
		return pnext;
	}
	void pause()
	{
		active = false;
	}
	void unpause()
	{
		active = true;
	}
	bool check() const
	{
		return (norm < 1.0) && active; //true if currently slewing
	}
	bool slew(int len = 1)
	{
		if (check()) { //TODO <-- someday it would be rad to get this better optimized!
			norm = min(norm + len * delta, 1.0);
			const double bnorm = 1.0 - norm;
			for (int n = 0; n < N; ++n) {
				pval[n] = (ptype)(norm * pnext[n] + bnorm * plast[n]);
			}
			return true; //new values
		}
		else {
			return false; //no new values
		}
	}

	operator const ptype*() const
	{
		return pval;
	}
	const ptype* get() const
	{
		return pval;
	}

	bool converge()
	{
		norm = 1.0;
		for (int n = 0; n < N; ++n) {
			pval[n] = pnext[n];
		}
		return true; //no checks on this one, return true to keep parallelism with scalar
	}

private:
	int N = 1; //vector length
	int L = 1; //slew in samples
	double delta = 1.0; //norm step per sample
	double norm = 1.0; //current normalized value on 0 to 
	bool active = true; //whether slewing is active or paused

	ptype* pval = nullptr; //output values
	ptype* plast = nullptr; //previous values
	ptype* pnext = nullptr; //next values
};

//Scalar value that permits crossfading between value changes,
//similar to a slewed value, except both new and old values are
//accessible.
template<class ptype> class faded
{
public:
	faded(int samples = 1, ptype initval = (ptype)0) :
		L(max(samples, 1)), delta(1.0f / (float)max(samples, 1)), g(1.0f),
		curval(initval), lastval(initval), nextval(initval)
	{}
	~faded() {}

	void setFade(int samples)
	{
		L = max(samples, 1);
		delta = 1.0f / L;
	}

	void target(ptype set, bool convergeInstantly = false)
	{
		nextval = set; //save the value for the next fade
		if (convergeInstantly) {
			converge();
		}
	}
	//permit value swapping, starting a new crossfade,
	//returns true if a new crossfade has started
	bool swap()
	{
		bool ret = false;
		if (nextval != curval) {
			if (g == 1.0f) {
				lastval = curval;
				curval = nextval;
				g = 0.0f;
				ret = true;
			}
		}
		return ret;
	}
	//true if currently cross-fading
	bool check() const
	{
		return (g < 1.0f); //true if currently fading
	}
	//returns true if currently cross-fading
	bool fade(int len = 1)
	{
		if (g < 1.0f) {
			g = min(g + len * delta, 1.0f);
			return true; //new crossfade
		}
		else {
			return false; //no new crossfade
		}
	}

	ptype current() const
	{
		return curval;
	}
	ptype last() const
	{
		return lastval;
	}
	//not involved in active fade! sometimes helpful to see the value that is "on deck"
	ptype next() const
	{
		return nextval;
	}
	float gain() const
	{
		return g;
	}

	//fades a block sample by sample, populating a gain trajectory
	void fadeBlock(float* gcur, int len)
	{
		vincspace(gcur, g, delta, len);
		vmin(gcur, 1.0f, gcur, len);
		g = gcur[len - 1];
	}

	bool converge()
	{
		bool ret = false;
		if (nextval != curval) {
			lastval = curval;
			curval = nextval;
			ret = true;
		}
		g = 1.0f;
		return ret;
	}

private:
	int L = 1; //slew in samples
	float delta = 1.0; //g step per sample
	float g = 1.0; //current current value crossfade gain

	ptype curval = (ptype)0; //current value
	ptype lastval = (ptype)0; //last value
	ptype nextval = (ptype)0; //next (on deck) value
};

//Vector value that permits crossfading between value changes,
//similar to a slewed vector, except both new and old values are
//accessible.
template<class ptype> class fadedvec
{
public:
	fadedvec(int len, int samples = 1) : N(len),
		L(max(samples, 1)), delta(1.0f / (float)max(samples, 1)), g(1.0f),
		curval(len), lastval(len), nextval(len)
	{}
	~fadedvec() {}

	void setFade(int samples)
	{
		L = max(samples, 1);
		delta = 1.0f / L;
	}

	ptype* target()
	{
		hasNextVal = true;
		return nextval.ptr();
	}
	//permit value swapping, starting a new crossfade,
	//returns true if a new crossfade has started
	bool swap()
	{
		bool ret = false;
		if (hasNextVal) {
			if (g == 1.0f) {
				vcopy(curval.ptr(), lastval.ptr(), N);
				vcopy(nextval.ptr(), curval.ptr(), N);
				hasNextVal = false;
				g = 0.0f;
				ret = true;
			}
		}
		return ret;
	}
	//true if currently cross-fading
	bool check() const
	{
		return (g < 1.0f); //true if currently fading
	}
	//returns true if currently cross-fading
	bool fade(int len = 1)
	{
		if (g < 1.0f) {
			g = min(g + len * delta, 1.0f);
			return true; //new crossfade
		}
		else {
			return false; //no new crossfade
		}
	}

	const ptype* current() const
	{
		return curval.ptr();
	}
	const ptype* last() const
	{
		return lastval.ptr();
	}
	float gain() const
	{
		return g;
	}

	//fades a block sample by sample, populating a gain trajectory
	void fadeBlock(float* gcur, int len)
	{
		vincspace(gcur, g, delta, len);
		vmin(gcur, 1.0f, gcur, len);
		g = gcur[len - 1];
	}

	bool converge()
	{
		bool ret = false;
		if (hasNextVal) {
			vcopy(curval.ptr(), lastval.ptr(), N);
			vcopy(nextval.ptr(), curval.ptr(), N);
			hasNextVal = false;
			ret = true;
		}
		g = 1.0f;
		return ret;
	}

private:
	int N = 1; //vector length
	int L = 1; //slew in samples
	float delta = 1.0; //g step per sample
	float g = 1.0; //current current value crossfade gain

	vec<ptype> curval; //current value
	vec<ptype> lastval; //last value
	vec<ptype> nextval; //next (on deck) value
	bool hasNextVal = false; //flags when nextval may be different from curval
};


}