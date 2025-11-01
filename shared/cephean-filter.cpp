
//--------------------------------------------------
// Defines first, second, and higher order filters
//--------------------------------------------------

#include "cephean-filter.h"

namespace cephean
{

//==================================================

fof::fof(int slew, int maxNch, bool blockTransposed) :
	multialg(maxNch, blockTransposed),
	sc(3, slew), vmem(maxNch)
{
	setCoefs(coefs{ 0.0, 1.0, 0.0 }, true); //init to bypass filter
	clear();
}
fof::~fof() {}

void fof::reallocate(int slew, int maxNch, bool blockTransposed)
{
	multialg::reallocate(maxNch, blockTransposed);
	sc.setSlew(slew);
	vmem.reset(maxNch);
	clear();
}

void fof::setCoefs(coefs newCoefs, bool converge)
{
	double* pc = sc.target();
	c = newCoefs; //save a reference copy
	pc[0] = newCoefs.na1;
	pc[1] = newCoefs.b0;
	pc[2] = newCoefs.b1;
	if (converge) {
		sc.converge();
	}
}

fof::coefs fof::getCoefs() const
{
	return c;
}

void fof::clear()
{
	memL = 0.0; //mono/stereo
	memR = 0.0;
	vset(vmem.ptr(), 0.0, vmem.size()); //multichannel
}

float fof::step(float x)
{
	sc.slew();
	const double* ct = sc.get();
	double temp = x + ct[0] * memL;
	x = (float)(ct[1] * temp + ct[2] * memL);
	memL = temp;
	return x;
}
void fof::stepBlock(const float* x, float* y, int len)
{
	if (sc.check()) { //slew every sample
		monoalg::stepBlock(x, y, len);
	}
	else { //run faster fixed coefficient alg
		double temp;
		const double* ct = sc.get();
		for (int i = 0; i < len; i++) {
			temp = x[i] + ct[0] * memL;
			y[i] = (float)(ct[1] * temp + ct[2] * memL);
			memL = temp;
		}
	}
}

void fof::step(float xL, float xR, float& yL, float& yR)
{
	sc.slew();
	const double* ct = sc.get();
	double tempL = xL + ct[0] * memL;
	double tempR = xR + ct[0] * memR;
	yL = (float)(ct[1] * tempL + ct[2] * memL);
	yR = (float)(ct[1] * tempR + ct[2] * memR);
	memL = tempL;
	memR = tempR;
}
void fof::stepBlock(const float* xL, const float* xR, float* yL, float* yR, int len)
{
	if (sc.check()) { //slew every sample
		stereoalg::stepBlock(xL, xR, yL, yR, len);
	}
	else { //run faster fixed coefficient alg
		double tempL, tempR;
		const double* ct = sc.get();
		for (int i = 0; i < len; i++) {
			tempL = xL[i] + ct[0] * memL;
			tempR = xR[i] + ct[0] * memR;
			yL[i] = (float)(ct[1] * tempL + ct[2] * memL);
			yR[i] = (float)(ct[1] * tempR + ct[2] * memR);
			memL = tempL;
			memR = tempR;
		}
	}
}

void fof::step(const float* x, float* y, int chan)
{
	sc.slew();
	const double* ct = sc.get();
	double temp = 0.0;
	for (int i = 0; i < chan; ++i) {
		temp = x[i] + ct[0] * vmem[i];
		y[i] = (float)(ct[1] * temp + ct[2] * vmem[i]);
		vmem[i] = temp;
	}
}
void fof::stepBlock(const float* const* x, float* const* y, int chan, int len)
{
	if (sc.check()) {
		multialg::stepBlock(x, y, chan, len);
	}
	else { //run faster fixed coefficient alg
		const double* ct = sc.get();
		double temp = 0.0;
		const float* px = nullptr;
		float* py = nullptr;
		if (transposed) {
			for (int i = 0; i < len; ++i) {
				px = x[i];
				py = y[i];
				for (int j = 0; j < chan; ++j) {
					temp = px[j] + ct[0] * vmem[j];
					py[j] = (float)(ct[1] * temp + ct[2] * vmem[j]);
					vmem[j] = temp;
				}
			}
		}
		else {
			double tmem = 0.0;
			for (int i = 0; i < chan; ++i) {
				px = x[i];
				py = y[i];
				tmem = vmem[i];
				for (int j = 0; j < len; ++j) {
					temp = px[j] + ct[0] * tmem;
					py[j] = (float)(ct[1] * temp + ct[2] * tmem);
					tmem = temp;
				}
				vmem[i] = tmem;
			}
		}
	}
}

namespace filt
{

// First order filter definitions
//====================================================

//Inverts the phase of the passed in coefficients structure
fof::coefs invert(const fof::coefs& coef)
{
	return { coef.na1, -coef.b0, -coef.b1 };
}

//Returns the filter inverse of a stable minimum phase coefs structure
fof::coefs inverse(const fof::coefs& coef)
{
	return { -coef.b1 / coef.b0, 1.0 / coef.b0, -coef.na1 / coef.b0 };
}

//Returns an allpass filter that phase matches a double application of
//the provided filter with one stage's numerator coefficients flipped
fof::coefs matchphase(const fof::coefs& coef)
{
	return { coef.na1, -coef.na1, 1.0 };
}

//Flips the numerator coefficients of the provided filter, switching
//min phase <-> max phase, usually for phase aligned double filtering
fof::coefs mp2mp(const fof::coefs& coef)
{
	return { coef.na1, coef.b1, coef.b0 };
}

//====================================================

//First order unity filter, returns default
fof::coefs unity1()
{
	return { 0.0, 1.0, 0.0 };
}

//First order trivial gain filter with gain multiplier g
fof::coefs gain1(float g)
{
	return { 0.0, (double)g, 0.0 };
}

//====================================================

//First order lowpass with rolloff frequency f
fof::coefs lowpass1(float f)
{
	fof::coefs coef;
	double w = (double)tan(f * constants.pi / 2.0f);
	coef.na1 = (1.0 - w) / (w + 1.0);
	coef.b0 = 0.5 * (1.0 - coef.na1);
	coef.b1 = coef.b0;
	return coef;
}

//First order compensated lowpass with rolloff frequency f and
//bounds minf, maxf for stability
fof::coefs lowpass1(float f, float minf, float maxf)
{
	fof::coefs coef;
	double w = (double)tan(bound(f, minf, maxf) * constants.pi / 2.0f);
	double glo = (double)min(f / minf, 1.0f);
	double ghi = (double)min(f / maxf, 1.0f);
	coef.na1 = (1.0 - w) / (w + 1.0);
	coef.b0 = 0.5 * ((glo + ghi) - coef.na1 * (glo - ghi));
	coef.b1 = 0.5 * ((glo - ghi) - coef.na1 * (glo + ghi));
	return coef;
}

//First order highpass with rolloff frequency f
fof::coefs highpass1(float f)
{
	fof::coefs coef;
	double w = (double)tan(f * constants.pi / 2.0f);
	coef.na1 = (1.0 - w) / (w + 1.0);
	coef.b0 = 0.5 * (1.0 + coef.na1);
	coef.b1 = -coef.b0;
	return coef;
}

//First order compensated highpass with rolloff frequency f and
//bounds minf, maxf for stability
fof::coefs highpass1(float f, float minf, float maxf)
{
	fof::coefs coef;
	double w = (double)tan(bound(f, minf, maxf) * constants.pi / 2.0f);
	double glo = 1.0 / max(f / minf, 1.0f);
	double ghi = 1.0 / max(f / maxf, 1.0f);
	coef.na1 = (1.0 - w) / (w + 1.0);
	coef.b0 = 0.5 * ((glo + ghi) - coef.na1 * (glo - ghi));
	coef.b1 = 0.5 * ((glo - ghi) - coef.na1 * (glo + ghi));
	return coef;
}

//====================================================

//First order low shelf with crossover frequency f and skirt gain multiplier g
fof::coefs lowshelf1(float f, float g)
{
	fof::coefs coef;
	double w = (double)tan(f * constants.pi / 2.0f);
	if (g <= 1) {
		coef.na1 = (1.0 - w) / (w + 1.0);
		coef.b0 = 0.5 * ((g + 1.0) - coef.na1 * (g - 1.0));
		coef.b1 = 0.5 * ((g - 1.0) - coef.na1 * (g + 1.0));
	}
	else {
		g = 1.0f / g;
		coef.b1 = (w - 1.0) / (w + 1.0);
		coef.b0 = 0.5 * ((g + 1.0) + coef.b1 * (g - 1.0));
		coef.na1 = -0.5 * ((g - 1.0) + coef.b1 * (g + 1.0));
		coef.b1 = coef.b1 / coef.b0;
		coef.na1 = coef.na1 / coef.b0;
		coef.b0 = 1.0 / coef.b0;
	}
	return coef;
}

//First order compensated low shelf with crossover f, skirt gain multiplier g,
//and bounds minf, maxf for stability
fof::coefs lowshelf1(float f, float g, float minf, float maxf)
{
	fof::coefs coef;
	double w = (double)tan(bound(f, minf, maxf) * constants.pi / 2);
	if (g <= 1) {
		double glo = (double)max(1.0f / max(f / minf, 1.0f), g);
		double ghi = (double)max(1.0f / max(f / maxf, 1.0f), g);
		coef.na1 = (1.0 - w) / (w + 1.0);
		coef.b0 = 0.5 * ((glo + ghi) - coef.na1 * (glo - ghi));
		coef.b1 = 0.5 * ((glo - ghi) - coef.na1 * (glo + ghi));
	}
	else {
		g = 1.0f / g;
		double glo = (double)max(1.0f / max(f / minf, 1.0f), g);
		double ghi = (double)max(1.0f / max(f / maxf, 1.0f), g);
		coef.b1 = (w - 1.0) / (w + 1.0);
		coef.b0 = 0.5 * ((glo + ghi) + coef.b1 * (glo - ghi));
		coef.na1 = -0.5 * ((glo - ghi) + coef.b1 * (glo + ghi));
		coef.b1 = coef.b1 / coef.b0;
		coef.na1 = coef.na1 / coef.b0;
		coef.b0 = 1.0 / coef.b0;
	}
	return coef;
}

//First order high shelf with crossover frequency f and skirt gain multiplier g
fof::coefs highshelf1(float f, float g)
{
	fof::coefs coef;
	double w = (double)tan(f * constants.pi / 2);
	if (g <= 1) {
		coef.na1 = (1.0 - w) / (w + 1.0);
		coef.b0 = 0.5 * ((1 + g) - coef.na1 * (1 - g));
		coef.b1 = 0.5 * ((1 - g) - coef.na1 * (1 + g));
	}
	else {
		g = 1.0f / g;
		coef.b1 = (w - 1.0) / (w + 1.0);
		coef.b0 = 0.5 * ((1.0 + g) + coef.b1 * (1.0 - g));
		coef.na1 = -0.5 * ((1.0 - g) + coef.b1 * (1.0 + g));
		coef.b1 = coef.b1 / coef.b0;
		coef.na1 = coef.na1 / coef.b0;
		coef.b0 = 1.0 / coef.b0;
	}
	return coef;
}

//First order compensated high shelf with crossover f, skirt gain multiplier g,
//and bounds minf, maxf for stability
fof::coefs highshelf1(float f, float g, float minf, float maxf)
{
	fof::coefs coef;
	double w = (double)tan(bound(f, minf, maxf) * constants.pi / 2);
	if (g <= 1) {
		double glo = bound(f / minf, g, 1.0f);
		double ghi = bound(f / maxf, g, 1.0f);
		coef.na1 = (1.0 - w) / (w + 1.0);
		coef.b0 = 0.5 * ((glo + ghi) - coef.na1 * (glo - ghi));
		coef.b1 = 0.5 * ((glo - ghi) - coef.na1 * (glo + ghi));
	}
	else {
		g = 1.0f / g;
		double glo = bound(f / minf, g, 1.0f);
		double ghi = bound(f / maxf, g, 1.0f);
		coef.b1 = (w - 1.0) / (w + 1.0);
		coef.b0 = 0.5 * ((glo + ghi) + coef.b1 * (glo - ghi));
		coef.na1 = -0.5 * ((glo - ghi) + coef.b1 * (glo + ghi));
		coef.b1 = coef.b1 / coef.b0;
		coef.na1 = coef.na1 / coef.b0;
		coef.b0 = 1.0 / coef.b0;
	}
	return coef;
}

//====================================================

//First order Pade delay approximation for del in samples. Del
//must be >0 for stability.
fof::coefs pade1(float del)
{
	fof::coefs coef;
	coef.na1 = (del - 1.0) / (1.0 + del);
	coef.b0 = -coef.na1;
	coef.b1 = 1.0;
	return coef;
}

//First order allpass with crossover frequency f, phase >f is 180 deg
fof::coefs allpass1(float f)
{
	fof::coefs coef;
	double w = (double)tan(f * constants.pi / 2.0);
	coef.na1 = (1.0 - w) / (w + 1.0);
	coef.b0 = -coef.na1;
	coef.b1 = 1.0;
	return coef;
}

//Simple one-pole smoother applied as a first order filter
fof::coefs smooth1(float N)
{
	fof::coefs coef;
	coef.na1 = exp(-1.0 / N);
	coef.b0 = 1.0 - coef.na1;
	coef.b1 = 0.0f;
	return coef;
}

}


//==================================================

sof::sof(int slew, int maxNch, bool blockTransposed) :
	multialg(maxNch, blockTransposed),
	sc(5, slew), vmem1(maxNch), vmem2(maxNch)
{
	setCoefs(coefs{ 0.0, 0.0, 1.0, 0.0, 0.0 }, true); //init to bypass filter
	clear();
}
sof::~sof() {}

void sof::reallocate(int slew, int maxNch, bool blockTransposed)
{
	multialg::reallocate(maxNch, blockTransposed);
	sc.setSlew(slew);
	vmem1.reset(maxNch);
	vmem2.reset(maxNch);
	clear();
}

void sof::setCoefs(coefs newCoefs, bool converge)
{
	double* pc = sc.target();
	c = newCoefs; //save a reference copy
	pc[0] = c.na1;
	pc[1] = c.na2;
	pc[2] = c.b0;
	pc[3] = c.b1;
	pc[4] = c.b2;
	if (converge) {
		sc.converge();
	}
}

sof::coefs sof::getCoefs() const
{
	return c;
}

void sof::clear()
{
	mem1L = 0.0; //mono/stereo
	mem1R = 0.0;
	mem2L = 0.0;
	mem2R = 0.0;

	vset(vmem1.ptr(), 0.0, vmem1.size()); //multichannel
	vset(vmem2.ptr(), 0.0, vmem1.size());
}

float sof::step(float x)
{
	sc.slew();
	const double* ct = sc.get();
	double temp = x + ct[0] * mem1L + ct[1] * mem2L;
	x = (float)(ct[2] * temp + ct[3] * mem1L + ct[4] * mem2L);
	mem2L = mem1L;
	mem1L = temp;
	return x;
}
void sof::stepBlock(const float* x, float* y, int len)
{
	if (sc.check()) { //slew every sample
		monoalg::stepBlock(x, y, len);
	}
	else { //run faster fixed coefficient alg
		double temp;
		const double* ct = sc.get();
		for (int i = 0; i < len; ++i) {
			temp = x[i] + ct[0] * mem1L + ct[1] * mem2L;
			y[i] = (float)(ct[2] * temp + ct[3] * mem1L + ct[4] * mem2L);
			mem2L = mem1L;
			mem1L = temp;
		}
	}
}

void sof::step(float xL, float xR, float& yL, float& yR)
{
	sc.slew();
	const double* ct = sc.get();
	double tempL = xL + ct[0] * mem1L + ct[1] * mem2L;
	double tempR = xR + ct[0] * mem1R + ct[1] * mem2R;
	yL = (float)(ct[2] * tempL + ct[3] * mem1L + ct[4] * mem2L);
	yR = (float)(ct[2] * tempR + ct[3] * mem1R + ct[4] * mem2R);
	mem2L = mem1L;
	mem1L = tempL;
	mem2R = mem1R;
	mem1R = tempR;
}
void sof::stepBlock(const float* xL, const float* xR, float* yL, float* yR, int len)
{
	if (sc.check()) { //slew every sample
		stereoalg::stepBlock(xL, xR, yL, yR, len);
	}
	else { //run faster fixed coefficient alg
		double tempL, tempR;
		const double* ct = sc.get();
		for (int i = 0; i < len; i++) {
			tempL = xL[i] + ct[0] * mem1L + ct[1] * mem2L;
			tempR = xR[i] + ct[0] * mem1R + ct[1] * mem2R;
			yL[i] = (float)(ct[2] * tempL + ct[3] * mem1L + ct[4] * mem2L);
			yR[i] = (float)(ct[2] * tempR + ct[3] * mem1R + ct[4] * mem2R);
			mem2L = mem1L;
			mem1L = tempL;
			mem2R = mem1R;
			mem1R = tempR;
		}
	}
}

void sof::step(const float* x, float* y, int chan)
{
	sc.slew();
	const double* ct = sc.get();
	double temp = 0.0;
	for (int i = 0; i < chan; ++i) {
		temp = x[i] + ct[0] * vmem1[i] + ct[1] * vmem2[i];
		y[i] = (float)(ct[2] * temp + ct[3] * vmem1[i] + ct[4] * vmem2[i]);
		vmem2[i] = vmem1[i];
		vmem1[i] = temp;
	}
}
void sof::stepBlock(const float* const* x, float* const* y, int chan, int len)
{
	if (sc.check()) { //slew every sample
		multialg::stepBlock(x, y, chan, len);
	}
	else { //run faster fixed coefficient alg
		const double* ct = sc.get();
		double temp = 0.0;
		const float* px = nullptr;
		float* py = nullptr;
		if (transposed) {
			for (int i = 0; i < len; ++i) {
				px = x[i];
				py = y[i];
				for (int j = 0; j < chan; ++j) {
					temp = px[j] + ct[0] * vmem1[j] + ct[1] * vmem2[j];
					py[j] = (float)(ct[2] * temp + ct[3] * vmem1[j] + ct[4] * vmem2[j]);
					vmem2[j] = vmem1[j];
					vmem1[j] = temp;
				}
			}
		}
		else {
			double tmem1 = 0.0;
			double tmem2 = 0.0;
			for (int i = 0; i < chan; ++i) {
				px = x[i];
				py = y[i];
				tmem1 = vmem1[i];
				tmem2 = vmem2[i];
				for (int j = 0; j < len; ++j) {
					temp = px[j] + ct[0] * tmem1 + ct[1] * tmem2;
					py[j] = (float)(ct[2] * temp + ct[3] * tmem1 + ct[4] * tmem2);
					tmem2 = tmem1;
					tmem1 = temp;
				}
				vmem1[i] = tmem1;
				vmem2[i] = tmem2;
			}
		}
	}
}

//====================================================

namespace filt {

// Second order filter definitions
//====================================================
// 
//Inverts the phase of the passed in coefficients structure
sof::coefs invert(const sof::coefs& coef)
{
	return { coef.na1, coef.na2, -coef.b0, -coef.b1, -coef.b2 };
}

//Returns the filter inverse of a stable minimum phase coefs structure
sof::coefs inverse(const sof::coefs& coef)
{
	return { -coef.b1 / coef.b0, -coef.b2 / coef.b0, 1.0 / coef.b0, -coef.na1 / coef.b0, -coef.na2 / coef.b0 };
}

//Merge two first order filters into a single second order one
sof::coefs merge(const fof::coefs& coef1, const fof::coefs& coef2)
{
	sof::coefs coef;
	coef.na1 = coef1.na1 + coef2.na1;
	coef.na2 = -coef1.na1 * coef2.na1;
	coef.b0 = coef1.b0 * coef2.b0;
	coef.b1 = coef1.b0 * coef2.b1 + coef1.b1 * coef2.b0;
	coef.b2 = coef1.b1 * coef2.b1;
	return coef;
}

//Applies bilinear warping so that content at frequency f1 appears at
//frequency f2 in the output coefficient set
sof::coefs warp(const sof::coefs& coef, float f1, float f2)
{
	sof::coefs ret;

	double A = tan(0.5 * constants.pi * f1) / tan(0.5 * constants.pi * f2);
	A = (1.0 - A) / (1.0 + A);
	double AA = A * A;
	double A2 = 2.0 * A;
	double AAp1 = AA + 1.0;
	double a0 = 1.0 - A * coef.na1 - AA * coef.na2;

	ret.b0 = (coef.b0 + A * coef.b1 + AA * coef.b2) / a0;
	ret.b1 = (A2 * coef.b0 + AAp1 * coef.b1 + A2 * coef.b2) / a0;
	ret.b2 = (AA * coef.b0 + A * coef.b1 + coef.b2) / a0;
	ret.na1 = (-A2 + AAp1 * coef.na1 + A2 * coef.na2) / a0;
	ret.na2 = (-AA + A * coef.na1 + coef.na2) / a0;

	return ret;
}

//Returns an allpass filter that phase matches a double application of
//the provided filter with one stage's numerator coefficients flipped
sof::coefs matchphase(const sof::coefs& coef)
{
	return { coef.na1, coef.na2, -coef.na2, -coef.na1, 1.0 };
}

//Flips the numerator coefficients of the provided filter, switching
//min phase <-> max phase, usually for phase aligned double filtering
sof::coefs mp2mp(const sof::coefs& coef)
{
	return { coef.na1, coef.na2, coef.b2, coef.b1, coef.b0 };
}

//====================================================

//Second order unity filter, returns default
sof::coefs unity2()
{
	return { 0.0, 0.0, 1.0, 0.0, 0.0 };
}

//Second order trivial gain filter with gain multiplier g
sof::coefs gain2(float g)
{
	return { 0.0, 0.0, g, 0.0, 0.0 };
}

//2-pole (second order) smoother with time constant of roughly L
//samples and resonance of r >= 1. 
sof::coefs smoother2(float L, float r)
{
	sof::coefs coef;

	L *= 0.5f;
	float temp1 = expf(-1.0f / (L * r));
	float temp2 = sqrtf(1.0f - 1.0f / (r * r));

	coef.na1 = 2.0 * temp1 * cosf(temp2 / L);
	coef.na2 = (double)(-temp1 * temp1);
	coef.b0 = (1.0 - coef.na1 - coef.na2);
	coef.b1 = 0.0; //2-pole smoother, no zeros
	coef.b2 = 0.0;

	return coef;
}

//====================================================

//Second order lowpass with rolloff frequency f and resonance
//multiplier r (r = Q).
sof::coefs lowpass2(float f, float r)
{
	sof::coefs coef;

	double wC = (double)cosf(constants.pi * f);
	double alpha = (double)(sinf(constants.pi * f) / (2.0f * r));
	double a0 = 1.0 + alpha;

	coef.b1 = (1.0 - wC) / a0;
	coef.b0 = 0.5 * coef.b1;
	coef.b2 = coef.b0;
	coef.na1 = 2.0 * wC / a0;
	coef.na2 = (alpha - 1.0) / a0;

	return coef;
}

//Second order lowpass with rolloff frequency f, resonance
//multiplier r (r = Q), and skirt depth multiplier d
sof::coefs lowpass2(float f, float r, float d)
{
	sof::coefs coef;

	double wC = (double)cosf(constants.pi * f);
	double alpha = (double)(sinf(constants.pi * f) / (2.0f * r));

	double wCb = d * (1.0 + wC) / (1.0 - wC);
	wCb = (wCb - 1.0) / (wCb + 1.0);
	double beta = sqrt(1.0 - wCb * wCb) / (2.0 * r);
	double g = (1.0 - wC) / (1.0 - wCb);

	double a0 = 1 + alpha;
	coef.na1 = 2.0 * wC / a0;
	coef.na2 = (alpha - 1.0) / a0;
	coef.b0 = g * (1.0 + beta) / a0;
	coef.b1 = -2.0 * g * wCb / a0;
	coef.b2 = g * (1.0 - beta) / a0;

	return coef;
}

//Second order compensated lowpass with rolloff frequency f, resonance
//multiplier r, bounds minf, maxf for stability, and optional skirt depth
//multiplier d
sof::coefs lowpass2(float f, float r, float minf, float maxf, float d)
{
	sof::coefs coef;

	float w = bound(f, minf, maxf);
	double wC = (double)cosf(constants.pi * w);
	double alpha = (double)(sinf(constants.pi * w) / (2.0f * r));

	double g = (double)min(f / minf, 1.0f);
	g = max(g * g, (double)d);
	double wCb = (double)(w / maxf);
	wCb = max(wCb * wCb, d / (g + constants.eps));

	wCb *= (1.0 + wC) / (1.0 - wC);
	wCb = (wCb - 1.0) / (wCb + 1.0);
	double beta = sqrt(1.0 - wCb * wCb) / (2.0 * r);
	g *= (1.0 - wC) / (1.0 - wCb);

	double a0 = 1 + alpha;
	coef.na1 = 2.0 * wC / a0;
	coef.na2 = (alpha - 1.0) / a0;
	coef.b0 = g * (1.0 + beta) / a0;
	coef.b1 = -2.0 * g * wCb / a0;
	coef.b2 = g * (1.0 - beta) / a0;

	return coef;
}

//Second order highpass with rolloff frequency f and resonance
//multiplier r (r = Q).
sof::coefs highpass2(float f, float r)
{
	sof::coefs coef;

	double wC = (double)cosf(constants.pi * f);
	double alpha = (double)(sinf(constants.pi * f) / (2.0f * r));
	double a0 = 1.0 + alpha;

	coef.b1 = -(1 + wC) / a0;
	coef.b0 = -0.5 * coef.b1;
	coef.b2 = coef.b0;
	coef.na1 = 2.0 * wC / a0;
	coef.na2 = (alpha - 1.0) / a0;

	return coef;
}

//Second order highpass with rolloff frequency f, resonance
//multiplier r (r = Q), and skirt depth multiplier d
sof::coefs highpass2(float f, float r, float d)
{
	sof::coefs coef;

	double wC = (double)cosf(constants.pi * f);
	double alpha = (double)(sinf(constants.pi * f) / (2.0f * r));

	double wCb = d * (1.0 - wC) / (1.0 + wC);
	wCb = (wCb - 1.0) / (wCb + 1.0);
	double beta = sqrt(1.0 - wCb * wCb) / (2.0 * r);
	double g = (1.0 + wC) / (1.0 - wCb);

	double a0 = 1 + alpha;
	coef.na1 = 2.0 * wC / a0;
	coef.na2 = (alpha - 1.0) / a0;
	coef.b0 = g * (1.0 + beta) / a0;
	coef.b1 = 2.0 * g * wCb / a0;
	coef.b2 = g * (1.0 - beta) / a0;

	return coef;
}

//Second order compensated highpass with rolloff frequency f, resonance
//multiplier r, bounds minf, maxf for stability, and optional skirt depth
//multiplier d
sof::coefs highpass2(float f, float r, float minf, float maxf, float d)
{
	sof::coefs coef;

	float w = bound(f, minf, maxf);
	double wC = (double)cosf(constants.pi * w);
	double alpha = (double)(sinf(constants.pi * w) / (2.0f * r));

	double g = 1.0 / max(f / maxf, 1.0f);
	g = max(g * g, (double)d);
	double wCb = (double)(minf / w);
	wCb = max(wCb * wCb, d / (g + constants.eps));

	wCb *= (1.0 - wC) / (1.0 + wC);
	wCb = (wCb - 1.0) / (wCb + 1.0);
	double beta = sqrt(1.0 - wCb * wCb) / (2.0 * r);
	g *= (1.0 + wC) / (1.0 - wCb);

	double a0 = 1 + alpha;
	coef.na1 = 2 * wC / a0;
	coef.na2 = (alpha - 1.0) / a0;
	coef.b0 = g * (1.0 + beta) / a0;
	coef.b1 = 2.0 * g * wCb / a0;
	coef.b2 = g * (1.0 - beta) / a0;

	return coef;
}

//Second order lowpass/highpass pair with shared frequency and r of 0.5,
//forming a Linkwitz-Riley crossover that sums coherently to flat.
void crossover2(sof::coefs& lpf, sof::coefs& hpf, float f)
{
	double wC = (double)cosf(constants.pi * f);
	double alpha = (double)sinf(constants.pi * f); //r = 0.5

	double a0 = 1.0 + alpha;

	lpf.b1 = (1.0 - wC) / a0;
	lpf.b0 = 0.5 * lpf.b1;
	lpf.b2 = lpf.b0;
	lpf.na1 = 2.0 * wC / a0;
	lpf.na2 = (alpha - 1.0) / a0;

	hpf.b1 = (1 + wC) / a0;
	hpf.b0 = -0.5 * hpf.b1;
	hpf.b2 = hpf.b0;
	hpf.na1 = lpf.na1;
	hpf.na2 = lpf.na2;
}

//Second order compensated lowpass/highpass pair with shared frequency
//and bounds and r of 0.5, forming a Linkwitz-Riley crossover that sums
//coherently to flat, and that approaches true single-band operation as
//the frequency approaches the bounds.
void crossover2(sof::coefs& lpf, sof::coefs& hpf, float f, float minf, float maxf)
{
	//common denominator
	f = bound(f, minf, maxf);
	double gamma = (double)tan(f * constants.pi / 2.0f);
	double na1 = (1.0 - gamma) / (1.0 + gamma); //of 1st order filter
	lpf.na1 = 2.0 * na1;
	lpf.na2 = -na1 * na1;
	hpf.na1 = lpf.na1;
	hpf.na2 = lpf.na2;

	//find extreme root gains for each filter
	double gl2h = (double)(f / maxf);
	double gh2l = (double)(minf / f);
	double gl2l = 1.0f;
	double gh2h = 1.0f;
	if (gl2h >= gh2l) { //f is nearer to maxf, lpf holds more power
		gh2h = sqrt(1.0 - gl2h * gl2h);
		gh2l = min(gh2l, gh2h);
		gl2l = sqrt(1.0 - gh2l * gh2l);
	}
	else { //f is nearer to minf, hpf holds more power
		gl2l = sqrt(1.0f - gh2l * gh2l);
		gl2h = min(gl2h, gl2l);
		gh2h = sqrt(1.0f - gl2h * gl2h);
	}

	//LPF linear phase numerator
	double b0 = 0.5 * ((gl2l + gl2h) - na1 * (gl2l - gl2h));
	double b1 = 0.5 * ((gl2l - gl2h) - na1 * (gl2l + gl2h));
	lpf.b0 = b0 * b1;
	lpf.b1 = b0 * b0 + b1 * b1;
	lpf.b2 = lpf.b0;

	//HPF linear phase numerator
	b0 = 0.5 * ((gh2l + gh2h) - na1 * (gh2l - gh2h));
	b1 = 0.5 * ((gh2l - gh2h) - na1 * (gh2l + gh2h));
	hpf.b0 = b0 * b1;
	hpf.b1 = b0 * b0 + b1 * b1;
	hpf.b2 = hpf.b0;
}

//====================================================

//Second order bandpass with center frequency f and octave bandwidth oct
sof::coefs bandpass2(float f, float oct)
{
	sof::coefs coef;

	oct = powf(2.0, -oct / 2.0f); //convert octaves to Q
	double wC = (double)cosf(constants.pi * f);
	double alpha = (double)(sinf(constants.pi * f) / (2.0f * oct));
	double a0 = 1.0 + alpha;

	coef.b0 = alpha / a0;
	coef.b1 = 0.0;
	coef.b2 = -coef.b0;
	coef.na1 = 2 * wC / a0;
	coef.na2 = (alpha - 1.0) / a0;

	return coef;
}

//TODO -- COMPENSATED BANDPASS

//Second order resonator with center frequency f and resonance multiplier r
sof::coefs resonator2(float f, float r)
{
	sof::coefs coef;

	double wC = (double)cosf(constants.pi * f);
	double alpha = (double)(sinf(constants.pi * f) / (2.0f * r));
	double a0 = 1.0 + alpha;

	coef.b0 = r * alpha / a0;
	coef.b1 = 0.0;
	coef.b2 = -coef.b0;
	coef.na1 = 2 * wC / a0;
	coef.na2 = (alpha - 1.0) / a0;

	return coef;
}

//TODO -- COMPENSATED RESONATOR

//====================================================

//Second order low shelf filter with crossover frequency f, quality factor
//Q, and low frequency linear gain g
sof::coefs lowshelf2(float f, float Q, float g)
{
	sof::coefs coef;

	double wC = (double)cosf(constants.pi * f);
	double alpha = (double)(sinf(constants.pi * f) / (2.0f * Q));
	double A = (double)sqrtf(g);
	double As = sqrt(A);
	double a0 = (A + 1.0) + (A - 1.0) * wC + 2 * As * alpha;

	coef.b0 = A * ((A + 1.0) - (A - 1.0) * wC + 2.0 * As * alpha) / a0;
	coef.b1 = A * (2.0 * ((A - 1.0) - (A + 1.0) * wC)) / a0;
	coef.b2 = A * ((A + 1.0) - (A - 1.0) * wC - 2.0 * As * alpha) / a0;
	coef.na1 = (2.0 * ((A - 1.0) + (A + 1.0) * wC)) / a0;
	coef.na2 = -((A + 1.0) + (A - 1.0) * wC - 2.0 * As * alpha) / a0;

	return coef;
}

//TODO -- COMPENSATED LOWSHELF

//Second order high shelf filter with crossover frequency f, quality factor
//Q, and high frequency linear gain g
sof::coefs highshelf2(float f, float Q, float g)
{
	sof::coefs coef;

	double wC = (double)cosf(constants.pi * f);
	double alpha = (double)(sinf(constants.pi * f) / (2.0f * Q));
	double A = (double)sqrtf(g);
	double As = sqrt(A);
	double a0 = (A + 1.0) - (A - 1.0) * wC + 2 * As * alpha;

	coef.b0 = A * ((A + 1.0) + (A - 1.0) * wC + 2.0 * As * alpha) / a0;
	coef.b1 = A * (-2.0 * ((A - 1.0) + (A + 1.0) * wC)) / a0;
	coef.b2 = A * ((A + 1.0) + (A - 1.0) * wC - 2.0 * As * alpha) / a0;
	coef.na1 = (-2.0 * ((A - 1.0) - (A + 1.0) * wC)) / a0;
	coef.na2 = -((A + 1.0) - (A - 1.0) * wC - 2.0 * As * alpha) / a0;

	return coef;
}

//TODO -- COMPENSATED HIGHSHELF

//====================================================

//Second order boost/cut peaking filter with center frequency f, quality
//factor Q, and linear peak gain g
sof::coefs peaking2(float f, float Q, float g)
{
	sof::coefs coef;

	double wC = (double)cosf(constants.pi * f);
	double alpha = (double)(sinf(constants.pi * f) / (2.0f * Q));
	double A = (double)sqrtf(g);
	double a0 = 1.0 + alpha / A;

	coef.b0 = (1.0 + alpha * A) / a0;
	coef.b1 = -2.0 * wC / a0;
	coef.b2 = (1.0 - alpha * A) / a0;
	coef.na1 = -coef.b1;
	coef.na2 = -(1.0 - alpha / A) / a0;

	return coef;
}

//TODO -- COMPENSATED... PEAKING...?

//====================================================

//TODO -- PHASE

}

//==================================================


namespace filt {

// Second order filter cascade definitions
//====================================================

//Fourth order lowpass with rolloff frequency f, resonance
//multiplier r (r = Q), and spread factor s where 0 is tightest
//and 1 is widest (s of 1 is a double filter)
const sofcasc<2>::coefs& lowpass4(sofcasc<2>::coefs& coef, float f, float r, float s)
{
	double wC = (double)cosf(constants.pi * f); //constants
	double wS = (double)sinf(constants.pi * f);

	//derive the Q factors for each filter
	double alpha1 = (double)sqrtf(r * (s + (1.0f - s) * (2.0f * r + sqrtf(4.0f * r * r - 1.0f))));
	double alpha2 = r / alpha1;
	alpha1 = wS / (2.0 * alpha1);
	alpha2 = wS / (2.0 * alpha2);

	double a0 = 1.0 + alpha1; //filter 1
	coef[0].b1 = (1.0 - wC) / a0;
	coef[0].b0 = 0.5 * coef[0].b1;
	coef[0].b2 = coef[0].b0;
	coef[0].na1 = 2.0 * wC / a0;
	coef[0].na2 = (alpha1 - 1.0) / a0;

	a0 = 1.0 + alpha2; //filter 2
	coef[1].b1 = (1.0 - wC) / a0;
	coef[1].b0 = 0.5 * coef[1].b1;
	coef[1].b2 = coef[1].b0;
	coef[1].na1 = 2.0 * wC / a0;
	coef[1].na2 = (alpha2 - 1.0) / a0;

	return coef;
}

//Fourth order lowpass with rolloff frequency f, resonance
//multiplier r (r = Q), spread factor s on [0, 1], and skirt
//depth multiplier d
const sofcasc<2>::coefs& lowpass4(sofcasc<2>::coefs& coef, float f, float r, float s, float d)
{
	double wC = (double)cosf(constants.pi * f); //constants
	double wS = (double)sinf(constants.pi * f);

	//derive the Q factors for each filter
	double alpha1 = (double)sqrtf(r * (s + (1.0f - s) * (2.0f * r + sqrtf(4.0f * r * r - 1.0f))));
	double alpha2 = r / alpha1;

	d = sqrtf(max(d, 0.0f));
	double wCb = d * (1.0 + wC) / (1.0 - wC);
	wCb = (wCb - 1.0) / (wCb + 1.0);
	double g = (1.0 - wC) / (1.0 - wCb);

	double beta2 = 0.5 * sqrt(1.0 - wCb * wCb);
	double beta1 = beta2 / alpha1;
	beta2 /= alpha2;
	alpha1 = wS / (2.0 * alpha1);
	alpha2 = wS / (2.0 * alpha2);

	double a0 = 1 + alpha1; //filter 1
	coef[0].na1 = 2 * wC / a0;
	coef[0].na2 = (alpha1 - 1.0) / a0;
	coef[0].b0 = g * (1.0 + beta1) / a0;
	coef[0].b1 = -2.0 * g * wCb / a0;
	coef[0].b2 = g * (1.0 - beta1) / a0;

	a0 = 1 + alpha2; //filter 2
	coef[1].na1 = 2 * wC / a0;
	coef[1].na2 = (alpha2 - 1.0) / a0;
	coef[1].b0 = g * (1.0 + beta2) / a0;
	coef[1].b1 = -2.0 * g * wCb / a0;
	coef[1].b2 = g * (1.0 - beta2) / a0;

	return coef;
}

//Fourth order compensated lowpass with rolloff frequency f, resonance
//multiplier r, spread factor s on [0, 1], bounds minf, maxf for stability,
//and optional skirt depth multiplier d
const sofcasc<2>::coefs& lowpass4(sofcasc<2>::coefs& coef, float f, float r, float s,
	float minf, float maxf, float d)
{
	double w = bound(f, minf, maxf);
	double wC = (double)cosf(constants.pi * (float)w);
	double wS = (double)sinf(constants.pi * (float)w);

	//derive the Q factors for each filter
	double alpha1 = (double)sqrtf(r * (s + (1.0f - s) * (2.0f * r + sqrtf(4.0f * r * r - 1.0f))));
	double alpha2 = r / alpha1;

	d = sqrtf(max(d, 0.0f));
	double g = (double)min(f / minf, 1.0f);
	g = max(g * g, (double)d);
	double wCb = (double)(w / maxf);
	wCb = max(wCb * wCb, d / (g + constants.eps));

	wCb *= (1.0 + wC) / (1.0 - wC);
	wCb = (wCb - 1.0) / (wCb + 1.0);
	g *= (1.0 - wC) / (1.0 - wCb);
	double beta2 = 0.5 * sqrt(1.0 - wCb * wCb);
	double beta1 = beta2 / alpha1;
	beta2 /= alpha2;
	alpha1 = wS / (2.0 * alpha1);
	alpha2 = wS / (2.0 * alpha2);

	double a0 = 1 + alpha1; //filter 1
	coef[0].na1 = 2 * wC / a0;
	coef[0].na2 = (alpha1 - 1.0) / a0;
	coef[0].b0 = g * (1.0 + beta1) / a0;
	coef[0].b1 = -2.0 * g * wCb / a0;
	coef[0].b2 = g * (1.0 - beta1) / a0;

	a0 = 1 + alpha2; //filter 2
	coef[1].na1 = 2 * wC / a0;
	coef[1].na2 = (alpha2 - 1.0) / a0;
	coef[1].b0 = g * (1.0 + beta2) / a0;
	coef[1].b1 = -2.0 * g * wCb / a0;
	coef[1].b2 = g * (1.0 - beta2) / a0;

	return coef;
}

//====================================================

//Fourth order bandpass with center frequency f and octave bandwidth oct derived
//from a frequency matched pair of lowpass/highpass filters with unity gain at
//the center frequency
const sofcasc<2>::coefs& bandpass4(sofcasc<2>::coefs& coef, float f, float oct)
{
	oct = powf(2.0, -oct / 2.0f); //convert octaves to Q
	double wC = (double)cosf(constants.pi * f);
	double alpha = (double)(sinf(constants.pi * f) / (2.0f * oct));
	double a0 = 1.0 + alpha;

	coef[0].b1 = (1.0 - wC) / (a0 * oct); //lowpass section
	coef[0].b0 = 0.5 * coef[0].b1;
	coef[0].b2 = coef[0].b0;
	coef[0].na1 = 2.0 * wC / a0;
	coef[0].na2 = (alpha - 1.0) / a0;

	coef[1].b1 = -(1 + wC) / (a0 * oct); //highpass section
	coef[1].b0 = -0.5 * coef[1].b1;
	coef[1].b2 = coef[1].b0;
	coef[1].na1 = coef[0].na1;
	coef[1].na2 = coef[0].na2;

	return coef;
}

//====================================================

//8th order (4 section) antialiasing lowpass filter with a bandwidth
//of 1/DSR the Nyquist rate
const sofcasc<4>::coefs& antialias8(sofcasc<4>::coefs& coef, int DSR)
{
	//8th order elliptic filter predesigned at 1/4 sampling rate
	coef[0] = { 0.890565605450, -0.252020135666, 0.010764012966, 0.019135272819, 0.010764012966 };
	coef[1] = { 0.704498870114, -0.490393099038, 1.000000000000, 0.816073294622, 1.000000000000 };
	coef[2] = { 0.523747331801, -0.743421743924, 1.000000000000, 0.224358160926, 1.000000000000 };
	coef[3] = { 0.442544252241, -0.923764924741, 1.000000000000, 0.003227023719, 1.000000000000 };

	if (DSR > 2) { //then warp it in place
		warp<4>(coef, 0.5f, 1.0f / DSR);
	}
	return coef;
}

}

}