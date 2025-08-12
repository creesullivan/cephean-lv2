
//--------------------------------------------------
// Defines DSP classes and functions
//--------------------------------------------------

#include "cephean-dsp.h"

namespace cephean
{

//==================================================

//Periodic triangle wave from -1 to +1 as a function of phase x
//where triwave(0) = triwave(2*pi) = 0.0f ala a sine wave.
float triwave(float x)
{
	x /= (2.0f * constants.pi);
	x -= 0.25f;
	x -= floorf(x); //x between 0 and 1
	x = (x <= 0.5f) ? (1.0f - 4.0f * x) : (4.0f * x - 3.0f);
	return x;
}

//Periodic quadratic wave from -1 to +1 as a function of phase x
//where quadwave(0) = quadwave(2*pi) = 0.0f ala a sine wave.
float quadwave(float x)
{
	x /= (2.0f * constants.pi);
	x += 0.25f;
	x -= floorf(x); //x between 0 and 1
	x = fabsf(2.0f*x - 1.0f); //x between 0 and 1
	if (x <= 0.5f) {
		x = 1.0f - 4.0f * x * x;
	}
	else {
		x = 1.0f - x;
		x = 4.0f * x * x - 1.0f;
	}
	return x;
}

//=====================================================

fastcircint::fastcircint(int len, int vset) : wrap(len)
{
	set(vset);
}
fastcircint::fastcircint(const fastcircint& other) : fastcircint(other.wrap, other.val)
{}
fastcircint::~fastcircint() {}

void fastcircint::set(int vset)
{
	if (vset >= wrap) {
		val = vset - wrap;
	}
	else if (vset < 0) {
		val = vset + wrap;
	}
	else {
		val = vset;
	}
}
void fastcircint::setlength(int len, setlenmode mode)
{
	wrap = len;
	switch (mode) {
	case setlenmode::WRAP:
		val %= wrap; //this behavior might break our fast assumptions
		if (val < 0) {
			val += wrap;
		}
		break;
	case setlenmode::BOUND:
		if (val >= wrap) {
			val = 0; //clear only if out of range
		}
		break;
	default: //ZERO
		val = 0; //simply clear
	}
	
}
void fastcircint::operator=(int vset)
{
	set(vset);
}

int fastcircint::length() const
{
	return wrap;
}
int fastcircint::get() const
{
	return val;
}
fastcircint::operator int() const
{
	return val;
}

void fastcircint::operator+=(int i)
{
	set(val + i);
}
void fastcircint::operator-=(int i)
{
	set(val - i);
}
void fastcircint::operator+=(unsigned int i)
{
	val += i;
	if (val >= wrap) {
		val -= wrap;
	}
}
void fastcircint::operator-=(unsigned int i)
{
	val -= i;
	if (val < 0) {
		val += wrap;
	}
}
int fastcircint::operator++()
{
	if (++val == wrap) {
		val = 0;
	} //preincrement
	return val;
}
int fastcircint::operator++(int)
{
	int ret = val++;
	if (val == wrap) {
		val = 0;
	}  //post-increment
	return ret;
}
int fastcircint::operator--()
{
	val = (val == 0) ? (wrap - 1) : (val - 1);
	return val; //pre-decrement
}
int fastcircint::operator--(int)
{
	int ret = val;
	val = (val == 0) ? (wrap - 1) : (val - 1);
	return ret; //post-decrement
}

fastcircint fastcircint::operator+(int i) const
{
	return fastcircint(wrap, val + i);
}
fastcircint fastcircint::operator-(int i) const
{
	return fastcircint(wrap, val - i);
}


fcivec::fcivec(int len, int wrap) : vec<fastcircint>(len)
{
	for (int i = 0; i < len; ++i) {
		p[i].setlength(wrap);
	}
}
fcivec::~fcivec() {}

void fcivec::setwrap(int wrap, fastcircint::setlenmode mode)
{
	for (int i = 0; i < N; ++i) {
		p[i].setlength(wrap, mode);
	}
}


fastcircfloat::fastcircfloat(float wrapset, float valset) : wrap(wrapset)
{
	set(valset);
}
fastcircfloat::fastcircfloat(const fastcircfloat& other) : fastcircfloat(other.wrap, other.val)
{
}
fastcircfloat::~fastcircfloat() {}

void fastcircfloat::set(float newval)
{
	if (newval >= wrap) {
		val = newval - wrap;
	}
	else if (newval < 0.0f) {
		val = newval + wrap;
	}
	else {
		val = newval;
	}
}
void fastcircfloat::setwrap(float newwrap, setwrapmode mode)
{
	wrap = newwrap;
	switch (mode) {
	case setwrapmode::WRAP:
		val = fmodf(val, wrap); //this behavior might break our fast assumptions
		if (val < 0.0f) {
			val += wrap;
		}
		break;
	case setwrapmode::BOUND:
		if (val >= wrap) {
			val = 0.0f; //clear only if out of range
		}
		break;
	default: //ZERO
		val = 0.0f; //simply clear
	}

	
}
void fastcircfloat::operator=(float valset)
{
	set(valset);
}

float fastcircfloat::getwrap() const
{
	return wrap;
}
float fastcircfloat::get() const
{
	return val;
}
fastcircfloat::operator float() const
{
	return val;
}

//inc must be positive! use safeinc() for mixed sign
void fastcircfloat::operator+= (float inc)
{
	val += inc;
	if (val >= wrap) {
		val -= wrap;
	}
}
//dec must be positive! use safeinc() for mixed sign
void fastcircfloat::operator-=(float dec)
{
	val -= dec;
	if (val < 0.0f) {
		val += wrap;
	}
}
//incdec may be positive (increment) or negative (decrement)
void fastcircfloat::safeinc(float incdec)
{
	set(val + incdec);
}

fastcircfloat fastcircfloat::operator+(float inc) const
{
	return fastcircfloat(wrap, val + inc);
}
fastcircfloat fastcircfloat::operator-(float dec) const
{
	return fastcircfloat(wrap, val - dec);
}


fcfvec::fcfvec(int len, float wrap) : vec<fastcircfloat>(len)
{
	for (int i = 0; i < len; ++i) {
		p[i].setwrap(wrap);
	}
}
fcfvec::~fcfvec() {}

void fcfvec::setwrap(float wrap, fastcircfloat::setwrapmode mode)
{
	for (int i = 0; i < N; ++i) {
		p[i].setwrap(wrap, mode);
	}
}


//==================================================

flatbuffer::flatbuffer(int len) : vec<float>(len)
{
	clear();
}
flatbuffer::~flatbuffer() {}

int flatbuffer::size() const
{
	return N;
}
void flatbuffer::reset(int newlen)
{
	vec<float>::reset(newlen);
	clear();
}

void flatbuffer::clear(float val)
{
	vset(p, val, N);
}

void flatbuffer::put(const float* x, int len)
{
	if (len >= N) { //fills up the whole buffer
		vcopy(x + (len - N), p, N);
	}
	else { //shift data back and put new data at the end
		vcopy(p + len, p, N - len);
		vcopy(x, p + (N - len), len);
	}
}

float flatbuffer::get(int del, int hostlen, int hostn) const
{
	return *(p + (N - hostlen - del + hostn));
}
void flatbuffer::get(float* y, int len, int del) const
{
	vcopy(p + (N - len - del), y, len); //copy the data out
}

float* flatbuffer::get(int len, int del)
{
	return (p + (N - len - del)); //return data pointer
}
const float* flatbuffer::get(int len, int del) const
{
	return (p + (N - len - del)); //return data pointer
}

circbuffer::circbuffer(int len) : vec<float>(len),
	ind(len)
{
	clear();
}
circbuffer::~circbuffer() {}

int circbuffer::size() const
{
	return N;
}

void circbuffer::reset(int newlen)
{
	vec<float>::reset(newlen);
	ind.setlength(newlen);
	clear();
}

void circbuffer::clear(float val)
{
	vset(p, val, N);
	ind = 0;
}

void circbuffer::put(float x)
{
	p[++ind] = x;
}
void circbuffer::put(const float* x, int len)
{
	if (len >= N) { //fills up the whole buffer
		vcopy(x + (len - N), p, N);
		ind = N-1; //always points to most recent sample
	}
	else { //overwrite oldest data
		int remlen = N - 1 - ind.get();
		if (len > remlen) { //overruns the boundaries
			vcopy(x, p + N - remlen, remlen);
			vcopy(x + remlen, p, len - remlen);
			ind = len - remlen - 1;
		}
		else { //treat like flat buffer
			ind += len;
			vcopy(x, p + ind.get() - len + 1, len);
		}
	}
}

void circbuffer::set(float x, int del, int hostlen, int hostn)
{
	fastcircint pind = (ind - del);
	pind -= (hostlen - 1 - hostn);
	p[pind] = x;
}
void circbuffer::set(const float* x, int len, int del)
{
	fastcircint pind = (ind - del);
	pind -= (len - 1);

	int remlen = N - 1 - pind.get();
	if (len > remlen) { //overruns the boundaries
		vcopy(x, p + N - remlen, remlen);
		vcopy(x + remlen, p, len - remlen);
	}
	else { //treat like a flat buffer
		vcopy(x, p + pind.get(), len);
	}
}

float circbuffer::get(int del, int hostlen, int hostn) const
{
	fastcircint gind = (ind - del);
	gind -= (hostlen - 1 - hostn);
	return p[gind];
}
void circbuffer::get(float* y, int len, int del) const
{
	fastcircint gind = (ind - del);
	gind -= (len - 1);

	int remlen = N - gind.get();
	if (len > remlen) { //overruns the boundaries
		vcopy(p + N - remlen, y, remlen);
		vcopy(p, y + remlen, len - remlen);
	}
	else { //treat like a flat buffer
		vcopy(p + gind.get(), y, len);
	}
}

float circbuffer::getNN(float del, int hostlen, int hostn) const
{
	return get((int)ceilf(del), hostlen, hostn); //ceil for consistency w/ others
}
//interpolated from a starting delay, adjusting ddel samp/samp for pitch shifting
void circbuffer::getNN(float* y, int len, float del0, float ddel, int hostlen, int hostn) const
{
	if (hostlen < 0) { //automatic rebuffering helpers
		hostlen = len;
	}

	fastcircfloat gindf((float)ind.length(), (ind.get() - hostlen + 1 + hostn) - del0);

	ddel = 1.0f - ddel; //convert to sample step
	if (ddel >= 0.0f) { //then all steps still move forward in time
		for (int n = 0; n < len; ++n) {
			y[n] = p[(int)floorf(gindf)]; //floor for speed and safety
			gindf += ddel;
		}
	}
	else{ //then all steps move backward in time
		ddel = -ddel;
		for (int n = 0; n < len; ++n) {
			y[n] = p[(int)floorf(gindf)];
			gindf -= ddel;
		}
	}
}
void circbuffer::getNN(float* y, const float* del, int len) const
{
	fastcircint gind = ind;
	int gind0 = ind.get() - len + 1;

	for (int n = 0; n < len; ++n) { //no easy way to copy fast for arbitrary del vector
		gind.set(gind0 + n - (int)ceilf(del[n])); //ceil for consistency w/ others
		y[n] = p[gind];
	}
}


//==================================================

gain::gain(int slew) : g(slew, 1.0f)
{}
gain::~gain() {}

void gain::setLevel(float mult, bool converge)
{
	g.target(mult, converge);
}

float gain::step(float x)
{
	g.slew();
	return x * g;
}
void gain::stepBlock(const float* x, float* y, int len)
{
	if (g.check()) { //slew sample by sample
		monoalg::stepBlock(x, y, len);
	}
	else { //high efficiency unslewed operation
		vmult(x, g.get(), y, len);
	}
}


//==================================================

attrel::attrel(float attackInSamples, float releaseInSamples)
{
	setAttack(attackInSamples); //init
	setRelease(releaseInSamples);
	clear();
}
attrel::~attrel() {}

void attrel::setAttack(float samples)
{
	aalpha = (samples > 0.0f) ? expf(-1.0f / samples) : 0.0f;
	abeta = 1.0f - aalpha;
}
void attrel::setRelease(float samples)
{
	ralpha = (samples > 0.0f) ? expf(-1.0f / samples) : 0.0f;
	rbeta = 1.0f - ralpha;
}

void attrel::clear(float val)
{
	env = val;
}
float attrel::step(float x)
{
	if (x >= env) { //attack
		env = aalpha * env + abeta * x;
	}
	else { //release
		env = ralpha * env + rbeta * x;
	}
	return env;
}

vsmo::vsmo(int Nch, bool blockTransposed, float smoothingInSamples) :
	multialg(Nch, blockTransposed), env(Nch)
{
	setSmooth(smoothingInSamples); //init
	clear();
}
vsmo::~vsmo() {}

void vsmo::setSmooth(float samples)
{
	alpha = (samples > 0.0f) ? expf(-1.0f / samples) : 0.0f;
	beta = 1.0f - alpha;
}

void vsmo::clear(float val)
{
	vset(env.ptr(), val, env.size());
}
void vsmo::step(const float* x, float* y, int chan)
{
	for (int i = 0; i < chan; ++i) {
		env[i] = env[i] * alpha + beta * x[i];
	}
	vcopy(env.ptr(), y, chan);
}

vatt::vatt(int Nch, bool blockTransposed, float attackInSamples) :
	multialg(Nch, blockTransposed), env(Nch)
{
	setAttack(attackInSamples); //init
	clear();
}
vatt::~vatt() {}

void vatt::setAttack(float samples)
{
	aalpha = (samples > 0.0f) ? expf(-1.0f / samples) : 0.0f;
	abeta = 1.0f - aalpha;
}

void vatt::clear(float val)
{
	vset(env.ptr(), val, env.size());
}
void vatt::step(const float* x, float* y, int chan)
{
	for (int i = 0; i < chan; ++i) {
		env[i] = min(x[i], env[i] * aalpha + abeta * x[i]);
	}
	vcopy(env.ptr(), y, chan);
}

vrel::vrel(int Nch, bool blockTransposed, float releaseInSamples) :
	multialg(Nch, blockTransposed), env(Nch)
{
	setRelease(releaseInSamples); //init
	clear();
}
vrel::~vrel() {}

void vrel::setRelease(float samples)
{
	ralpha = (samples > 0.0f) ? expf(-1.0f / samples) : 0.0f;
	rbeta = 1.0f - ralpha;
}

void vrel::clear(float val)
{
	vset(env.ptr(), val, env.size());
}
void vrel::step(const float* x, float* y, int chan)
{
	for (int i = 0; i < chan; ++i) {
		env[i] = max(x[i], env[i] * ralpha + rbeta * x[i]);
	}
	vcopy(env.ptr(), y, chan);
}



//==================================================

lma::lma(int slew, double leakInSamples, int maxLengthInSamples, int maxBlockSize) :
	L(leakInSamples), buff(maxLengthInSamples + maxBlockSize), N(slew, 1),
	scr(maxBlockSize)
{
	alpha = exp(-1.0 / L);
	setLength(1, true);
	clear();
}
lma::~lma() {}

void lma::setLength(int samples, bool converge)
{
	N.target(samples, converge);
	if (converge) {
		update();
	}
}

void lma::clear()
{
	buff.clear();
	mem = 0.0;
}
float lma::step(float x)
{
	if (N.slew()) {
		update();
	}
	mem = alpha * mem + (double)x;	//apply fixed leak
	buff.put((float)mem);			//delay output
	mem -= betaN * buff.get(N);		//cancel tail
	return (float)(gain * mem);		//makeup gain
}
void lma::stepBlock(const float* x, float* y, int len)
{
	if (N.check()) { //slew sample by sample
		monoalg::stepBlock(x, y, len);
	}
	else { //run faster block-based algorithm
		for (int i = 0; i < len; ++i) {
			mem = alpha * mem + (double)x[i];
			y[i] = (float)mem; //apply fixed leak
		}
		buff.put(y, len); //delay output
		buff.get(scr, len, N);
		for (int i = 0; i < len; ++i) { //cancel tail with makeup gain
			y[i] = (float)(gain * (y[i] - betaN * scr[i]));
		}
	}
}

void lma::update()
{
	betaN = exp(-N / L);
	gain = (1.0 - alpha) / (1.0 - betaN);
}


hold4::hold4(float clearval) : x0(clearval), n(1, 0)
{
	clear();
}
hold4::~hold4() {}

//Set the quarter-length of the hold process, the full hold = 4*Nq
void hold4::setQuarterLength(int Nq)
{
	n.setlength(Nq, fastcircint::setlenmode::BOUND);
}

void hold4::clear()
{
	n = 0;
	vset(mem, x0, 4);
}

float hold4::step(float x)
{
	//keep track of maxes
	vmax(mem, x, mem, 4);

	float y = mem[0]; //output global max
	if (y == x) {
		n = 0; //reset on new global max
	}

	//step for next time
	if (++n == 0) { //on overrun, shift the memory back
		mem[0] = mem[1];
		mem[1] = mem[2];
		mem[2] = mem[3];
		mem[3] = x0;
	}

	return y;
}

vhold4::vhold4(int Nch, bool blockTransposed, float clearval) : multialg(Nch, blockTransposed),
	x0(clearval), n(Nch, 1), mem(Nch * 4), M(Nch)
{
	clear();
}
vhold4::~vhold4() {}

//Set the quarter-length of the hold process, the full hold = 4*Nq
void vhold4::setQuarterLength(int Nq)
{
	n.setwrap(Nq, fastcircint::setlenmode::BOUND);
}

void vhold4::clear()
{
	vset(mem.ptr(), x0, mem.size());
	vset(n.ptr(), 0, n.size());
}

void vhold4::step(const float* x, float* y, int chan)
{
	UNUSED(chan); //chan better equal M, or we will have issues

	float* mptr = mem.ptr();
	for (int i = 0; i < M; ++i) {
		vmax(mptr, x[i], mptr, 4); //keep track of maxes

		y[i] = mptr[0]; //copy out the global maxes on each channel

		if (y[i] == x[i]) {
			n[i] = 0; //reset on new global max
		}

		if (++n[i] == 0) { //on overrun, shift the memory back
			mptr[0] = mptr[1];
			mptr[1] = mptr[2];
			mptr[2] = mptr[3];
			mptr[3] = x0;
		}

		mptr += 4; //shift
	}
}

//==================================================

fastratio::fastratio(int slew, int maxNch, bool blockTransposed) :
	multialg(maxNch, blockTransposed),
	prop(2, slew)
{
	set(1.0f, 1.0f, true); //init
}

void fastratio::set(float newRatio, float newDrive, bool converge)
{
	newRatio = bound(newRatio, 0.1f, 0.99f); //perfect linearity not supported
	newDrive = max(newDrive, 0.0f);

	prop.target()[0] = newRatio;
	prop.target()[1] = newDrive;
	if (converge) {
		prop.converge();
		update();
	}
}

float fastratio::step(float x)
{
	if (prop.swap()) {
		update();
	}

	x = fabsf(x);
	float y = coefA.a1 / (1.0f + coefA.b1 * x) + coefA.a2 / (1.0f + coefA.b2 * x);

	if (prop.fade()) {
		x = coefB.a1 / (1.0f + coefB.b1 * x) + coefB.a2 / (1.0f + coefB.b2 * x);
		y *= prop.gain();
		y += (1.0f - prop.gain()) * x;
	}
	return y;
}

void fastratio::step(float xL, float xR, float& yL, float& yR)
{
	if (prop.swap()) {
		update();
	}

	xL = fabsf(xL);
	xR = fabsf(xR);
	yL = coefA.a1 / (1.0f + coefA.b1 * xL) + coefA.a2 / (1.0f + coefA.b2 * xL);
	yR = coefA.a1 / (1.0f + coefA.b1 * xR) + coefA.a2 / (1.0f + coefA.b2 * xR);

	if (prop.fade()) {
		xL = coefB.a1 / (1.0f + coefB.b1 * xL) + coefB.a2 / (1.0f + coefB.b2 * xL);
		xR = coefB.a1 / (1.0f + coefB.b1 * xR) + coefB.a2 / (1.0f + coefB.b2 * xR);
		yL *= prop.gain();
		yR *= prop.gain();
		yL += (1.0f - prop.gain()) * xL;
		yR += (1.0f - prop.gain()) * xR;
	}
}

void fastratio::step(const float* x, float* y, int chan)
{
	if (prop.swap()) {
		update();
	}
	bool fading = prop.fade();

	float xtemp = 0.0f;
	float ytemp = 0.0f;
	for (int i = 0; i < chan; ++i) {
		xtemp = fabsf(x[i]);
		y[i] = coefA.a1 / (1.0f + coefA.b1 * xtemp) + coefA.a2 / (1.0f + coefA.b2 * xtemp);

		if (fading) {
			ytemp = coefB.a1 / (1.0f + coefB.b1 * xtemp) + coefB.a2 / (1.0f + coefB.b2 * xtemp);
			y[i] *= prop.gain();
			y[i] += (1.0f - prop.gain()) * ytemp;
		}
	}
}

void fastratio::update()
{
	coefB = coefA; //copy out old values for crossfade

	const float r = prop.current()[0];
	const float g0 = powf(x0, r - 1.0f); //gain as x-> 0
	const float y1 = powf(x1, r);
	const float y2 = powf(x2, r);
	//y3 = x3 = 1

	//form and solve the SOLE
	float eq11 = y1 - x1;
	float eq12 = x1 * (y1 - 1.0f);
	float eq21 = y2 - x2;
	float eq22 = x2 * (y2 - 1.0f);
	float det = eq11 * eq22 - eq12 * eq21;

	float sol1 = g0 * (1.0f - x1) + x1 - (y1 / x1);
	float sol2 = g0 * (1.0f - x2) + x2 - (y2 / x2);

	coefA.b1 = (sol1 * eq22 - sol2 * eq12) / det;
	coefA.b2 = (sol2 * eq11 - sol1 * eq21) / det;

	//convert SOLE solution to waveshaper parameters
	coefA.b2 = 0.5f * (coefA.b1 + sqrtf(coefA.b1 * coefA.b1 - 4.0f * coefA.b2));
	coefA.b1 = coefA.b1 - coefA.b2;
	coefA.a1 = (1.0f + coefA.b1) * ((1.0f + coefA.b2) - g0) / (coefA.b2 - coefA.b1);
	coefA.a2 = g0 - coefA.a1;

	//modify to include drive gain
	const float d = prop.current()[1];
	coefA.b1 *= d;
	coefA.b2 *= d;
	float g3 = coefA.a1 / (1.0f + coefA.b1) + coefA.a2 / (1.0f + coefA.b2);
	coefA.a1 /= g3;
	coefA.a2 /= g3; //enforce g3 = 1, so 1 = 1*1
}

//==================================================

//Linear interpolation at fractional sample indices ind, forming
// y = x[ind]. Extrapolation is unsafe! Pre-bound ind to [0, xlen-1].
//
// Interpolation is unsafe in place, but ind and y can share memory.
void interp(const float* x, const float* ind, float* y, int ylen)
{
	int ind1, ind2;
	float g1, g2;
	for (int i = 0; i < ylen; ++i) {
		g1 = floorf(ind[i]);
		ind1 = (int)g1;
		ind2 = (int)ceilf(ind[i]);
		g2 = ind[i] - g1;
		g1 = 1.0f - g2;

		y[i] = x[ind1] * g1 + x[ind2] * g2;
	}
}


interpolator::interpolator(int ylen, int slew) :
	N(ylen), ind(ylen, slew), g1(ylen), g2(ylen), i1(ylen), i2(ylen)
{
	vset(ind.target(), 0.0f, ind.size()); //init simply
	update();
}

interpolator::~interpolator() {}

//Sets the index map so that y = x[index], len == ylen and the indices
//should be bounded by 0 and xlen - 1 for safety
void interpolator::set(const float* newIndex, int len, bool converge)
{
	vcopy(newIndex, ind.target(), len);
	if (converge) {
		ind.converge();
		update();
	}
}

//Returns a pointer to the index map for direct modification. Converge
//manually if desired with the converge() function.
float* interpolator::set()
{
	return ind.target();
}

void interpolator::slewpause(bool paused)
{
	paused ? ind.pause() : ind.unpause();
}

void interpolator::converge()
{
	ind.converge();
	update();
}

//Applies the previously set index map. x == y is UNSAFE
void interpolator::apply(const float* x, float* y, int xlen, int ylen)
{
	if (ind.slew()) {
		update();
	}
	for (int i = 0; i < ylen; ++i) {
		y[i] = g1[i] * x[i1[i]] + g2[i] * x[i2[i]];
	}
}

void interpolator::update()
{
	float temp;
	for (int i = 0; i < N; ++i) {
		temp = ind[i];

		temp = floorf(ind[i]);
		i1[i] = (int)temp;
		i2[i] = (int)ceilf(ind[i]);
		g2[i] = ind[i] - temp;
		g1[i] = 1.0f - g2[i];
	}
}

//==================================================

upsampler::upsampler(int rate) : R(rate)
{
	clear();
}
upsampler::~upsampler() {}

void upsampler::setRate(int newRate)
{
	R = newRate;
	clear();
}
int upsampler::getRate() const
{
	return R;
}

void upsampler::clear()
{
	r = R - 1;
}

void upsampler::step(float x, float* y, int ylen)
{
	vset(y, 0.0f, ylen);
	bool xused = false;
	for (int i = 0; i < ylen; ++i) {
		if (++r >= R) {
			r = 0;
			if (!xused) {
				y[i] = R * x;
				xused = true;
			}
		}
	}
}
void upsampler::stepBlock(const float* x, float* y, int xlen, int ylen)
{
	vset(y, 0.0f, ylen);
	int j = 0;
	for (int i = 0; i < ylen; ++i) {
		if (++r >= R) {
			r = 0;
			if (j < xlen) {
				y[i] = R * x[j++];
			}
		}
	}
}

downsampler::downsampler(int rate) : R(rate)
{
	clear();
}
downsampler::~downsampler() {}

void downsampler::setRate(int newRate)
{
	R = newRate;
	clear();
}
int downsampler::getRate() const
{
	return R;
}

void downsampler::clear()
{
	r = R - 1;
}

int downsampler::step(float x, float* y)
{
	int ylen = 0;
	if (++r >= R) {
		r = 0;
		*y = x;
		ylen = 1;
	}
	return ylen;
}
int downsampler::stepBlock(const float* x, float* y, int xlen)
{
	int ylen = 0;
	for (int i = 0; i < xlen; ++i) {
		if (++r >= R) {
			r = 0;
			y[ylen++] = x[i];
		}
	}
	return ylen;
}

//==================================================
/*
waveverb::waveverb(int slew, int maxdel) :
	M(16), N(maxdel + 1),
	alpha(16, slew), delay(16), buffer((maxdel + 1) * 16),
	pind((maxdel + 1) * 16), gind(16, (maxdel + 1) * 16),
	scratch(16), scratch2(4)
{
	vset(alpha.target(), 0.0f, M);
	alpha.converge(); //initialize the alpha parameters
}
waveverb::~waveverb() {}

//Set the -20 dB decay point in samples
void waveverb::setDecay(float samples, bool converge)
{
	decay = samples; //save a reference value
	float* setalpha = alpha.target();
	for (int m = 0; m < M; ++m) { //precompute and slew the sub-parameters
		setalpha[m] = expf(-(constants.common2nats * delay[m]) / (decay + constants.eps)) / 4.0f;
	}
	if (converge) {
		alpha.converge();
	}
}

//Retrieves the M-length waveguide delay pointer for direct adjustment
int* waveverb::getDelayPointer()
{
	return delay;
}

void waveverb::clear()
{
	//clear buffer
	vset(buffer.ptr(), 0.0f, buffer.size());
	
	//reset indices
	pind = 0;
	for (int m = 0; m < M; ++m) {
		gind[m] = m - M * delay[m];
	}
}

float waveverb::step(float x)
{
	alpha.slew();

	for (int m = 0; m < M; ++m) {
		scratch[m] = alpha[m] * buffer[gind[m]]; //voice output with leak and compensation
		gind[m] += M; //step over other voices
	}

	//collect group sums with power adjustment
	const float* p1 = scratch;
	for (int n = 0; n < 4; ++n) {
		scratch2[n] = sum(p1, 4);
		p1 += 4;
	}
	float y = sum(scratch2.ptr(), 4) + x / M; //output sum

	//apply feedback
	p1 = scratch;
	const float* p2 = &(scratch[4]);
	const float* p3 = &(scratch[8]);
	const float* p4 = &(scratch[12]);
	float* p = &(buffer[pind]);
	for (int m = 0; m < 4; ++m) {
		*p++ = y - 2.0f * (scratch2[0] - p1[m] + p2[m] + p3[m] + p4[m]);
	}
	for (int m = 0; m < 4; ++m) {
		*p++ = y - 2.0f * (scratch2[1] + p1[m] - p2[m] + p3[m] + p4[m]);
	}
	for (int m = 0; m < 4; ++m) {
		*p++ = y - 2.0f * (scratch2[2] + p1[m] + p2[m] - p3[m] + p4[m]);
	}
	for (int m = 0; m < 4; ++m) {
		*p++ = y - 2.0f * (scratch2[3] + p1[m] + p2[m] + p3[m] - p4[m]);
	}

	pind += M; //step over other voices

	return 4.0f * (y - x / M);
}
*/
waveverb::waveverb(int slew, int maxbsize, int nominalRate) :
	M(16), N(maxdel0 * nominalRate + maxbsize), maxlen(maxbsize), rate(nominalRate),
	alpha(M, slew), delay(M), fscr(maxbsize), dryscr(maxbsize)
{
	//allocate memory
	for (int m = 0; m < M; ++m) {
		buffer[m].reset(N);
		xscr[m].reset(maxlen);
	}
	for (int m = 0; m < 4; ++m) {
		sscr[m].reset(maxlen);
	}

	setSeed(0); //initialize the delays
	setDecay(0.0f, true); //initialize the alpha parameters
	clear();
}
waveverb::~waveverb() {}

const int waveverb::del0[16][16] = { {765, 918, 911, 4252, 2657, 2605, 1279, 3676, 4185, 3086, 1498, 3650, 2822, 1510, 1319, 2509},
			{658, 2319, 2001, 562, 1019, 3302, 2787, 4730, 2156, 1869, 2033, 2204, 1826, 2661, 2894, 3589},
			{626, 1139, 520, 3404, 4072, 3873, 2662, 3048, 911, 1369, 3476, 1775, 1842, 1077, 2391, 2921},
			{695, 1806, 3468, 2377, 4285, 4105, 3105, 2167, 1426, 2094, 4696, 2698, 3676, 2320, 1021, 3481},
			{848, 2244, 1853, 3087, 1585, 3158, 2714, 2634, 951, 2172, 2054, 4771, 4531, 3480, 4384, 3026},
			{600, 1100, 3737, 4287, 3202, 2534, 3552, 669, 2045, 3292, 4665, 1420, 922, 1932, 1041, 1728},
			{652, 1177, 1624, 2521, 621, 891, 1772, 1107, 4097, 3766, 2286, 2129, 2763, 4543, 3740, 4114},
			{891, 2873, 3696, 2943, 928, 2051, 4046, 3693, 2599, 500, 4463, 1020, 1124, 1257, 1665, 2392},
			{597, 1238, 2566, 2532, 3683, 3748, 1090, 2985, 1944, 3835, 4782, 1386, 689, 1803, 767, 974},
			{850, 3593, 2633, 2816, 2704, 909, 4355, 1090, 1018, 3250, 1976, 4794, 3246, 3566, 3127, 3529},
			{796, 3584, 2001, 3705, 1205, 998, 583, 1810, 1076, 887, 4414, 4068, 1288, 3051, 4612, 2926},
			{713, 2258, 1059, 550, 2278, 4765, 3469, 2852, 4309, 1801, 3857, 2948, 3900, 1080, 1547, 2410},
			{765, 3821, 2875, 1079, 1297, 2087, 869, 4062, 2008, 3896, 1517, 3645, 4223, 3648, 2714, 4544},
			{825, 582, 2310, 3429, 2954, 1641, 2566, 1487, 2220, 1931, 3088, 4463, 1040, 4506, 521, 2839},
			{902, 836, 4450, 3346, 1724, 2551, 1755, 2299, 3037, 1026, 3956, 1171, 3570, 1814, 4294, 1343},
			{841, 1172, 3050, 4168, 1427, 3515, 4788, 4488, 682, 2181, 2366, 4678, 1032, 3138, 773, 1013} };

//Set the -20 dB decay point in samples
void waveverb::setDecay(float samples, bool converge)
{
	decay = samples; //save a reference value
	float* setalpha = alpha.target();
	for (int m = 0; m < M; ++m) { //precompute and slew the sub-parameters
		setalpha[m] = expf(-(constants.common2nats * delay[m]) / (decay + constants.eps)) / 4.0f;
	}
	if (converge) {
		alpha.converge();
	}
}

//Retrieves the M-length waveguide delay pointer for direct adjustment
int* waveverb::getDelayPointer()
{
	return delay;
}

//Set the seed value on [0, 15], defining the waveverb's timbre
void waveverb::setSeed(unsigned int newSeed)
{
	for (int m = 0; m < M; ++m) {
		delay[m] = del0[newSeed][m] * rate;
	}
	setDecay(decay, true); //on seed change, automatically refresh
}

void waveverb::clear()
{
	for (int m = 0; m < M; ++m) {
		buffer[m].clear();
	}
}
void waveverb::stepBlock(const float* x, float* y, int len)
{
	//preserve dry signal in case x == y
	vmult(x, 1.0f / M, dryscr.ptr(), len);

	//get delayed signals from each buffer first
	for (int m = 0; m < M; ++m) {
		buffer[m].get(xscr[m], len, delay[m] - len);
	}

	//apply feedback gain
	if (alpha.check()) { //slew sample by sample
		for (int l = 0; l < len; ++l) {
			alpha.slew();
			for (int m = 0; m < M; ++m) {
				xscr[m][l] *= alpha[m];
			}
		}
	}
	else { //no slew
		for (int m = 0; m < M; ++m) {
			vmult(xscr[m].ptr(), alpha[m], xscr[m].ptr(), len);
		}
	}

	//collect group sums with power adjustment
	vadd(xscr[0].ptr(), xscr[1].ptr(), sscr[0].ptr(), len); //group 1
	vadd(sscr[0].ptr(), xscr[2].ptr(), sscr[0].ptr(), len);
	vadd(sscr[0].ptr(), xscr[3].ptr(), sscr[0].ptr(), len);
	vadd(xscr[4].ptr(), xscr[5].ptr(), sscr[1].ptr(), len); //group 2
	vadd(sscr[1].ptr(), xscr[6].ptr(), sscr[1].ptr(), len);
	vadd(sscr[1].ptr(), xscr[7].ptr(), sscr[1].ptr(), len);
	vadd(xscr[8].ptr(), xscr[9].ptr(), sscr[2].ptr(), len); //group 3
	vadd(sscr[2].ptr(), xscr[10].ptr(), sscr[2].ptr(), len);
	vadd(sscr[2].ptr(), xscr[11].ptr(), sscr[2].ptr(), len);
	vadd(xscr[12].ptr(), xscr[13].ptr(), sscr[3].ptr(), len); //group 4
	vadd(sscr[3].ptr(), xscr[14].ptr(), sscr[3].ptr(), len);
	vadd(sscr[3].ptr(), xscr[15].ptr(), sscr[3].ptr(), len);
	vadd(sscr[0].ptr(), sscr[1].ptr(), y, len); //total
	vadd(y, sscr[2].ptr(), y, len);
	vadd(y, sscr[3].ptr(), y, len);
	vadd(y, dryscr.ptr(), y, len);

	//generate and put feedback signals to the buffer
	vmult(y, 0.5f, y, len);
	vsub(y, sscr[0].ptr(), sscr[0].ptr(), len); //group 1
	for (int m = 0; m < 4; ++m) {
		vadd(sscr[0].ptr(), xscr[m].ptr(), fscr.ptr(), len);
		vsub(fscr.ptr(), xscr[m + 4].ptr(), fscr.ptr(), len);
		vsub(fscr.ptr(), xscr[m + 8].ptr(), fscr.ptr(), len);
		vsub(fscr.ptr(), xscr[m + 12].ptr(), fscr.ptr(), len);
		vmult(fscr.ptr(), 2.0f, fscr.ptr(), len);
		buffer[m].put(fscr.ptr(), len);
	}
	vsub(y, sscr[1].ptr(), sscr[1].ptr(), len); //group 2
	for (int m = 0; m < 4; ++m) {
		vsub(sscr[1].ptr(), xscr[m].ptr(), fscr.ptr(), len);
		vadd(fscr.ptr(), xscr[m + 4].ptr(), fscr.ptr(), len);
		vsub(fscr.ptr(), xscr[m + 8].ptr(), fscr.ptr(), len);
		vsub(fscr.ptr(), xscr[m + 12].ptr(), fscr.ptr(), len);
		vmult(fscr.ptr(), 2.0f, fscr.ptr(), len);
		buffer[m + 4].put(fscr.ptr(), len);
	}
	vsub(y, sscr[2].ptr(), sscr[2].ptr(), len); //group 3
	for (int m = 0; m < 4; ++m) {
		vsub(sscr[2].ptr(), xscr[m].ptr(), fscr.ptr(), len);
		vsub(fscr.ptr(), xscr[m + 4].ptr(), fscr.ptr(), len);
		vadd(fscr.ptr(), xscr[m + 8].ptr(), fscr.ptr(), len);
		vsub(fscr.ptr(), xscr[m + 12].ptr(), fscr.ptr(), len);
		vmult(fscr.ptr(), 2.0f, fscr.ptr(), len);
		buffer[m + 8].put(fscr.ptr(), len);
	}
	vsub(y, sscr[3].ptr(), sscr[3].ptr(), len); //group 4
	for (int m = 0; m < 4; ++m) {
		vsub(sscr[3].ptr(), xscr[m].ptr(), fscr.ptr(), len);
		vsub(fscr.ptr(), xscr[m + 4].ptr(), fscr.ptr(), len);
		vsub(fscr.ptr(), xscr[m + 8].ptr(), fscr.ptr(), len);
		vadd(fscr.ptr(), xscr[m + 12].ptr(), fscr.ptr(), len);
		vmult(fscr.ptr(), 2.0f, fscr.ptr(), len);
		buffer[m + 12].put(fscr.ptr(), len);
	}

	//remove dry signal and gain adjust output
	vmultaccum(dryscr.ptr(), -0.5f, y, len);
	vmult(y, 8.0f, y, len);
}

/*
miniverb::miniverb(int slew, int maxdel) :
	M(9), N(maxdel + 1),
	alpha(M, slew), delay(M), buffer(N * M),
	pind(N * M), gind(M, N * M),
	scratch(M), scratch2(3)
{
	vset(alpha.target(), 0.0f, M);
	alpha.converge(); //initialize the alpha parameters
}
miniverb::~miniverb() {}

//Set the -20 dB decay point in samples
void miniverb::setDecay(float samples, bool converge)
{
	decay = samples; //save a reference value
	float* setalpha = alpha.target();
	for (int m = 0; m < M; ++m) { //precompute and slew the sub-parameters
		setalpha[m] = expf(-(constants.common2nats * delay[m]) / (decay + constants.eps)) / 3.0f;
	}
	if (converge) {
		alpha.converge();
	}
}

//Retrieves the M-length waveguide delay pointer for direct adjustment
int* miniverb::getDelayPointer()
{
	return delay;
}

void miniverb::clear()
{
	//clear buffer
	vset(buffer.ptr(), 0.0f, buffer.size());

	//reset indices
	pind = 0;
	for (int m = 0; m < M; ++m) {
		gind[m] = m - M * delay[m];
	}
}

float miniverb::step(float x)
{
	alpha.slew();

	for (int m = 0; m < M; ++m) {
		scratch[m] = alpha[m] * buffer[gind[m]]; //voice output with leak and compensation
		gind[m] += M; //step over other voices
	}

	//collect group sums with power adjustment
	const float* p1 = scratch;
	for (int n = 0; n < 3; ++n) {
		scratch2[n] = sum(p1, 3);
		p1 += 3;
	}
	float y = sum(scratch2.ptr(), 3) + x / M; //output sum

	//apply feedback
	p1 = scratch;
	const float* p2 = &(scratch[3]);
	const float* p3 = &(scratch[6]);
	float* p = &(buffer[pind]);
	for (int m = 0; m < 3; ++m) {
		*p++ = (4.0f / 3.0f) * y - 2.0f * (scratch2[0] - 0.5f * p1[m] + p2[m] + p3[m]);
	}
	for (int m = 0; m < 3; ++m) {
		*p++ = (4.0f / 3.0f) * y - 2.0f * (scratch2[1] + p1[m] - 0.5f * p2[m] + p3[m]);
	}
	for (int m = 0; m < 3; ++m) {
		*p++ = (4.0f / 3.0f) * y - 2.0f * (scratch2[2] + p1[m] + p2[m] - 0.5f * p3[m]);
	}

	pind += M; //step over other voices

	return 3.0f * (y - x / M);
}
*/
miniverb::miniverb(int slew, int maxbsize, int nominalRate) :
	M(9), N(maxdel0 * nominalRate + maxbsize), maxlen(maxbsize), rate(nominalRate),
	alpha(M, slew), delay(M), fscr(maxbsize), dryscr(maxbsize)
{
	//allocate memory
	for (int m = 0; m < M; ++m) {
		buffer[m].reset(N);
		xscr[m].reset(maxlen);
	}
	for (int m = 0; m < 3; ++m) {
		sscr[m].reset(maxlen);
	}

	setSeed(0); //initialize the delays
	setDecay(0.0f, true); //initialize the alpha parameters
	clear();
}
miniverb::~miniverb() {}

const int miniverb::del0[16][9] = {
		{ 781, 1508, 1651, 671, 1131, 4479, 3502, 2357, 2647 },
		{ 637, 587, 884, 1781, 2148, 3337, 3870, 2426, 1315 },
		{ 523, 3750, 3562, 1198, 3072, 1106, 3952, 4563, 3371 },
		{ 651, 1246, 3199, 1643, 862, 3659, 541, 2265, 1887 },
		{ 495, 946, 3643, 3499,4043, 1499, 2790, 4114, 2261 },
		{ 772, 1991, 4033, 3170, 2815, 2248, 3525, 1821, 4454 },
		{ 539, 2922, 1947, 3084, 4506, 2463, 1444, 3630, 2977 },
		{ 543, 2716, 1543, 1396, 3784, 2802, 817, 4477, 1968 },
		{ 707, 2603, 4612, 484, 3480, 1320, 1634, 4321, 1897 },
		{ 888, 4428, 561, 4361, 2769, 3782, 1995, 1303, 4675 },
		{ 631, 863, 2862, 2592, 3110, 1445, 3266, 2737, 1113 },
		{ 572, 2616, 3818, 2268, 4173, 2927, 2390, 4540, 935 },
		{ 611, 1966, 2135, 1517, 2638, 2516, 850, 3981, 4626 },
		{ 604, 893, 3333, 2598, 2480, 2873, 2925, 1516, 4357 },
		{ 527, 3543, 3158, 3829, 2350, 1445, 3325, 1090, 1362 },
		{ 790, 3929, 1233, 3639, 2408, 1873, 3068, 3279, 2983 } };

//Set the -20 dB decay point in samples
void miniverb::setDecay(float samples, bool converge)
{
	decay = samples; //save a reference value
	float* setalpha = alpha.target();
	for (int m = 0; m < M; ++m) { //precompute and slew the sub-parameters
		setalpha[m] = expf(-(constants.common2nats * delay[m]) / (decay + constants.eps)) / 3.0f;
	}
	if (converge) {
		alpha.converge();
	}
}

//Retrieves the M-length waveguide delay pointer for direct adjustment
int* miniverb::getDelayPointer()
{
	return delay;
}

//Set the seed value on [0, 15], defining the miniverb's timbre
void miniverb::setSeed(unsigned int newSeed)
{
	for (int m = 0; m < M; ++m) {
		delay[m] = del0[newSeed][m] * rate;
	}
	setDecay(decay, true); //on seed change, automatically refresh
}

void miniverb::clear()
{
	for (int m = 0; m < M; ++m) {
		buffer[m].clear();
	}
}
void miniverb::stepBlock(const float* x, float* y, int len)
{
	//preserve dry signal in case x == y
	vcopy(x, dryscr.ptr(), len);

	//get delayed signals from each buffer first
	for (int m = 0; m < M; ++m) {
		buffer[m].get(xscr[m], len, delay[m] - len);
	}

	//apply feedback gain
	if (alpha.check()) { //slew sample by sample
		for (int l = 0; l < len; ++l) {
			alpha.slew();
			for (int m = 0; m < M; ++m) {
				xscr[m][l] *= alpha[m];
			}
		}
	}
	else { //no slew
		for (int m = 0; m < M; ++m) {
			vmult(xscr[m].ptr(), alpha[m], xscr[m].ptr(), len);
		}
	}

	//collect group sums with power adjustment
	vadd(xscr[0].ptr(), xscr[1].ptr(), sscr[0].ptr(), len); //group 1
	vadd(sscr[0].ptr(), xscr[2].ptr(), sscr[0].ptr(), len);
	vadd(xscr[3].ptr(), xscr[4].ptr(), sscr[1].ptr(), len); //group 2
	vadd(sscr[1].ptr(), xscr[5].ptr(), sscr[1].ptr(), len);
	vadd(xscr[6].ptr(), xscr[7].ptr(), sscr[2].ptr(), len); //group 3
	vadd(sscr[2].ptr(), xscr[8].ptr(), sscr[2].ptr(), len);
	vadd(sscr[0].ptr(), sscr[1].ptr(), y, len); //total
	vadd(y, sscr[2].ptr(), y, len);
	vmultaccum(dryscr.ptr(), 1.0f / M, y, len);

	//generate and put feedback signals to the buffer
	vmult(y, 2.0f / 3.0f, y, len);
	vsub(y, sscr[0].ptr(), sscr[0].ptr(), len); //group 1
	for (int m = 0; m < 3; ++m) {
		vcopy(sscr[0].ptr(), fscr.ptr(), len);
		vmultaccum(xscr[m].ptr(), 0.5f, fscr.ptr(), len);
		vsub(fscr.ptr(), xscr[m + 3].ptr(), fscr.ptr(), len);
		vsub(fscr.ptr(), xscr[m + 6].ptr(), fscr.ptr(), len);
		vmult(fscr.ptr(), 2.0f, fscr.ptr(), len);
		buffer[m].put(fscr.ptr(), len);
	}
	vsub(y, sscr[1].ptr(), sscr[1].ptr(), len); //group 2
	for (int m = 0; m < 3; ++m) {
		vcopy(sscr[1].ptr(), fscr.ptr(), len);
		vsub(fscr.ptr(), xscr[m].ptr(), fscr.ptr(), len);
		vmultaccum(xscr[m + 3].ptr(), 0.5f, fscr.ptr(), len);
		vsub(fscr.ptr(), xscr[m + 6].ptr(), fscr.ptr(), len);
		vmult(fscr.ptr(), 2.0f, fscr.ptr(), len);
		buffer[m + 3].put(fscr.ptr(), len);
	}
	vsub(y, sscr[2].ptr(), sscr[2].ptr(), len); //group 3
	for (int m = 0; m < 3; ++m) {
		vcopy(sscr[2].ptr(), fscr.ptr(), len);
		vsub(fscr.ptr(), xscr[m].ptr(), fscr.ptr(), len);
		vsub(fscr.ptr(), xscr[m + 3].ptr(), fscr.ptr(), len);
		vmultaccum(xscr[m + 6].ptr(), 0.5f, fscr.ptr(), len);
		vmult(fscr.ptr(), 2.0f, fscr.ptr(), len);
		buffer[m + 6].put(fscr.ptr(), len);
	}

	//remove dry signal and gain adjust output
	vmult(y, 9.0f / 2.0f, y, len);
	vmultaccum(dryscr.ptr(), -1.0f / 3.0f, y, len);
}

diffuser::diffuser(int slew, int maxbsize, int nominalRate) :
	rate(nominalRate),
	core(slew, maxbsize, nominalRate),
	level(slew),
	scratch(maxbsize)
{
	setSeed(0); //init
	setDecay(0.0f, true);
}
diffuser::~diffuser() {}

const float diffuser::g[16] =
{ 2.9f, 3.4f, 1.9f, 3.8f, 2.1f, 1.7f, 1.8f, 2.5f, 2.5f, 2.0f, 2.9f, 1.9f, 2.4f, 2.2f, 2.7f, 1.9f };

//Set the -20 dB decay point in samples
void diffuser::setDecay(float samples, bool converge)
{
	core.setDecay(samples, converge);
	tdecay = samples / (48000.0f * rate); //favor power accuracy over timing accuracy
	level.setLevel(getNormGain(), converge);
}

//Set the seed value on [0, 15], defining the miniverb's timbre
void diffuser::setSeed(unsigned int newSeed)
{
	seed = newSeed;
	core.setSeed(seed);
	level.setLevel(getNormGain(), true); //always converge on seed update
}

void diffuser::clear()
{
	core.clear();
}
void diffuser::stepBlock(const float* x, float* y, int len)
{
	vcopy(x, scratch.ptr(), len); //preserve dry in case x == y
	core.stepBlock(x, y, len); //generate wet path
	vmultaccum(scratch.ptr(), drymix, y, len); //mix with low level dry
	level.stepBlock(y, y, len); //gain adjust for power norm
}

float diffuser::getNormGain() const
{
	return sqrtf(1.0f / (powf(tdecay, r) * g[seed] + drymix * drymix));
}

//==================================================

corrbuffer::corrbuffer(int len, int chunk, int dsr, int maxbsize) :
	Nu(len), Nd(len/dsr), Md(chunk/dsr), DSR(dsr), maxlen(maxbsize),
	Rlen(Nd - Md + 1), K(getKForNpts(Nd)),
	F(Nd), D(dsr),
	T(constants.eps),
	ave(1, (chunk/dsr)*100.0, chunk/dsr, maxbsize/dsr),
	buffd(Nd), bpowd(Nd - Md + 1),
	fscr1(K), fscr2(K), tscr(max(Nd, maxbsize))
{
	sofcasc<4>::coefs C; //design anti-aliasing filter
	aa.setCoefs(filt::antialias8(C, dsr), true);
	ave.setLength(Md, true);
	clear(); //pre-clear buffer state
}
corrbuffer::~corrbuffer() {}

int corrbuffer::size() const
{
	return Nu;
}
int corrbuffer::dsize() const
{
	return Nd;
}
int corrbuffer::Rsize() const
{
	return Rlen;
}

//Power threshold for division stability
void corrbuffer::setThresh(float powerThreshold)
{
	T = powerThreshold;
	clear();
}

//Correlation samples less than this are ignored by peak finding
void corrbuffer::setMinPeriod(int samples)
{
	minperd = samples / DSR;
}

void corrbuffer::clear()
{
	buffd.clear();
	bpowd.clear(T);
	D.clear();
	ave.clear();
}

void corrbuffer::put(const float* x, int len)
{
	aa.stepBlock(x, tscr, len);
	len = D.stepBlock(tscr, tscr, len);
	buffd.put(tscr, len); //put downsampled buffer

	vmult(tscr.ptr(), tscr.ptr(), tscr.ptr(), len);
	vmult(tscr.ptr(), (float)Md, tscr.ptr(), len); //convert ave -> sum
	ave.stepBlock(tscr, tscr, len);
	vmax(tscr.ptr(), T, tscr.ptr(), len); //guarantee norm stability
	bpowd.put(tscr, len); //update chunk power sequence vector
}

//Computes the downsampled normalized autocorrelation sequence and copies to R
void corrbuffer::corr(float* R)
{
	vrevcopy(buffd.get(Nd), tscr.ptr(), Nd); //instead of freq domain conj
	F.rfwd(tscr, fscr1, Nd); //flipped full buffer
	F.rfwd(buffd.get(Md), fscr2, Md); //front chunk

	vmult(fscr1.ptr(), fscr2.ptr(), fscr1.ptr(), K); //perform autocorrelation
	F.rinv(fscr1, tscr, Nd);

	vrevcopy(bpowd.get(Rlen), R, Rlen); //normalize
	vmult(R, R[0], R, Rlen);
	vsqrt(R, R, Rlen);
	vdiv(tscr.ptr() + (Md - 1), R, R, Rlen);
}

//Searches R over its valid range for the largest positive value
float corrbuffer::corrpeak(const float* R, int& delay) const
{
	float Rpeak = vfindmax(R + minperd, Rlen - minperd, delay);
	delay += minperd;
	return Rpeak;
}

downshift::downshift(float sampleStep, int blendLength, const corrbuffer& Robj) :
	Nu(Robj.Nu), DSR(Robj.DSR), maxlen(Robj.maxlen), Rlen(Robj.Rlen),
	ddel(1.0f - sampleStep), blendlen(blendLength),
	maxdeld(min((int)floorf((Nu - maxlen - blendLength * ddel) / DSR), Rlen - 1) - 1),
	maxdel(maxdeld * DSR),
	buff(Robj.Nu), blend(blendlen, false),
	scr(maxlen)
{
	//design the post-filter
	sofcasc<2>::coefs C;
	lpf.setCoefs(filt::lowpass4(C, 0.9f * sampleStep, 0.707f), true);
}
downshift::~downshift() {}

//Sets the absolute correlation threshold that flags a transient and snaps to front
void downshift::setSnapThresh(float thresh)
{
	Tsnap = thresh;
}
//Sets the relative correlation threshold that flags a period to blend across
void downshift::setBlendThresh(float thresh)
{
	Tblend = thresh;
}

void downshift::clear()
{
	buff.clear();
	lpf.clear();
	blend.converge();
	del1 = 0.0f;
	del2 = 0.0f;
	intransient = false;
}

void downshift::stepBlock(const float* x, float* y, int len, const float* R, float Rmax, int Rmaxdel)
{
	//put new data to the buffer
	buff.put(x, len);

	int nbs = -1; //countdown to a new blend, -1 for none
	float del1bs = 0.0f; //del1 assignment on new blend
	float del2bs = 0.0f; //del2 assignment on new blend

	if (!(blend.check())) { //if no active blend, check to see if we need to start one
		if (Rmax < Tsnap) { //transient detected, flag
			intransient = true;
		}
		if(intransient){
			del1bs = 0.0f; //snap to front
			del2bs = del1;
			nbs = 0; //immediately
			if (Rmax >= Tsnap) {
				intransient = false; //clear after snap to a periodic region
			}
		}
		else {
			//search upcoming values for a valid period blend point
			int curmindeld = (int)floorf(del1 / DSR);
			int curmaxdeld = min((int)ceilf((del1 + (len - 1) * ddel) / DSR), maxdeld);
			float curbestdeld = -1.0f; //-1 to indicate no period blend
			for (int i = curmindeld; i < curmaxdeld; ++i) {
				if (curbestdeld < 0.0f) { //no blend found yet
					if (R[i] >= (Rmax * Tblend)) { //adequate blend correlation
						if ((R[i] > R[i - 1]) && (R[i] > R[i + 1])) { //local maxima
							curbestdeld = (float)i + parabsolve(R[i - 1], R[i], R[i + 1]); //assign
						}
					}
				}
			}
			if (curbestdeld >= 0.0f) { //probably found a blend, try to snap
				nbs = (int)ceilf((curbestdeld * DSR - del1) / ddel);
				if (nbs < len) { //avoid small chance we picked something out of reach of this block
					del2bs = del1 + nbs * ddel;
					del1bs = del2bs - curbestdeld * DSR; //guaranteed small >= 0
				}
				else {
					nbs = -1; //skip if we picked something that wouldn't start blending this block
				}
			}
			if (nbs < 0) { //did not find a blend
				nbs = (int)ceilf((maxdel - del1) / ddel); //samples remaining until we run out of buffer
				if (nbs < len) { //if we will reach the end of the buffer this block
					del2bs = del1 + nbs * ddel;
					del1bs = del2bs - maxdel; //guaranteed small >= 0
				}
				else {
					nbs = -1; //no danger of hitting the end of the buffer
				}
			}
		}
	}

	if (blend.check()) { //if already blending, use two voices over the whole range and manage gains
		buff.getNN(y, len, del1, ddel); //primary voice
		buff.getNN(scr, len, del2, ddel); //secondary voice
		del1 += (ddel * len);
		del2 += (ddel * len);
		for (int n = 0; n < len; ++n) { //apply blending gains, forming output
			blend.fade();
			y[n] *= blend.gain();
			y[n] += scr[n] * (1.0f - blend.gain());
		}
	}
	else {
		if (nbs >= 0) { //if starting a new blend, split the block into one voice | two voices
			buff.getNN(y, nbs, del1, ddel, len, 0); //one voice portion

			blend.target(!(blend.current())); //kick off new blend
			blend.swap();
			del1 = del1bs;
			del2 = del2bs;

			buff.getNN(y + nbs, len - nbs, del1, ddel, len, nbs); //primary voice
			buff.getNN(scr.ptr() + nbs, len - nbs, del2, ddel, len, nbs); //secondary voice
			del1 += (ddel * (len - nbs));
			del2 += (ddel * (len - nbs));
			for (int n = nbs; n < len; ++n) { //apply blending gains, forming output
				blend.fade();
				y[n] *= blend.gain();
				y[n] += scr[n] * (1.0f - blend.gain());
			}
		}
		else { //no active or new blends, free run one voice
			buff.getNN(y, len, del1, ddel);
			del1 += (ddel * len);
		}
	}

	lpf.stepBlock(y, y, len); //post-filter
}

upshift::upshift(float sampleStep, int blendLength, const corrbuffer& Robj) :
	Nu(Robj.Nu), DSR(Robj.DSR), maxlen(Robj.maxlen), Rlen(Robj.Rlen),
	ddel(1.0f - sampleStep), blendlen(blendLength),
	mindeld((int)ceilf((ceilf(-blendLength * ddel / maxlen) * maxlen - blendLength * ddel) / DSR) + 1), //ddel <= 0
	mindel(mindeld * DSR),
	buff(Robj.Nu), blend(blendlen, false),
	scr(maxlen)
{
	//design the pre-filter
	sofcasc<2>::coefs C;
	lpf.setCoefs(filt::lowpass4(C, 0.9f / sampleStep, 0.707f), true);
}
upshift::~upshift() {}

//Sets the absolute correlation threshold that flags a transient and snaps to front
void upshift::setSnapThresh(float thresh)
{
	Tsnap = thresh;
}
//Sets the relative correlation threshold that flags a period to blend across
void upshift::setBlendThresh(float thresh)
{
	Tblend = thresh;
}
//Minimum sample shift corresponding roughly to highest pitch period
void upshift::setMinPeriod(int samples)
{
	minperd = max(samples / DSR, 1); //floor to downsampling rate
	minper = minperd * DSR;
}

void upshift::clear()
{
	buff.clear();
	lpf.clear();
	blend.converge();
	del1 = (float)mindel;
	del2 = 0.0f;
}

void upshift::stepBlock(const float* x, float* y, int len, const float* R, float Rmax, int Rmaxdel)
{
	lpf.stepBlock(x, y, len); //pre-filter

	//put new data to the buffer
	buff.put(y, len);

	bool newblend = false; //flags a new blend
	float del1bs = 0.0f; //del1 assignment on new blend
	float del2bs = 0.0f; //del2 assignment on new blend

	if (!(blend.check())) { //if no active blend, check to see if we need to start one
		if (Rmax < Tsnap) { //transient detected,
			del1bs = (float)mindel; //snap to front
			del2bs = del1;
			newblend = true;
		}
		else {
			if (del1 <= (float)mindel) { //nearing the end of the buffer,
				int curmindeld = (int)ceilf(mindeld - del1 / DSR); //search for furthest backward acceptable blend point
				int curmaxdeld = min((int)floorf((Nu - maxlen + 1 - del1) / DSR), Rlen - 1);
				bool isdone = false;
				int i = curmaxdeld;
				float curbestdeld = (float)curmindeld;
				float curbestR = 0.0f;
				while (!isdone) {
					if (R[i] >= (Rmax * Tblend)) { //adequate blend correlation
						if ((R[i] > R[i - 1]) && (R[i] > R[i + 1])) { //local maxima
							curbestdeld = (float)i + parabsolve(R[i - 1], R[i], R[i + 1]); //assign
							isdone = true; //and escape
						}
						else if(R[i] > curbestR){ //not local maxima, but still keep track of best found
							curbestdeld = (float)i; //assign
							curbestR = R[i]; //but don't escape, wait for something better
						}
					}
					else if(R[i] > curbestR){ //inadequate correlation, but still keep track of best found
						if ((R[i] > R[i - 1]) && (R[i] > R[i + 1])) { //local maxima
							curbestdeld = (float)i + parabsolve(R[i - 1], R[i], R[i + 1]); //assign accurately
							curbestR = R[i]; //but don't escape, wait for something better
						}
						else { //not local maxima
							curbestdeld = (float)i; //assign directly
							curbestR = R[i];
						}
					}
					if (--i < curmindeld) {
					//if (++i > curmaxdeld) {
						isdone = true; //escape with our best found R value
					}
				}
				del1bs = del1 + curbestdeld * DSR; //assign the best blend
				del2bs = del1;
				newblend = true;
			}
			//if not near the end of the buffer, no need to blend
		}
	}

	if (newblend) { //kick off new blend
		blend.target(!(blend.current()));
		blend.swap();
		del1 = del1bs;
		del2 = del2bs;
	}

	if (blend.check()) { //if already blending, use two voices and manage gains
		buff.getNN(y, len, del1, ddel); //primary voice
		buff.getNN(scr, len, del2, ddel); //secondary voice
		del1 += (ddel * len);
		del2 += (ddel * len);
		for (int n = 0; n < len; ++n) { //apply blending gains, forming output
			blend.fade();
			y[n] *= blend.gain();
			y[n] += scr[n] * (1.0f - blend.gain());
		}
	}
	else { //no active blend, free run one voice
		buff.getNN(y, len, del1, ddel);
		del1 += (ddel * len);
	}
}

} //namespace cephean

