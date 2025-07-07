
//--------------------------------------------------
// Defines classes and functions helpful for
// creating LV2 plugins and other very high level
// general purpose stuff
//--------------------------------------------------

#include "cephean-lv2.h"

namespace cephean
{

//==================================================

disableDenormals::disableDenormals()
{
#if ISWINDOWS
	MXCSR = _mm_getcsr();
	_mm_setcsr(MXCSR |= 0x8040); //PUSH
#endif
#if ISARM
	fegetenv(&fenv);
	fesetenv(FE_DFL_DISABLE_SSE_DENORMS_ENV); //PUSH
#endif
}
disableDenormals::~disableDenormals()
{
#if ISWINDOWS
	_mm_setcsr(MXCSR); //POP
#endif
#if ISARM
	fesetenv(&fenv); //POP
#endif
}

//==================================================

plugin::plugin(int numControls, int numInputs, int numOutputs,
	float setFs,
	const char* bundle_path,
	const LV2_Feature* const* features) : Fs(setFs),
	vctrl(numControls), vin(numInputs), vout(numOutputs)
{
	bundle_path; //unused
	features; //unused
}

plugin::~plugin() {}

void plugin::connect_port(uint32_t port, void* data)
{
	for (int i = 0; i < vctrl.size(); ++i) {
		if (vctrl[i]->getIndex() == port) {
			vctrl[i]->connect(data);
		}
	}
	for (int i = 0; i < vin.size(); ++i) {
		if (vin[i]->getIndex() == port) {
			vin[i]->connect(data);
		}
	}
	for (int i = 0; i < vout.size(); ++i) {
		if (vout[i]->getIndex() == port) {
			vout[i]->connect(data);
		}
	}
}

void plugin::commit()
{
	for (int i = 0; i < vctrl.size(); ++i) {
		vctrl[i]->commit();
	}
}

plugin::control::control(int portIndex) : index(portIndex) {}
plugin::control::~control() {}
int plugin::control::getIndex() const
{
	return index;
}

bool plugin::control::isnew() const
{
	return vnew;
}
void plugin::control::clearnew()
{
	vnew = false;
}

plugin::input::input(plugin* host, int portIndex) : index(portIndex)
{
	host->vin[host->in_count++] = this;
}
plugin::input::~input() {}
int plugin::input::getIndex() const
{
	return index;
}
void plugin::input::connect(void* data)
{
	p = static_cast<const float*>(data);
}
plugin::input::operator const float*() const
{
	return p;
}
const float* plugin::input::get() const
{
	return p;
}

plugin::output::output(plugin* host, int portIndex) : index(portIndex)
{
	host->vout[host->out_count++] = this;
}
plugin::output::~output() {}
int plugin::output::getIndex() const
{
	return index;
}
void plugin::output::connect(void* data)
{
	p = static_cast<float*>(data);
}
plugin::output::operator float* ()
{
	return p;
}
float* plugin::output::get()
{
	return p;
}

void plugin::init() {}
void plugin::deactivate() {}

//==================================================

//Returns 1, 2, or 4 to represet 1x, 2x, and 4x nominal
//sampling rates. Rates < 48k yield 1, and > 192 yield 4,
//otherwise the nearest rate is used:
//	1 - 48k
//	2 - 96k
//	4 - 192k
int getNominalUSR(float Fs)
{
	Fs = bound(Fs, 48000.0f, 192000.0f);
	return (int)roundf(Fs / 48000.0f);
}

//==================================================

//Collection of mathematical constants
const Constants constants;

//==================================================

//Converts a sorted index vector ind of length len into a place vector
//of the same length. Ex. [2, 0, 1] -> [1, 2, 0]. ind == place is UNSAFE
void unsort(const int* ind, int* place, int len)
{
	for (int i = 0; i < len; ++i) {
		place[ind[i]] = i;
	}
}

//==================================================

float boundlinmap(float x, float xmin, float xmax, float ymin, float ymax)
{
	return bound(linmap(x, xmin, xmax, ymin, ymax), ymin, ymax);
}

float boundlogmap(float x, float xmin, float xmax, float ymin, float ymax)
{
	return bound(expf(linmap(x, xmin, xmax, logf(ymin), logf(ymax))), ymin, ymax);
}

//==================================================

float rms(const float* x, int len)
{
	return sqrtf(ssq(x, len) / len);
}
double rms(const double* x, int len)
{
	return sqrt(ssq(x, len) / len);
}
float rssq(const float* x, int len)
{
	return sqrtf(ssq(x, len));
}
double rssq(const double* x, int len)
{
	return sqrt(ssq(x, len));
}

//==================================================

void vmagdB(const float* x, float* y, int len)
{
	for (int i = 0; i < len; ++i) {
		y[i] = 20.0f * log10f(x[i]);
	}
}
void vpowdB(const float* x, float* y, int len)
{
	for (int i = 0; i < len; ++i) {
		y[i] = 10.0f * log10f(x[i]);
	}
}

void vlog(const float* x, float* y, int len)
{
	for (int i = 0; i < len; ++i) {
		y[i] = logf(x[i]);
	}
}
void vexp(const float* x, float* y, int len)
{
	for (int i = 0; i < len; ++i) {
		y[i] = expf(x[i]);
	}
}
void vsqrt(const float* x, float* y, int len)
{
	for (int i = 0; i < len; ++i) {
		y[i] = sqrtf(x[i]);
	}
}

void vexpi(const float* ph, complex<float>* y, int len)
{
	for (int i = 0; i < len; ++i) {
		y[i] = expi(ph[i]);
	}
}

//==================================================

rebuffer::rebuffer(int maxlen) : N(maxlen)
{
	reset();
}
rebuffer::~rebuffer() {}

//Clears internal sample count to 0. Typically called on init or
//some other whole plugin reinitialization.
void rebuffer::reset()
{
	n = 0;
}

//Returns the next application block size when remlen samples
//remain in the I/O buffers. Typical usage:
//	//pointer p holds data, len holds block length
// while(len > 0){
//	curlen = rebuffer.next(len);
//	<DO PROCESSING ON p OVER curlen SAMPLES>
//	p += curlen; //step the data pointer forward
//  len -= curlen; //reduce the remaining samples
// }
int rebuffer::next(int remlen)
{
	int ret = N - n;
	if (ret <= remlen) {
		n = 0;
	}
	else {
		ret = remlen;
		n += remlen;
	}
	return ret;
}

//==================================================

monoalg::monoalg() {}
monoalg::~monoalg() {}

void monoalg::stepBlock(const float* x, float* y, int len)
{
	for (int i = 0; i < len; ++i) {
		y[i] = step(x[i]);
	}
}

stereoalg::stereoalg() {}
stereoalg::~stereoalg() {}

void stereoalg::stepBlock(const float* xL, const float* xR,
	float* yL, float* yR, int len)
{
	for (int i = 0; i < len; ++i) {
		step(xL[i], xR[i], yL[i], yR[i]);
	}
}

multialg::multialg(int maxNch, bool blockTransposed) :
	xscr(maxNch), yscr(maxNch), transposed(blockTransposed)
{
}
multialg::~multialg() {}

void multialg::reallocate(int maxNch, bool blockTransposed)
{
	xscr.reset(maxNch);
	yscr.reset(maxNch);
	transposed = blockTransposed;
}

void multialg::stepBlock(const float* const* x, float* const* y, int chan, int len)
{
	if (transposed) {
		for (int i = 0; i < len; ++i) {
			step(x[i], y[i], chan);
		}
	}
	else {
		for (int i = 0; i < len; ++i) {
			for (int j = 0; j < chan; ++j) {
				xscr[j] = x[j][i];
			}
			step(xscr, yscr, chan);
			for (int j = 0; j < chan; ++j) {
				y[j][i] = yscr[j];
			}
		}
	}
	
}

//==================================================

range::range(float setxmin, float setxmax, float setymin, float setymax, shape setshape)
{
	set(setxmin, setxmax, setymin, setymax, setshape);
}
range::~range() {}

void range::set(float setxmin, float setxmax, float setymin, float setymax, shape setshape)
{
	xmin = setxmin;
	xmax = setxmax;
	rshape = setshape;

	switch (rshape)
	{
	case shape::LOG:
		ymin = logf(setymin);
		ymax = logf(setymax);
		break;

	default: //LIN
		ymin = setymin;
		ymax = setymax;
		break;
	}
}

float range::map(float x) const
{
	x = (x - xmin) / (xmax - xmin);
	switch (rshape)
	{
	case shape::LOG:
		return expf(bound(x * (ymax - ymin) + ymin, ymin, ymax)); //range values already log'd

	default: //LIN
		return bound(x * (ymax - ymin) + ymin, ymin, ymax);
	}
}
float range::unmap(float y) const
{
	float x = 0.0f;
	switch (rshape)
	{
	case shape::LOG:
		x = (logf(y) - ymin) / (ymax - ymin); //range values already log'd
		break;

	default: //LIN
		x = (y - ymin) / (ymax - ymin);
		break;
	}
	return bound(x*(xmax - xmin) + xmin, xmin, xmax);
}


//==================================================

multirange::multirange(int stops) :
	xstop(stops), ystop(stops), r(stops - 1)
{}
multirange::~multirange() {}

//the shape applies from this x to the next one, unused at the list end
void multirange::setStop(int index, float x, float y, shape setshape)
{
	xstop[index] = x;
	ystop[index] = y;

	if (index < r.size()) {
		r[index].set(xstop[index], xstop[index+1],
			ystop[index], ystop[index+1], setshape);
	}
	if (index > 0) {
		r[index-1].set(xstop[index-1], xstop[index],
			ystop[index-1], ystop[index], r[index-1].rshape);
	}
}

float multirange::map(float x) const
{
	for (int i = 0; i < r.size(); ++i) {
		if (x <= r[i].xmax) {
			return r[i].map(x);
		}
	}
	return r[r.size() - 1].map(x);
}


} //namespace cephean

