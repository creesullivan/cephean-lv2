
//--------------------------------------------------
// Defines classes and functions for offline testing
// and debugging
//--------------------------------------------------

#include "cephean-test.h"

namespace cephean
{

//==================================================

wavreader::wavreader(const std::string& filename, double duration, double startTime, double endTime)
{
	fr.shouldLogErrorsToConsole(true);
	fr.load(filename);
	std::cout << "Read from " + filename << std::endl;
	const double Fs = (double)fr.getSampleRate();

	Lfile = fr.getNumSamplesPerChannel();
	lf2 = (endTime < 0.0) ? (Lfile - 1) : std::min((unsigned long int)round(endTime * Fs), Lfile - 1);
	lf1 = std::min((unsigned long int)round(startTime * Fs), lf2);
	lf = lf1;
	
	L = (unsigned long int)round(duration * Fs);
	l = 0;
}
wavreader::~wavreader() {}

//Returns the sampling rate of the read audio file
double wavreader::getSampleRate() const
{
	return (double)fr.getSampleRate();
}

//Fills a mono buffer block of len samples from the file.
//If we have exceeded the duration, pads with 0s.
void wavreader::get(float* buff, int len)
{
	for (int i = 0; i < len; ++i) {
		if (l < L) {
			buff[i] = fr.samples[0][lf++];
			if (lf > lf2) {
				lf = lf1; //loop
			}
			++l;
		}
		else {
			buff[i] = 0.0f;
		}
	}
}

//Fills a stereo buffer block of len samples from the file.
//If we have exceeded the duration, pads with 0s, and if the
//file is mono, duplicates to stereo.
void wavreader::get(float* buffL, float* buffR, int len)
{
	for (int i = 0; i < len; ++i) {
		if (l < L) {
			buffL[i] = fr.samples[0][lf];
			if (fr.getNumChannels() >= 2) {
				buffR[i] = fr.samples[1][lf++];
			}
			else { //duplicate L -> R if a mono file
				buffR[i] = fr.samples[0][lf++];
			}
			if (lf > lf2) {
				lf = lf1; //loop
			}
			++l;
		}
		else {
			buffL[i] = 0.0f;
			buffR[i] = 0.0f;
		}
	}
}

//==================================================

wavwriter::wavwriter(const std::string& filename, int chan,
	double Fs, double duration) : f(filename)
{
	fw.shouldLogErrorsToConsole(true);
	fw.setBitDepth(24);
	fw.setNumChannels(chan);
	fw.setSampleRate((uint32_t)Fs);
	
	L = (int)round(duration * Fs);
	fw.setNumSamplesPerChannel(L);

	//preclear the data so it writes safely if we fail to fill it
	fw.samples.resize(chan);
	for (int i = 0; i < chan; i++) {
		fw.samples[i].resize(L, 0.0f);
	}

	l = 0;
}
wavwriter::~wavwriter() {}

//Stores a mono block of size len samples to the 
//writer output. If we exceed the duration, discards.
void wavwriter::put(const float* buff, int len)
{
	for (int i = 0; i < len; ++i) {
		if (l < L) {
			fw.samples[0][l++] = buff[i];
		}
	}
}

//Stores a stereo block of size len samples to the 
//writer output. If we exceed the duration, discards.
void wavwriter::put(const float* buffL, const float* buffR, int len)
{
	for (int i = 0; i < len; ++i) {
		if (l < L) {
			if (fw.samples.size() >= 2) {
				fw.samples[1][l] = buffR[i];
			}
			fw.samples[0][l++] = buffL[i];
		}
	}
}

//Writes the stored data to file
void wavwriter::write()
{
	fw.save(f, AudioFileFormat::Wave);
	std::cout << "Wrote to " + f << std::endl;
}

//==================================================

automator::automator(float value) :
	x({ 0 }), v({ value }), l(0), n(0), loops(false)
{
	init();
}

automator::automator(std::vector<long unsigned int> samples,
	std::vector<float> values,
	bool loop) :
	x(samples), v(values), l(0), n(0), loops(loop)
{
	init();
}

automator::~automator() {}

float* automator::connect()
{
	return &p;
}

void automator::init()
{
	p = v[0];
}

void automator::step(int len)
{
	if (loops) {
		bool newval = true;
		while (newval) {
			if (n < x.size()) {
				if (x[n] < l) {
					++n;
					if (n == x.size()) { //loop point
						l -= (x[x.size() - 1] + 1);
						n = 0;
					}
				}
				else {
					newval = false;
				}
			}
			else {
				newval = false;
			}
		}

		if (n == 0) {
			float temp = (float)(((double)(l + 1)) / ((double)(x[0] + 1)));
			p = (v[v.size() - 1] * (1.0f - temp) + v[0] * temp);
		}
		else {
			float temp = (float)(((double)(l - x[n - 1])) / ((double)(x[n] - x[n - 1])));
			p = (v[n - 1] * (1.0f - temp) + v[n] * temp);
		}

		l += len;
	}
	else {
		bool newval = true;
		while (newval) {
			if (n < x.size()) {
				if (x[n] <= l) {
					++n;
				}
				else {
					newval = false;
				}
			}
			else {
				newval = false;
			}
		}

		if (n == 0) {
			p = v[0];
		}
		else if (n == x.size()) {
			p = v[x.size() - 1];
		}
		else {
			float temp = (float)(((double)(l - x[n - 1])) / ((double)(x[n] - x[n - 1])));
			p = (v[n - 1] * (1.0f - temp) + v[n] * temp);
		}

		l += len;
	}
}

iautomator::iautomator(int value) :
	x({ 0 }), v({ value }), l(0), n(0), loops(false)
{
	init();
}

iautomator::iautomator(std::vector<long unsigned int> samples,
	std::vector<int> values,
	bool loop) :
	x(samples), v(values), l(0), n(0), loops(loop)
{
	init();
}

iautomator::~iautomator() {}

int* iautomator::connect()
{
	return &p;
}

void iautomator::init()
{
	p = v[0];
}

void iautomator::step(int len)
{
	if (loops) {
		bool newval = true;
		while (newval) {
			if (n < x.size()) {
				if (x[n] < l) {
					++n;
					if (n == x.size()) { //loop point
						l -= (x[x.size() - 1] + 1);
						n = 0;
					}
				}
				else {
					newval = false;
				}
			}
			else {
				newval = false;
			}
		}

		if (n == 0) {
			p = v[x.size() - 1];
		}
		else {
			p = v[n - 1];
		}
		l += len;
	}
	else {
		bool newval = true;
		while (newval) {
			if (n < x.size()) {
				if (x[n] <= l) {
					++n;
				}
				else {
					newval = false;
				}
			}
			else {
				newval = false;
			}
		}

		if (n == 0) {
			p = v[0];
		}
		else {
			p = v[n - 1];
		}
		l += len;
	}
}


//==================================================

impulse::impulse() : sent(false) {}
impulse::~impulse() {}

//Resets the impulse generator to output an impulse on the next sample
void impulse::reset()
{
	sent = false;
}

//Fills a mono buffer block of len samples from the impulse
//generator
void impulse::get(float* buff, int len)
{
	vset(buff, 0.0f, len);
	if ((len > 0) && !sent) {
		buff[0] = 1.0f;
		sent = true;
	}
}


//==================================================

sinewave::sinewave(float setFrequency, float setGain) :
	g(setGain), dph(setFrequency * constants.pi)
{}
sinewave::~sinewave() {}

//Fills a mono buffer block of len samples from the sine generator
void sinewave::get(float* buff, int len)
{
	for (int i = 0; i < len; ++i) {
		buff[i] = g * sinf(ph);
		ph += dph;
	}
}


//==================================================

testbuff::testbuff(int L) : x(L), l(L, 0)
{
	clear();
}
testbuff::~testbuff() {}

void testbuff::clear()
{
	vset<float>(x, 0.0f, x.size());
	l = 0;
}
void testbuff::put(const float* buff, int len)
{
	for (int i = 0; i < len; ++i) {
		x[l++] = buff[i]; //circularly addressed
	}
}

int testbuff::size() const
{
	return x.size();
}
const float* testbuff::data() const
{
	return x.ptr();
}
float* testbuff::data()
{
	return x.ptr();
}

//==================================================

//Generates a plottable vector of the smoothed power of a signal
//in dB over time
fvec plotSmoothedPower(const float* x, int len, float Nsmo)
{
	fvec p(len); //generate power vector
	vmult(x, x, p.ptr(), len);

	fof smooth; //one-pole power smoother
	smooth.setCoefs(filt::smooth1(Nsmo), true);
	smooth.stepBlock(p, p, len); //apply smoother

	return p;
}

//Generates a plottable vector of the smoothed power of a signal
//in dB over time
fvec plotSmoothedPowerDB(const float* x, int len, float Nsmo, float floordB)
{
	fvec p = plotSmoothedPower(x, len, Nsmo);
	vmax(p.ptr(), powf(10.0f, floordB / 10.0f), p.ptr(), len); //bound
	vpowdB(p.ptr(), p.ptr(), len); //cast

	return p;
}

//==================================================

//Generates a plottable vector of full linear dB magnitude spectra
//computed from the time domain vector x -- creating a temporary FFT
//for testing simplicity.
fvec plotDBSpectrum(const float* x, int len)
{
	int npts = (int)roundf(powf(2.0f, ceil(log2f((float)len))));
	int K = getKForNpts(npts);

	fft F(npts);
	fvec X(K);

	F.rmag(x, X, len);
	vadd(X.ptr(), constants.eps, X.ptr(), K); //stabilize the dB cast
	vmagdB(X, X, K);

	return X;
}

//Generates a plottable vector of log-spaced dB magnitude spectra
//computed from the time domain vector x -- creating a temporary FFT
//and spectral smoother for testing simplicity.
fvec plotLogSmoothedDBSpectrum(const float* x, int len,
	float f1, float f2, float Noct, float OS)
{
	int npts = (int)roundf(powf(2.0f, ceil(log2f((float)len))));
	int K = getKForNpts(npts);
	int Klog = getKlog(K, f1, f2, Noct / OS);

	fft F(npts);
	fvec X(K);
	fvec Y(Klog);
	specsmooth smo(K, Noct);

	F.rmag(x, X, len);
	smo.wtsmooth(X, X, K);
	vmagdB(X, X, K);
	lin2log(X, Y, K, Klog, f1, f2);
	
	return Y;
}

}