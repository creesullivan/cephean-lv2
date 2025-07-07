#pragma once

//--------------------------------------------------
// Defines classes and functions for offline testing
// and debugging
//--------------------------------------------------

#include "cephean-lv2.h"
#include "cephean-dsp.h"
#include "cephean-fft.h"

#include "extern/AudioFile/AudioFile.h"
#include <iostream>
#include <vector>

namespace cephean
{

//==================================================

//Easy wrapper around the AudioFile library that constructs
//with a filename, a duration, and an optional time range, and
//manages block-based sample retrievals with automatic looping.
class wavreader
{
public:
	wavreader(const std::string& filename,
		double duration = 10.0, double startTime = 0.0, double endTime = -1.0);
	~wavreader();

	//Returns the sampling rate of the read audio file
	double getSampleRate() const;

	//Fills a mono buffer block of len samples from the file.
	//If we have exceeded the duration, pads with 0s.
	void get(float* buff, int len);

	//Fills a stereo buffer block of len samples from the file.
	//If we have exceeded the duration, pads with 0s, and if the
	//file is mono, duplicates to stereo.
	void get(float* buffL, float* buffR, int len);

private:
	unsigned long int lf1 = 0; //starting sample index in the file
	unsigned long int lf2 = 0; //ending sample index in the file
	unsigned long int Lfile = 0; //length of the file, looping if overrun

	unsigned long int L = 0; //total length of the test duration

	AudioFile<float> fr; //file reader object
	unsigned long int l = 0; //current total sample location
	unsigned long int lf = 0; //current sample location in the file
};

//==================================================

//Easy wrapper around the AudioFile library that constructs
//with a filename, sampling rate, channel count, and duration,
//and manages block-based sample puts before a final write()
//command. 
class wavwriter
{
public:
	wavwriter(const std::string& filename, int setChan=1,
		double setFs = 48000.0, double duration = 10.0);
	~wavwriter();

	//Stores a mono block of size len samples to the 
	//writer output. If we exceed the duration, discards.
	void put(const float* buff, int len);

	//Stores a stereo block of size len samples to the 
	//writer output. If we exceed the duration, discards.
	void put(const float* buffL, const float* buffR, int len);

	//Writes the stored data to file
	void write();

private:
	const std::string f = "";
	int L = 0; //total file length in samples

	AudioFile<float> fw; //file writer object
	int l = 0; //current write location in the file
};


//==================================================

//Class that grabs hold of a parameter control and automates
//it with a piecewise linear drag throughout the recording test.
class automator
{
public:
	automator(float value);
	automator(std::vector<long unsigned int> samples,
		std::vector<float> values,
		bool loop = true);
	~automator();

	float* connect();

	void init();
	void step(int len);

private:
	float p = 0.0f; //automated parameter

	long unsigned int l = 0;
	size_t n = 0;
	std::vector<long unsigned int> x;
	std::vector<float> v;
	bool loops = true;
};

//Class that grabs hold of an int parameter control and automates
//it with a stairstep selection throughout the recording test.
class iautomator
{
public:
	iautomator(int value);
	iautomator(std::vector<long unsigned int> samples,
		std::vector<int> values,
		bool loop = true);
	~iautomator();

	int* connect();

	void init();
	void step(int len);

private:
	int p = 0; //automated parameter

	long unsigned int l = 0;
	size_t n = 0;
	std::vector<long unsigned int> x;
	std::vector<int> v;
	bool loops = true;
};


//==================================================

//Class that generates a single unity gain impulse
class impulse
{
public:
	impulse();
	~impulse();

	//Resets the impulse generator to output an impulse on the next sample
	void reset();

	//Fills a mono buffer block of len samples from the impulse
	//generator
	void get(float* buff, int len);
private:
	bool sent = false;
};


//==================================================

//Class that generates a pure tone with fixed frequency and gain
class sinewave
{
public:
	sinewave(float setFrequency, float setGain=1.0f);
	~sinewave();

	//Fills a mono buffer block of len samples from the sine generator
	void get(float* buff, int len);
private:
	const float dph = 0.01f; //phase step in radians
	const float g = 1.0f; //gain multiplier

	fastcircfloat ph{ 2.0f * constants.pi, 0.0f }; //phase
};


//==================================================

//Class that automatically stores mono blocks of data to a circular
//buffer for manipulation later.
class testbuff
{
public:
	testbuff(int L);
	~testbuff();

	void clear();
	void put(const float* buff, int len);

	int size() const;
	const float* data() const;
	float* data();

private:
	fvec x; //data vector
	fastcircint l;
};


//==================================================

//Generates a plottable vector of the smoothed power of a signal
//over time
fvec plotSmoothedPower(const float* x, int len, float Nsmo);

//Generates a plottable vector of the smoothed power of a signal
//in dB over time
fvec plotSmoothedPowerDB(const float* x, int len, float Nsmo, float floordB=-100.0f);

//==================================================

//Generates a plottable vector of full linear dB magnitude spectra
//computed from the time domain vector x -- creating a temporary FFT
//for testing simplicity.
fvec plotDBSpectrum(const float* x, int len);

//Generates a plottable vector of log-spaced dB magnitude spectra
//computed from the time domain vector x -- creating a temporary FFT
//and spectral smoother for testing simplicity.
fvec plotLogSmoothedDBSpectrum(const float* x, int len,
	float f1 = 0.0f, float f2 = 1.0f, float Noct = 1 / 3.0f, float OS = 16.0f);

}