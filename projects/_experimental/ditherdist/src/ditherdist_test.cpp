
//--------------------------------------------------
// Runs tests on the raw plugin object without an
// LV2 wrapper
//--------------------------------------------------

#include <iostream>
#include "cephean-lv2.h"
#include "cephean-test.h"

#include "ditherdist.h"

using namespace std;
using namespace cephean;


//==================================================

static void runSoundTest(bool automation)
{
	// File I/O setup ---------------------------
	const double dur = 10.0; //seconds
	const float Fs = 48000.0f; //samples/sec
	const int bsize = 128; //samples per block
	const int iter = (int)ceil(dur * Fs / bsize); // to fill the duration

	//wavreader in("C:\\Users\\crees\\Music\\Cephean tests\\doo_riffing.wav", dur, 5.0);
	//wavreader in("C:\\Users\\crees\\Music\\Cephean tests\\guitar_chugs.wav", dur, 5.0);
	//wavreader in("C:\\Users\\crees\\Music\\Cephean tests\\dry_guitar_lead.wav", dur, 48.0);
	wavreader in("C:\\Users\\crees\\Music\\Cephean tests\\long_guitar_chords.wav", dur, 17.0);
	//wavreader in("C:\\Users\\crees\\Music\\Cephean tests\\finger_bass_riff.wav", dur, 25.0);
	assert(in.getSampleRate() == Fs);
	//sinewave in(1289.3f / Fs, 1.0f);

	wavwriter out("C:\\Users\\crees\\Music\\Cephean tests\\test_output.wav", 1, Fs, dur);

	// Plugin setup -----------------------------
	ditherdist plug(Fs);

	plugintester test(PLUG_CONTROLS, PLUG_INPUTS, PLUG_OUTPUTS, bsize);

	//controls
	test.setControlValue(0, 25.0f);		//drive
	test.setControlValue(1, 0.0f);		//bias
	test.setControlValue(2, 0.0f);		//gate
	test.setControlValue(3, 0.0f);		//low
	test.setControlValue(4, 100.0f);	//high
	test.setControlValue(5, 0.0f);		//tone
	test.setControlValue(6, 0.0f);		//dry
	test.setControlValue(7, 50.0f);		//wet
	if (automation) {

	}

	test.connect(&plug);
	plug.init();

	cout << "Loaded plugin, stepping..." << endl;

	// Running
	for (int it = 0; it < iter; ++it) {
		in.get(test.inputData(), bsize);

		test.step(bsize);

		if (it == (iter / 2)) {
			float bpdummy = 0.0f;
		}

		plug.step(bsize);

		out.put(test.outputData(), bsize);
	}

	//Export the audio file
	out.write();

	cout << "done." << endl;
}

//==================================================

static void runStressTest(bool automation)
{
	// File I/O setup ---------------------------
	const double dur = 1000.0; //seconds
	const float Fs = 48000.0f; //samples/sec
	const int bsize = 128; //samples per block
	const int iter = (int)ceil(dur * Fs / bsize); // to fill the duration

	wavreader in("C:\\Users\\crees\\Music\\Cephean tests\\doo_riffing.wav", dur, 5.0);
	assert(in.getSampleRate() == Fs);

	// Plugin setup -----------------------------
	ditherdist plug(Fs);

	plugintester test(PLUG_CONTROLS, PLUG_INPUTS, PLUG_OUTPUTS, bsize);

	//controls
	test.setControlValue(0, 50.0f);		//drive
	test.setControlValue(1, 0.0f);		//bias
	test.setControlValue(2, 0.0f);		//gate
	test.setControlValue(3, 0.0f);		//low
	test.setControlValue(4, 100.0f);	//high
	test.setControlValue(5, 0.0f);		//tone
	test.setControlValue(6, 0.0f);		//dry
	test.setControlValue(7, 50.0f);		//wet
	if (automation) {

	}

	test.connect(&plug);
	plug.init();

	cout << "Loaded plugin, stepping..." << endl;

	// Running
	float maxval = 0.0f;
	for (int it = 0; it < iter; ++it) {
		in.get(test.inputData(), bsize);

		test.step(bsize);
		plug.step(bsize);
	}

	cout << "done." << endl;
}


//==================================================

int main()
{
	//runSoundTest(false);
	runStressTest(false);

	/*
	//TESTING that the timedither object roughly preserves a sine wave, yep!
	float Fs = 48000.0f;
	int npts = 1024;
	sinewave tone(500.0f / Fs);
	fvec x(npts);
	fvec y(npts);
	fvec z(npts);

	timedither dither(npts);
	tone.get(x, npts);
	dither.apply(x, y, npts);
	//NO INTERNAL PROCESSING, testing
	dither.invert(y, z, npts);
	*/

	/*
	//TESTING the transistor waveshaper, yep!
	int N = 1001;
	fvec x(N);
	vincspace(x.ptr(), -1.0f, 2.0f / (N - 1), N);
	fvec y(N);

	transistor ws;
	ws.set(10.0f, 0.3f, 0.1f, true);
	ws.stepBlock(x, y, N);
	*/

	/*
	//TESTING the tone balance filter
	fof filter;
	filter.setCoefs(filt::inverse(filt::balance1(0.5f, 2.0f)), true);
	
	int npts = 4096;
	fvec x(npts);
	fvec y(npts);
	vset(x.ptr(), 0.0f, npts); x[0] = 1.0f; //impulse
	filter.stepBlock(x, y, npts);

	fvec Y = plotDBSpectrum(y, npts);
	*/

	return 0;
}