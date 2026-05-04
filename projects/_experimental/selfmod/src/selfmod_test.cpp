
//--------------------------------------------------
// Runs tests on the raw plugin object without an
// LV2 wrapper
//--------------------------------------------------

#include <iostream>
#include "cephean-lv2.h"
#include "cephean-test.h"

#include "selfmod.h"

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
	wavreader in("C:\\Users\\crees\\Music\\Cephean tests\\dry_guitar_lead.wav", dur, 48.0);
	//wavreader in("C:\\Users\\crees\\Music\\Cephean tests\\long_guitar_chords.wav", dur, 17.0);
	//wavreader in("C:\\Users\\crees\\Music\\Cephean tests\\finger_bass_riff.wav", dur, 25.0);
	assert(in.getSampleRate() == Fs);
	wavwriter out("C:\\Users\\crees\\Music\\Cephean tests\\test_output.wav", 1, Fs, dur);

	// Plugin setup -----------------------------
	selfmod plug(Fs);

	plugintester test(PLUG_CONTROLS, PLUG_INPUTS, PLUG_OUTPUTS, bsize);

	//controls
	test.setControlValue(0, 10.0f);		//depth
	test.setControlValue(1, 50.0f);		//bandwidth
	test.setControlValue(2, 0.0f);		//rolloff
	test.setControlValue(3, 50.0f);		//wet
	test.setControlValue(4, 0.0f);		//dry
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

	//wavreader in("C:\\Users\\crees\\Music\\Cephean tests\\doo_riffing.wav", dur, 5.0);
	//wavreader in("C:\\Users\\crees\\Music\\Cephean tests\\guitar_chugs.wav", dur, 5.0);
	wavreader in("C:\\Users\\crees\\Music\\Cephean tests\\dry_guitar_lead.wav", dur, 48.0);
	//wavreader in("C:\\Users\\crees\\Music\\Cephean tests\\long_guitar_chords.wav", dur, 17.0);
	//wavreader in("C:\\Users\\crees\\Music\\Cephean tests\\finger_bass_riff.wav", dur, 25.0);
	assert(in.getSampleRate() == Fs);

	// Plugin setup -----------------------------
	selfmod plug(Fs);

	plugintester test(PLUG_CONTROLS, PLUG_INPUTS, PLUG_OUTPUTS, bsize);

	//controls
	test.setControlValue(0, 50.0f);		//level
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
	runSoundTest(false);
	//runStressTest(false);

	/*
	//TESTING -- check the Hilbert transform coefficients
	sofcasc<2>::coefs coef0;
	sofcasc<4>::coefs coef90;
	filt::hilbert7(coef0, coef90);

	int K = 1000;
	fvec f(K);
	vlinspace(f.ptr(), logf(50.0f / 24000.0f), logf(20000.0f / 24000.0f), K);
	vexp(f, f, K);

	cfvec H0(K);
	cfvec H90(K);

	freqz<2>(coef0, f, H0, K);
	freqz<4>(coef90, f, H90, K);

	fvec dist0(K);
	fvec dist90(K);
	fvec ang0(K);
	fvec ang90(K);
	fvec dang(K);

	vabs(H0.ptr(), dist0.ptr(), K);
	vmagdB(dist0, dist0, K);
	vbound(dist0.ptr(), -20.0f, 20.0f, dist0.ptr(), K);

	vabs(H90.ptr(), dist90.ptr(), K);
	vmagdB(dist90, dist90, K);
	vbound(dist90.ptr(), -20.0f, 20.0f, dist90.ptr(), K);

	varg(H0, dang, K);
	vunwrap(dang, ang0, K);
	varg(H90, dang, K);
	vunwrap(dang, ang90, K);
	vsub(ang90.ptr(), ang0.ptr(), dang.ptr(), K);
	vmult(dang.ptr(), 180.0f / constants.pi, dang.ptr(), K); //to degrees
	*/
	return 0;
}