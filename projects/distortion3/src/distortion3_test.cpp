
//--------------------------------------------------
// Runs tests on the raw plugin object without an
// LV2 wrapper
//--------------------------------------------------

#include <iostream>
#include "cephean-lv2.h"
#include "cephean-test.h"

#include "distortion3.h"

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
	wavreader in("C:\\Users\\crees\\Music\\Cephean tests\\guitar_chugs.wav", dur, 5.0);
	//wavreader in("C:\\Users\\crees\\Music\\Cephean tests\\dry_guitar_lead.wav", dur, 48.0);
	//wavreader in("C:\\Users\\crees\\Music\\Cephean tests\\long_guitar_chords.wav", dur, 17.0);
	//wavreader in("C:\\Users\\crees\\Music\\Cephean tests\\finger_bass_riff.wav", dur, 25.0);
	assert(in.getSampleRate() == Fs);
	wavwriter out("C:\\Users\\crees\\Music\\Cephean tests\\test_output.wav", 1, Fs, dur);

	// Plugin setup -----------------------------
	distortion3 plug(Fs);

	plugintester test(PLUG_CONTROLS, PLUG_INPUTS, PLUG_OUTPUTS, bsize);

	//controls
	test.setControlValue(0, 50.0f);		//drive
	test.setControlValue(1, 50.0f);		//gain
	test.setControlValue(2, 100.0f);	//level
	test.setControlValue(3, 50.0f);		//shape
	test.setControlValue(4, 50.0f);		//tone
	test.setControlValue(5, 50.0f);		//color
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
	//wavreader in("C:\\Users\\crees\\Music\\Cephean tests\\dry_guitar_lead.wav", dur, 48.0);
	wavreader in("C:\\Users\\crees\\Music\\Cephean tests\\long_guitar_chords.wav", dur, 17.0);
	//wavreader in("C:\\Users\\crees\\Music\\Cephean tests\\finger_bass_riff.wav", dur, 25.0);
	assert(in.getSampleRate() == Fs);

	// Plugin setup -----------------------------
	distortion3 plug(Fs);

	plugintester test(PLUG_CONTROLS, PLUG_INPUTS, PLUG_OUTPUTS, bsize);

	//controls
	test.setControlValue(0, 50.0f);		//drive
	test.setControlValue(1, 50.0f);		//gain
	test.setControlValue(2, 100.0f);	//level
	test.setControlValue(3, 50.0f);		//shape
	test.setControlValue(4, 50.0f);		//tone
	test.setControlValue(5, 50.0f);		//color
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

	/* //TESTING DISTORTION CURVE
	const int L = 100;
	const float s = 10.0f;
	fvec x(L);
	fvec z(L);
	fvec y(L);
	fvec T(L);

	vincspace(x.ptr(), 0.0f, 1.0f / (L - 1), L);
	vmult(x.ptr(), s, z.ptr(), L);
	vset(T.ptr(), 1.0f, L);
	vmin(z.ptr(), T.ptr(), z.ptr(), L);

	simsaturator tester(1, L);
	tester.set(s, 0.0f, true);
	tester.clear();
	tester.stepBlock(x, T, y, L);

	float bpdummy = 0.0f;
	*/

	/* //CHECKING DISTORTION FILTER
	const float Fs = 48000.0f;
	const int L = 4800;
	fvec h(L);
	const float f0 = 20000.0f;
	const float alpha = expf(-2.0f * constants.pi * f0 / Fs);
	fvec x(L);
	fvec y(L);
	fvec T(L);

	vset(x.ptr(), 0.0f, L);
	x[0] = 1.0f;
	vset(T.ptr(), 100.0f, L); //really big to basically ignore nonlinearity

	simsaturator tester(1, L);
	tester.set(1.0f, alpha, true);
	tester.clear();
	tester.stepBlock(x, T, y, L);

	fvec Y = plotLogSmoothedDBSpectrum(y, L, 0.0f, 1.0f, 1.0f / 12.0f, 4.0f);
	float bpdummy = 0.0f;
	*/

	return 0;
}