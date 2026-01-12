
//--------------------------------------------------
// Runs tests on the raw plugin object without an
// LV2 wrapper
//--------------------------------------------------

#include <iostream>
#include "cephean-lv2.h"
#include "cephean-test.h"

#include "fourwayswitch.h"

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

	wavreader in("C:\\Users\\crees\\Music\\Cephean tests\\doo_riffing.wav", dur, 5.0);
	//wavreader in("C:\\Users\\crees\\Music\\Cephean tests\\guitar_chugs.wav", dur, 5.0);
	//wavreader in("C:\\Users\\crees\\Music\\Cephean tests\\dry_guitar_lead.wav", dur, 48.0);
	//wavreader in("C:\\Users\\crees\\Music\\Cephean tests\\long_guitar_chords.wav", dur, 17.0);
	//wavreader in("C:\\Users\\crees\\Music\\Cephean tests\\finger_bass_riff.wav", dur, 25.0);
	assert(in.getSampleRate() == Fs);
	wavwriter out("C:\\Users\\crees\\Music\\Cephean tests\\test_output.wav", 1, Fs, dur);

	// Plugin setup -----------------------------
	fourwayswitch plug(Fs);

	plugintester test(PLUG_CONTROLS, PLUG_INPUTS, PLUG_OUTPUTS, bsize);

	//controls
	test.setControlValue(0, 2.0f);		//position
	if (automation) {

	}

	test.connect(&plug);
	plug.init();

	cout << "Loaded plugin, stepping..." << endl;

	// Running
	for (int it = 0; it < iter; ++it) {
		vset(test.inputData(0), 0.0f, bsize);
		vset(test.inputData(1), 0.0f, bsize);
		vset(test.inputData(2), 0.0f, bsize);
		vset(test.inputData(3), 0.0f, bsize);
		in.get(test.inputData(1), bsize);

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
	//wavreader in("C:\\Users\\crees\\Music\\Cephean tests\\guitar_chugs.wav", dur, 5.0);
	//wavreader in("C:\\Users\\crees\\Music\\Cephean tests\\dry_guitar_lead.wav", dur, 48.0);
	//wavreader in("C:\\Users\\crees\\Music\\Cephean tests\\long_guitar_chords.wav", dur, 17.0);
	//wavreader in("C:\\Users\\crees\\Music\\Cephean tests\\finger_bass_riff.wav", dur, 25.0);
	assert(in.getSampleRate() == Fs);

	// Plugin setup -----------------------------
	fourwayswitch plug(Fs);

	plugintester test(PLUG_CONTROLS, PLUG_INPUTS, PLUG_OUTPUTS, bsize);

	//controls
	test.setControlValue(0, 1.0f);		//position
	if (automation) {

	}

	test.connect(&plug);
	plug.init();

	cout << "Loaded plugin, stepping..." << endl;

	// Running
	float maxval = 0.0f;
	for (int it = 0; it < iter; ++it) {
		vset(test.inputData(0), 0.0f, bsize);
		vset(test.inputData(1), 0.0f, bsize);
		vset(test.inputData(2), 0.0f, bsize);
		vset(test.inputData(3), 0.0f, bsize);
		in.get(test.inputData(1), bsize);

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


	return 0;
}