
//--------------------------------------------------
// Runs tests on the raw plugin object without an
// LV2 wrapper
//--------------------------------------------------

#include <iostream>
#include "cephean-lv2.h"
#include "cephean-test.h"

#include "wahwah.h"

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
	wavwriter out("C:\\Users\\crees\\Music\\Cephean tests\\test_output.wav", 1, Fs, dur);

	// Plugin setup -----------------------------
	wahwah plug(Fs);

	plugintester test(PLUG_CONTROLS, PLUG_INPUTS, PLUG_OUTPUTS, bsize);

	//controls
	test.setControlValue(0, 50.0f);				//control
	test.setControlValue(1, 50.0f);				//freq1
	test.setControlValue(2, 50.0f);				//freq2
	test.setControlValue(3, 25.0f);				//res1
	test.setControlValue(4, 75.0f);				//res2
	test.setControlValue(5, 25.0f);				//rolloff
	test.setControlValue(6, 50.0f);				//level
	if (automation) {
		test.setControlAutomator(0, automator({ 0, 24000, 48000, 148000, 172000, 196000 },
											  { 0.0f, 100.0f, 0.0f, 0.0f, 100.0f, 0.0f }, false)); //simulates a wah pedal
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
	//wavreader in("C:\\Users\\crees\\Music\\Cephean tests\\dry_guitar_lead.wav", dur, 48.0);
	wavreader in("C:\\Users\\crees\\Music\\Cephean tests\\long_guitar_chords.wav", dur, 17.0);
	//wavreader in("C:\\Users\\crees\\Music\\Cephean tests\\finger_bass_riff.wav", dur, 25.0);
	assert(in.getSampleRate() == Fs);

	// Plugin setup -----------------------------
	wahwah plug(Fs);

	plugintester test(PLUG_CONTROLS, PLUG_INPUTS, PLUG_OUTPUTS, bsize);

	//controls
	test.setControlValue(0, 50.0f);				//control
	test.setControlValue(1, 50.0f);				//freq1
	test.setControlValue(2, 50.0f);				//freq2
	test.setControlValue(3, 50.0f);				//res1
	test.setControlValue(4, 50.0f);				//res2
	test.setControlValue(5, 0.0f);				//rolloff
	test.setControlValue(6, 50.0f);				//level
	if (automation) {
		test.setControlAutomator(0, automator({ 0, 48000 }, { 0.0f, 100.0f })); //simulates a wah pedal
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
	runSoundTest(true);
	//runStressTest(true);

	/*
	//TESTING a way to smoothly embed gain and resonance into the wah design

	const int N = 2048;
	const int K = getKForNpts(N);
	const float f = 0.5f; //to avoid log plotting

	fvec h(N);

	float r = expf(constants.db2nats * 0.0f);
	r = sqrtf(2.0f * r);
	impz(filt::resonator2(f, r), h, N);
	vmult(h.ptr(), r / 2.0f, h.ptr(), N);
	fvec H0 = plotDBSpectrum(h, N);
	vbound(H0.ptr(), -20.0f, 20.0f, H0.ptr(), K);

	r = expf(constants.db2nats * 3.0f);
	r = sqrtf(2.0f * r);
	impz(filt::resonator2(f, r), h, N);
	vmult(h.ptr(), r / 2.0f, h.ptr(), N);
	fvec H3 = plotDBSpectrum(h, N);
	vbound(H3.ptr(), -20.0f, 20.0f, H3.ptr(), K);

	r = expf(constants.db2nats * 6.0f);
	r = sqrtf(2.0f * r);
	impz(filt::resonator2(f, r), h, N);
	vmult(h.ptr(), r / 2.0f, h.ptr(), N);
	fvec H6 = plotDBSpectrum(h, N);
	vbound(H6.ptr(), -20.0f, 20.0f, H6.ptr(), K);

	r = expf(constants.db2nats * 9.0f);
	r = sqrtf(2.0f * r);
	impz(filt::resonator2(f, r), h, N);
	vmult(h.ptr(), r / 2.0f, h.ptr(), N);
	fvec H9 = plotDBSpectrum(h, N);
	vbound(H9.ptr(), -20.0f, 20.0f, H9.ptr(), K);

	r = expf(constants.db2nats * 12.0f);
	r = sqrtf(2.0f * r);
	impz(filt::resonator2(f, r), h, N);
	vmult(h.ptr(), r / 2.0f, h.ptr(), N);
	fvec H12 = plotDBSpectrum(h, N);
	vbound(H12.ptr(), -20.0f, 20.0f, H12.ptr(), K);

	r = expf(constants.db2nats * 15.0f);
	r = sqrtf(2.0f * r);
	impz(filt::resonator2(f, r), h, N);
	vmult(h.ptr(), r / 2.0f, h.ptr(), N);
	fvec H15 = plotDBSpectrum(h, N);
	vbound(H15.ptr(), -20.0f, 20.0f, H15.ptr(), K);

	r = expf(constants.db2nats * 18.0f);
	r = sqrtf(2.0f * r);
	impz(filt::resonator2(f, r), h, N);
	vmult(h.ptr(), r / 2.0f, h.ptr(), N);
	fvec H18 = plotDBSpectrum(h, N);
	vbound(H18.ptr(), -20.0f, 20.0f, H18.ptr(), K);
	*/

	return 0;
}