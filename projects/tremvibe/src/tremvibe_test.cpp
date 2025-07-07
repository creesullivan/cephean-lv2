
//--------------------------------------------------
// Runs tests on the raw plugin object without an
// LV2 wrapper
//--------------------------------------------------

#include <iostream>
#include "cephean-lv2.h"
#include "cephean-test.h"

#include "tremvibe.h"

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
	tremvibe plug(Fs);

	automator tremelo(25.0f); //controls
	automator vibrato(50.0f);
	automator rate(100.0f);
	automator random(100.0f);
	if (automation) {
		tremelo = automator(	{ 0,		1 * 48000,	2 * 48000 },
								{ 0.0f,		100.0f,		0.0f });
		vibrato = automator(	{ 0,		1 * 48000,	2 * 48000 },
								{ 0.0f,		100.0f,		0.0f });
		rate = automator(		{ 0,		1 * 48000,	2 * 48000 },
								{ 0.0f,		100.0f,		0.0f });
		random = automator(		{ 0,		1 * 48000,	2 * 48000 },
								{ 0.0f,		100.0f,		0.0f });
	}
	fvec inbuff(bsize); //I/O
	fvec outbuff(bsize);

	plug.connect_port(0, (void*)(tremelo.connect())); //ports
	plug.connect_port(1, (void*)(vibrato.connect()));
	plug.connect_port(2, (void*)(rate.connect()));
	plug.connect_port(3, (void*)(random.connect()));
	plug.connect_port(4, (void*)inbuff.ptr());
	plug.connect_port(5, (void*)outbuff.ptr());

	plug.init();

	cout << "Loaded plugin, stepping..." << endl;

	// Running
	for (int it = 0; it < iter; ++it) {
		in.get(inbuff, bsize);

		tremelo.step(bsize);
		vibrato.step(bsize);
		rate.step(bsize);
		random.step(bsize);

		plug.step(bsize);

		out.put(outbuff, bsize);
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

	wavreader in("C:\\Users\\crees\\Music\\Cephean tests\\doo_riffing.wav", dur, 5.0, 15.0);
	assert(in.getSampleRate() == Fs);

	// Plugin setup -----------------------------
	tremvibe plug(Fs);

	automator tremelo(45.0f); //controls
	automator vibrato(10.0f);
	automator rate(40.0f);
	automator random(50.0f);
	if (automation) {
		tremelo = automator({ 0,		1 * 48000,	2 * 48000 },
			{ 0.0f,		100.0f,		0.0f });
		vibrato = automator({ 0,		1 * 48000,	2 * 48000 },
			{ 0.0f,		100.0f,		0.0f });
		rate = automator({ 0,		1 * 48000,	2 * 48000 },
			{ 0.0f,		100.0f,		0.0f });
		random = automator({ 0,		1 * 48000,	2 * 48000 },
			{ 0.0f,		100.0f,		0.0f });
	}
	fvec inbuff(bsize); //I/O
	fvec outbuff(bsize);

	plug.connect_port(0, (void*)(tremelo.connect())); //ports
	plug.connect_port(1, (void*)(vibrato.connect()));
	plug.connect_port(2, (void*)(rate.connect()));
	plug.connect_port(3, (void*)(random.connect()));
	plug.connect_port(4, (void*)inbuff.ptr());
	plug.connect_port(5, (void*)outbuff.ptr());

	plug.init();

	cout << "Loaded plugin, stepping..." << endl;

	// Running
	for (int it = 0; it < iter; ++it) {
		in.get(inbuff, bsize);

		tremelo.step(bsize);
		vibrato.step(bsize);
		rate.step(bsize);
		random.step(bsize);

		plug.step(bsize);
	}

	cout << "done." << endl;
}


//==================================================

int main()
{
	runSoundTest(false);
	//runStressTest(true);

	return 0;
}