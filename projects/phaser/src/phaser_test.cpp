
//--------------------------------------------------
// Runs tests on the raw plugin object without an
// LV2 wrapper
//--------------------------------------------------

#include <iostream>
#include "cephean-lv2.h"
#include "cephean-test.h"

#include "phaser.h"

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
	phaser plug(Fs);

	automator depth(0.0f); //controls
	automator stop1(25.0f);
	automator stop2(75.0f);
	automator rate(50.0f);
	automator tone(100.0f);
	automator level(50.0f);
	if (automation) {
		depth = automator({ 0,		1 * 48000,	2 * 48000 },
			{ 0.0f,		100.0f,		0.0f });
		stop1 = automator({ 0,		1 * 48000,	2 * 48000 },
			{ 0.0f,		0.0f,		0.0f });
		stop2 = automator({ 0,		1 * 48000,	2 * 48000 },
			{ 100.0f,		100.0f,		100.0f });
		rate = automator({ 0,		1 * 48000,	2 * 48000 },
			{ 50.0f,		50.0f,		50.0f });
		tone = automator({ 0,		1 * 48000,	2 * 48000 },
			{ 50.0f,		50.0f,		50.0f });
		level = automator({0,		1 * 48000,	2 * 48000},
			{ 50.0f,		50.0f,		50.0f });
	}
	fvec inbuff(bsize); //I/O
	fvec outbuff(bsize);

	plug.connect_port(0, (void*)(depth.connect())); //ports
	plug.connect_port(1, (void*)(stop1.connect()));
	plug.connect_port(2, (void*)(stop2.connect()));
	plug.connect_port(3, (void*)(rate.connect()));
	plug.connect_port(4, (void*)(tone.connect()));
	plug.connect_port(5, (void*)(level.connect()));
	plug.connect_port(6, (void*)inbuff.ptr());
	plug.connect_port(7, (void*)outbuff.ptr());

	plug.init();

	cout << "Loaded plugin, stepping..." << endl;

	// Running
	for (int it = 0; it < iter; ++it) {
		in.get(inbuff, bsize);

		depth.step(bsize);
		stop1.step(bsize);
		stop2.step(bsize);
		rate.step(bsize);
		tone.step(bsize);
		level.step(bsize);

		if (it == 1000)
		{
			float bpdummy = 0.0f;
		}

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

	wavreader in("C:\\Users\\crees\\Music\\Cephean tests\\doo_riffing.wav", dur, 5.0, 15.0f);
	assert(in.getSampleRate() == Fs);

	// Plugin setup -----------------------------
	phaser plug(Fs);

	automator depth(100.0f); //controls
	automator stop1(0.0f);
	automator stop2(100.0f);
	automator rate(20.0f);
	automator tone(100.0f);
	automator level(50.0f);
	if (automation) {
		depth = automator({ 0,		1 * 48000,	2 * 48000 },
			{ 0.0f,		100.0f,		0.0f });
		stop1 = automator({ 0,		1 * 48000,	2 * 48000 },
			{ 0.0f,		100.0f,		0.0f });
		stop2 = automator({ 0,		1 * 48000,	2 * 48000 },
			{ 0.0f,		100.0f,		0.0f });
		rate = automator({ 0,		1 * 48000,	2 * 48000 },
			{ 0.0f,		100.0f,		0.0f });
		tone = automator({ 0,		1 * 48000,	2 * 48000 },
			{ 0.0f,		100.0f,		0.0f });
		level = automator({ 0,		1 * 48000,	2 * 48000 },
			{ 0.0f,		100.0f,		0.0f });
	}
	fvec inbuff(bsize); //I/O
	fvec outbuff(bsize);

	plug.connect_port(0, (void*)(depth.connect())); //ports
	plug.connect_port(1, (void*)(stop1.connect()));
	plug.connect_port(2, (void*)(stop2.connect()));
	plug.connect_port(3, (void*)(rate.connect()));
	plug.connect_port(4, (void*)(tone.connect()));
	plug.connect_port(5, (void*)(level.connect()));
	plug.connect_port(6, (void*)inbuff.ptr());
	plug.connect_port(7, (void*)outbuff.ptr());

	plug.init();

	cout << "Loaded plugin, stepping..." << endl;

	// Running
	for (int it = 0; it < iter; ++it) {
		in.get(inbuff, bsize);

		depth.step(bsize);
		stop1.step(bsize);
		stop2.step(bsize);
		rate.step(bsize);
		tone.step(bsize);
		level.step(bsize);

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