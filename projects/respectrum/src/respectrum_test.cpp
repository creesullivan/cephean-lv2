
//--------------------------------------------------
// Runs tests on the raw plugin object without an
// LV2 wrapper
//--------------------------------------------------

#include <iostream>
#include "cephean-lv2.h"
#include "cephean-test.h"

#include "respectrum.h"

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
	respectrum plug(Fs);

	automator gateLow(0.0f); //controls
	automator gateHigh(0.0f);
	automator timbre(50.0f);
	automator attack(50.0f);
	automator time(50.0f);
	automator level(50.0f);
	if (automation) {
		gateLow = automator(	{ 0,		1 * 48000,	2 * 48000 },
								{ 0.0f,		0.0f,		0.0f });
		gateHigh = automator(	{ 0,		1 * 48000,	2 * 48000 },
								{ 0.0f,		0.0f,		0.0f });
		timbre = automator(		{ 0,		1 * 48000,	2 * 48000 },
								{ 50.0f,		50.0f,		50.0f });
		attack = automator(		{ 0,		5 * 48000,	10 * 48000 },
								{ 50.0f,		100.0f,		50.0f });
		time = automator(		{ 0,		1 * 48000,	2 * 48000 },
								{ 50.0f,		50.0f,		50.0f });
		level = automator(		{ 0,		1 * 48000,	2 * 48000 },
								{ 50.0f,		50.0f,		50.0f });
	}
	fvec inbuff(bsize); //I/O
	fvec outbuff(bsize);

	plug.connect_port(0, (void*)(gateLow.connect())); //ports
	plug.connect_port(1, (void*)(gateHigh.connect())); //ports
	plug.connect_port(2, (void*)(timbre.connect())); //ports
	plug.connect_port(3, (void*)(attack.connect())); //ports
	plug.connect_port(4, (void*)(time.connect())); //ports
	plug.connect_port(5, (void*)(level.connect())); //ports
	plug.connect_port(6, (void*)inbuff.ptr());
	plug.connect_port(7, (void*)outbuff.ptr());

	plug.init();

	cout << "Loaded plugin, stepping..." << endl;

	// Running
	for (int it = 0; it < iter; ++it) {
		in.get(inbuff, bsize);

		gateLow.step(bsize);
		gateHigh.step(bsize);
		timbre.step(bsize);
		attack.step(bsize);
		time.step(bsize);
		level.step(bsize);

		if (it == (iter / 6)) {
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
	const double dur = 100.0; //seconds
	const float Fs = 48000.0f; //samples/sec
	const int bsize = 128; //samples per block
	const int iter = (int)ceil(dur * Fs / bsize); // to fill the duration

	wavreader in("C:\\Users\\crees\\Music\\Cephean tests\\doo_riffing.wav", dur, 5.0, 15.0);
	assert(in.getSampleRate() == Fs);

	// Plugin setup -----------------------------
	respectrum plug(Fs);

	automator gateLow(30.0f); //controls
	automator gateHigh(10.0f);
	automator timbre(90.0f);
	automator attack(50.0f);
	automator time(50.0f);
	automator level(50.0f);
	if (automation) {
		gateLow = automator({ 0,		1 * 48000,	2 * 48000 },
			{ 0.0f,		100.0f,		0.0f });
		gateHigh = automator({ 0,		1 * 48000,	2 * 48000 },
			{ 0.0f,		100.0f,		0.0f });
		timbre = automator({ 0,		1 * 48000,	2 * 48000 },
			{ 0.0f,		100.0f,		0.0f });
		attack = automator({ 0,		1 * 48000,	2 * 48000 },
			{ 0.0f,		100.0f,		0.0f });
		time = automator({ 0,		1 * 48000,	2 * 48000 },
			{ 0.0f,		100.0f,		0.0f });
		level = automator({ 0,		1 * 48000,	2 * 48000 },
			{ 0.0f,		100.0f,		0.0f });
	}
	fvec inbuff(bsize); //I/O
	fvec outbuff(bsize);

	plug.connect_port(0, (void*)(gateLow.connect())); //ports
	plug.connect_port(1, (void*)(gateHigh.connect())); //ports
	plug.connect_port(2, (void*)(timbre.connect())); //ports
	plug.connect_port(3, (void*)(attack.connect())); //ports
	plug.connect_port(4, (void*)(time.connect())); //ports
	plug.connect_port(5, (void*)(level.connect())); //ports
	plug.connect_port(6, (void*)inbuff.ptr());
	plug.connect_port(7, (void*)outbuff.ptr());

	plug.init();

	cout << "Loaded plugin, stepping..." << endl;

	// Running
	for (int it = 0; it < iter; ++it) {
		in.get(inbuff, bsize);

		gateLow.step(bsize);
		gateHigh.step(bsize);
		timbre.step(bsize);
		attack.step(bsize);
		time.step(bsize);
		level.step(bsize);

		plug.step(bsize);
	}

	cout << "done." << endl;
}




//==================================================

int main()
{
	runSoundTest(true);
	//runStressTest(true);

	return 0;
}