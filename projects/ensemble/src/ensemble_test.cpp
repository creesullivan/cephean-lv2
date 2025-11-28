
//--------------------------------------------------
// Runs tests on the raw plugin object without an
// LV2 wrapper
//--------------------------------------------------

#include <iostream>
#include "cephean-lv2.h"
#include "cephean-test.h"

#include "ensemble.h"

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
	//wavreader in("C:\\Users\\crees\\Music\\Cephean tests\\long_guitar_chords.wav", dur, 17.0);
	wavreader in("C:\\Users\\crees\\Music\\Cephean tests\\finger_bass_riff.wav", dur, 25.0);
	assert(in.getSampleRate() == Fs);
	wavwriter out("C:\\Users\\crees\\Music\\Cephean tests\\test_output.wav", 1, Fs, dur);

	// Plugin setup -----------------------------
	ensemble plug(Fs);

	automator cmode(4.0f); //controls
	automator ctime(50.0f);
	automator cpitch(20.0f);
	automator crate(50.0f);
	automator camount(100.0f);
	automator crolloff(50.0f);
	automator clevel(50.0f);
	if (automation) {
		cmode = automator(		{ 0,		1 * 48000,	2 * 48000 },
								{ 0.0f,		3.0f,		3.0f });
	}
	fvec inbuff(bsize); //I/O
	fvec outbuff(bsize);

	plug.connect_port(0, (void*)(cmode.connect())); //ports
	plug.connect_port(1, (void*)(ctime.connect())); //ports
	plug.connect_port(2, (void*)(cpitch.connect())); //ports
	plug.connect_port(3, (void*)(crate.connect())); //ports
	plug.connect_port(4, (void*)(camount.connect())); //ports
	plug.connect_port(5, (void*)(crolloff.connect())); //ports
	plug.connect_port(6, (void*)(clevel.connect())); //ports
	plug.connect_port(7, (void*)inbuff.ptr());
	plug.connect_port(8, (void*)outbuff.ptr());

	plug.init();

	cout << "Loaded plugin, stepping..." << endl;

	// Running
	for (int it = 0; it < iter; ++it) {
		in.get(inbuff, bsize);

		cmode.step(bsize);
		ctime.step(bsize);
		cpitch.step(bsize);
		crate.step(bsize);
		camount.step(bsize);
		crolloff.step(bsize);
		clevel.step(bsize);

		plug.step(bsize);

		out.put(outbuff, bsize);
	}

	//Export the audio file
	out.write(); //LEFT OFF HERE <-- some sort of weird memory leak in the buff.get() function, hmm

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
	ensemble plug(Fs);

	automator cmode(4.0f); //controls
	automator ctime(50.0f);
	automator cpitch(20.0f);
	automator crate(50.0f);
	automator camount(100.0f);
	automator crolloff(50.0f);
	automator clevel(50.0f);
	/*
	if (automation) {
		level = automator(		{ 0,		1 * 48000,	2 * 48000 },
								{ 0.0f,		100.0f,		0.0f });
	}
	*/
	fvec inbuff(bsize); //I/O
	fvec outbuff(bsize);

	plug.connect_port(0, (void*)(cmode.connect())); //ports
	plug.connect_port(1, (void*)(ctime.connect())); //ports
	plug.connect_port(2, (void*)(cpitch.connect())); //ports
	plug.connect_port(3, (void*)(crate.connect())); //ports
	plug.connect_port(4, (void*)(camount.connect())); //ports
	plug.connect_port(5, (void*)(crolloff.connect())); //ports
	plug.connect_port(6, (void*)(clevel.connect())); //ports
	plug.connect_port(7, (void*)inbuff.ptr());
	plug.connect_port(8, (void*)outbuff.ptr());

	plug.init();

	cout << "Loaded plugin, stepping..." << endl;

	// Running
	for (int it = 0; it < iter; ++it) {
		in.get(inbuff, bsize);

		cmode.step(bsize);
		ctime.step(bsize);
		cpitch.step(bsize);
		crate.step(bsize);
		camount.step(bsize);
		crolloff.step(bsize);
		clevel.step(bsize);

		plug.step(bsize);
	}

	cout << "done." << endl;
}


//==================================================

int main()
{
	//runSoundTest(false);
	runStressTest(false);

	return 0;
}