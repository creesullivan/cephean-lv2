
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

	wavreader in("C:\\Users\\crees\\Music\\Cephean tests\\doo_riffing.wav", dur, 5.0);
	assert(in.getSampleRate() == Fs);
	wavwriter out("C:\\Users\\crees\\Music\\Cephean tests\\test_output.wav", 1, Fs, dur);

	// Plugin setup -----------------------------
	ensemble plug(Fs);

	iautomator mode(1); //controls
	automator time(50.0f);
	automator pitch(50.0f);
	automator level(50.0f);
	automator rate(50.0f);
	automator dry(50.0f);
	automator wet(50.0f);
	/*
	if (automation) {
		level = automator(		{ 0,		1 * 48000,	2 * 48000 },
								{ 0.0f,		100.0f,		0.0f });
	}
	*/
	fvec inbuff(bsize); //I/O
	fvec outbuff(bsize);

	plug.connect_port(0, (void*)(mode.connect())); //ports
	plug.connect_port(1, (void*)(time.connect())); //ports
	plug.connect_port(2, (void*)(pitch.connect())); //ports
	plug.connect_port(3, (void*)(level.connect())); //ports
	plug.connect_port(4, (void*)(rate.connect())); //ports
	plug.connect_port(5, (void*)(dry.connect())); //ports
	plug.connect_port(6, (void*)(wet.connect())); //ports
	plug.connect_port(7, (void*)inbuff.ptr());
	plug.connect_port(8, (void*)outbuff.ptr());

	plug.init();

	cout << "Loaded plugin, stepping..." << endl;

	// Running
	for (int it = 0; it < iter; ++it) {
		in.get(inbuff, bsize);

		mode.step(bsize);
		time.step(bsize);
		pitch.step(bsize);
		level.step(bsize);
		rate.step(bsize);
		dry.step(bsize);
		wet.step(bsize);

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

	wavreader in("C:\\Users\\crees\\Music\\Cephean tests\\doo_riffing.wav", dur, 5.0);
	assert(in.getSampleRate() == Fs);

	// Plugin setup -----------------------------
	ensemble plug(Fs);

	iautomator mode(1); //controls
	automator time(50.0f);
	automator pitch(50.0f);
	automator level(50.0f);
	automator rate(50.0f);
	automator dry(50.0f);
	automator wet(50.0f);
	/*
	if (automation) {
		level = automator(		{ 0,		1 * 48000,	2 * 48000 },
								{ 0.0f,		100.0f,		0.0f });
	}
	*/
	fvec inbuff(bsize); //I/O
	fvec outbuff(bsize);

	plug.connect_port(0, (void*)(mode.connect())); //ports
	plug.connect_port(1, (void*)(time.connect())); //ports
	plug.connect_port(2, (void*)(pitch.connect())); //ports
	plug.connect_port(3, (void*)(level.connect())); //ports
	plug.connect_port(4, (void*)(rate.connect())); //ports
	plug.connect_port(5, (void*)(dry.connect())); //ports
	plug.connect_port(6, (void*)(wet.connect())); //ports
	plug.connect_port(7, (void*)inbuff.ptr());
	plug.connect_port(8, (void*)outbuff.ptr());

	plug.init();

	cout << "Loaded plugin, stepping..." << endl;

	// Running
	for (int it = 0; it < iter; ++it) {
		in.get(inbuff, bsize);

		mode.step(bsize);
		time.step(bsize);
		pitch.step(bsize);
		level.step(bsize);
		rate.step(bsize);
		dry.step(bsize);
		wet.step(bsize);

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