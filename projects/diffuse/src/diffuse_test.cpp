
//--------------------------------------------------
// Runs tests on the raw plugin object without an
// LV2 wrapper
//--------------------------------------------------

#include <iostream>
#include "cephean-lv2.h"
#include "cephean-test.h"

#include "diffuse.h"

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
	diffuse plug(Fs);

	//LEFT OFF HERE <--- sounds like it all works! fine tune the control ranges
	//and verify automation performance, then export

	iautomator seed(0); //controls
	automator spread(50.0f);
	automator color(50.0f);
	automator wet(50.0f);
	automator dry(0.0f);
	if (automation) {
		//seed = iautomator({ 0,		1 * 48000,	2 * 48000 },
		//	{ 0, 1, 2 });
		//spread = automator(		{ 0,		1 * 48000,	2 * 48000 },
		//						{ 0.0f,		100.0f,		0.0f });
		//color = automator({ 0,		1 * 48000,	2 * 48000 },
		//	{ 0.0f,		100.0f,		0.0f });
		//wet = automator({ 0,		1 * 48000,	2 * 48000 },
		//	{ 0.0f,		100.0f,		0.0f });
		//dry = automator({ 0,		1 * 48000,	2 * 48000 },
		//	{ 100.0f,		0.0f,		100.0f });
	}
	fvec inbuff(bsize); //I/O
	fvec outbuff(bsize);

	plug.connect_port(0, (void*)(seed.connect())); //ports
	plug.connect_port(1, (void*)(spread.connect()));
	plug.connect_port(2, (void*)(color.connect()));
	plug.connect_port(3, (void*)(wet.connect()));
	plug.connect_port(4, (void*)(dry.connect()));
	plug.connect_port(5, (void*)inbuff.ptr());
	plug.connect_port(6, (void*)outbuff.ptr());

	plug.init();

	cout << "Loaded plugin, stepping..." << endl;

	// Running
	for (int it = 0; it < iter; ++it) {
		in.get(inbuff, bsize);

		seed.step(bsize);
		spread.step(bsize);
		color.step(bsize);
		wet.step(bsize);
		dry.step(bsize);

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
	diffuse plug(Fs);

	iautomator seed(0); //controls
	automator spread(50.0f);
	automator color(50.0f);
	automator wet(50.0f);
	automator dry(0.0f);
	if (automation) {
		//seed = iautomator({ 0,		1 * 48000,	2 * 48000 },
		//	{ 0, 1, 2 });
		spread = automator(		{ 0,		1 * 48000,	2 * 48000 },
								{ 0.0f,		100.0f,		0.0f });
		color = automator({ 0,		1 * 48000,	2 * 48000 },
			{ 0.0f,		100.0f,		0.0f });
		wet = automator({ 0,		1 * 48000,	2 * 48000 },
			{ 0.0f,		100.0f,		0.0f });
		dry = automator({ 0,		1 * 48000,	2 * 48000 },
			{ 100.0f,		0.0f,		100.0f });
	}
	fvec inbuff(bsize); //I/O
	fvec outbuff(bsize);

	plug.connect_port(0, (void*)(seed.connect())); //ports
	plug.connect_port(1, (void*)(spread.connect()));
	plug.connect_port(2, (void*)(color.connect()));
	plug.connect_port(3, (void*)(wet.connect()));
	plug.connect_port(4, (void*)(dry.connect()));
	plug.connect_port(5, (void*)inbuff.ptr());
	plug.connect_port(6, (void*)outbuff.ptr());

	plug.init();

	cout << "Loaded plugin, stepping..." << endl;

	// Running
	for (int it = 0; it < iter; ++it) {
		in.get(inbuff, bsize);

		seed.step(bsize);
		spread.step(bsize);
		color.step(bsize);
		wet.step(bsize);
		dry.step(bsize);

		plug.step(bsize);
	}

	cout << "done." << endl;
}


//==================================================

int main()
{
	//runSoundTest(false);
	runStressTest(true);

	return 0;
}