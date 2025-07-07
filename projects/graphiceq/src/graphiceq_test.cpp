
//--------------------------------------------------
// Runs tests on the raw plugin object without an
// LV2 wrapper
//--------------------------------------------------

#include <iostream>
#include "cephean-lv2.h"
#include "cephean-test.h"

#include "graphiceq.h"

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
	graphiceq plug(Fs);

	automator level(50.0f); //controls
	if (automation) {
		level = automator(		{ 0,		1 * 48000,	2 * 48000 },
								{ 0.0f,		100.0f,		0.0f });
	}
	fvec inbuff(bsize); //I/O
	fvec outbuff(bsize);

	plug.connect_port(0, (void*)(level.connect())); //ports
	plug.connect_port(1, (void*)inbuff.ptr());
	plug.connect_port(2, (void*)outbuff.ptr());

	plug.init();

	cout << "Loaded plugin, stepping..." << endl;

	// Running
	for (int it = 0; it < iter; ++it) {
		in.get(inbuff, bsize);

		level.step(bsize);

		plug.step(bsize);

		out.put(outbuff, bsize);
	}

	//Export the audio file
	out.write();

	cout << "done." << endl;
}


//==================================================

int main()
{
	runSoundTest(true);

	return 0;
}