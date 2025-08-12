
//--------------------------------------------------
// Runs tests on the raw plugin object without an
// LV2 wrapper
//--------------------------------------------------

#include <iostream>
#include "cephean-lv2.h"
#include "cephean-test.h"

#include "octaver.h"

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
	octaver plug(Fs);

	automator downone(50.0f); //controls
	automator dry(50.0f);
	automator uphalf(50.0f);
	automator upone(50.0f);
	automator uptwo(50.0f);
	automator color(50.0f);
	if (automation) {
		downone = automator({ 0,		1 * 48000,	2 * 48000 },
				{ 0.0f,		100.0f,		0.0f });
		dry = automator({ 0,		1 * 48000,	2 * 48000 },
			{ 0.0f,		100.0f,		0.0f });
		uphalf = automator({ 0,		1 * 48000,	2 * 48000 },
			{ 0.0f,		100.0f,		0.0f });
		upone = automator({ 0,		1 * 48000,	2 * 48000 },
			{ 0.0f,		100.0f,		0.0f });
		uptwo = automator({ 0,		1 * 48000,	2 * 48000 },
			{ 0.0f,		100.0f,		0.0f });
		color = automator({ 0,		1 * 48000,	2 * 48000 },
			{ 0.0f,		100.0f,		0.0f });
	}
	fvec inbuff(bsize); //I/O
	fvec outbuff(bsize);

	plug.connect_port(0, (void*)(downone.connect())); //ports
	plug.connect_port(1, (void*)(dry.connect()));
	plug.connect_port(2, (void*)(uphalf.connect()));
	plug.connect_port(3, (void*)(upone.connect()));
	plug.connect_port(4, (void*)(uptwo.connect()));
	plug.connect_port(5, (void*)(color.connect()));
	plug.connect_port(6, (void*)inbuff.ptr());
	plug.connect_port(7, (void*)outbuff.ptr());

	plug.init();

	cout << "Loaded plugin, stepping..." << endl;

	// Running
	for (int it = 0; it < iter; ++it) {
		in.get(inbuff, bsize);

		downone.step(bsize);
		dry.step(bsize);
		uphalf.step(bsize);
		upone.step(bsize);
		uptwo.step(bsize);
		color.step(bsize);

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
	octaver plug(Fs);

	automator downone(0.0f); //controls
	automator dry(50.0f);
	automator uphalf(50.0f);
	automator upone(50.0f);
	automator uptwo(50.0f);
	automator color(50.0f);
	if (automation) {
		downone = automator({ 0,		1 * 48000,	2 * 48000 },
			{ 0.0f,		100.0f,		0.0f });
		dry = automator({ 0,		1 * 48000,	2 * 48000 },
			{ 0.0f,		100.0f,		0.0f });
		uphalf = automator({ 0,		1 * 48000,	2 * 48000 },
			{ 0.0f,		100.0f,		0.0f });
		upone = automator({ 0,		1 * 48000,	2 * 48000 },
			{ 0.0f,		100.0f,		0.0f });
		uptwo = automator({ 0,		1 * 48000,	2 * 48000 },
			{ 0.0f,		100.0f,		0.0f });
		color = automator({ 0,		1 * 48000,	2 * 48000 },
			{ 0.0f,		100.0f,		0.0f });
	}
	fvec inbuff(bsize); //I/O
	fvec outbuff(bsize);

	plug.connect_port(0, (void*)(downone.connect())); //ports
	plug.connect_port(1, (void*)(dry.connect()));
	plug.connect_port(2, (void*)(uphalf.connect()));
	plug.connect_port(3, (void*)(upone.connect()));
	plug.connect_port(4, (void*)(uptwo.connect()));
	plug.connect_port(5, (void*)(color.connect()));
	plug.connect_port(6, (void*)inbuff.ptr());
	plug.connect_port(7, (void*)outbuff.ptr());

	plug.init();

	cout << "Loaded plugin, stepping..." << endl;

	// Running
	for (int it = 0; it < iter; ++it) {
		in.get(inbuff, bsize);

		downone.step(bsize);
		dry.step(bsize);
		uphalf.step(bsize);
		upone.step(bsize);
		uptwo.step(bsize);
		color.step(bsize);

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