
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

	automator band1(50.0f); //controls
	automator band2(50.0f);
	automator band3(50.0f);
	automator band4(50.0f);
	automator band5(50.0f);
	automator band6(50.0f);
	automator band7(50.0f);
	automator band8(50.0f);
	automator band9(50.0f);
	if (automation) {
		band6 = automator({ 0,		5 * 48000,	10 * 48000 },
			{ 0.0f,		100.0f,		0.0f });
		/*
		band2 = automator({ 0,		1 * 48000,	2 * 48000 },
			{ 0.0f,		100.0f,		0.0f });
		band3 = automator({ 0,		1 * 48000,	2 * 48000 },
			{ 0.0f,		100.0f,		0.0f });
		band4 = automator({ 0,		1 * 48000,	2 * 48000 },
			{ 0.0f,		100.0f,		0.0f });
		band5 = automator({ 0,		1 * 48000,	2 * 48000 },
			{ 0.0f,		100.0f,		0.0f });
		band6 = automator({ 0,		1 * 48000,	2 * 48000 },
			{ 0.0f,		100.0f,		0.0f });
		band7 = automator({ 0,		1 * 48000,	2 * 48000 },
			{ 0.0f,		100.0f,		0.0f });
		band8 = automator({ 0,		1 * 48000,	2 * 48000 },
			{ 0.0f,		100.0f,		0.0f });
		band9 = automator({ 0,		1 * 48000,	2 * 48000 },
			{ 0.0f,		100.0f,		0.0f });
			*/
	}
	fvec inbuff(bsize); //I/O
	fvec outbuff(bsize);

	plug.connect_port(0, (void*)(band1.connect())); //ports
	plug.connect_port(1, (void*)(band2.connect()));
	plug.connect_port(2, (void*)(band3.connect()));
	plug.connect_port(3, (void*)(band4.connect()));
	plug.connect_port(4, (void*)(band5.connect()));
	plug.connect_port(5, (void*)(band6.connect()));
	plug.connect_port(6, (void*)(band7.connect()));
	plug.connect_port(7, (void*)(band8.connect()));
	plug.connect_port(8, (void*)(band9.connect()));
	plug.connect_port(9, (void*)inbuff.ptr());
	plug.connect_port(10, (void*)outbuff.ptr());

	plug.init();

	cout << "Loaded plugin, stepping..." << endl;

	// Running
	for (int it = 0; it < iter; ++it) {
		in.get(inbuff, bsize);

		band1.step(bsize);
		band2.step(bsize);
		band3.step(bsize);
		band4.step(bsize);
		band5.step(bsize);
		band6.step(bsize);
		band7.step(bsize);
		band8.step(bsize);
		band9.step(bsize);

		plug.step(bsize);

		out.put(outbuff, bsize);
	}

	//Export the audio file
	out.write();

	cout << "done." << endl;
}

//==================================================

static void runImpulseTest()
{
	//Impulse setup -----------------------------
	const float Fs = 48000.0;
	const int bsize = 128;
	const int npts = 65536;
	const int iter = npts / bsize;

	impulse in;
	testbuff out(npts);

	// Plugin setup -----------------------------
	graphiceq plug((double)Fs);

	automator band1(50.0f); //controls
	automator band2(50.0f);
	automator band3(50.0f);
	automator band4(50.0f);
	automator band5(100.0f);
	automator band6(100.0f);
	automator band7(50.0f);
	automator band8(50.0f);
	automator band9(50.0f);

	fvec inbuff(bsize); //I/O
	fvec outbuff(bsize);

	plug.connect_port(0, (void*)(band1.connect())); //ports
	plug.connect_port(1, (void*)(band2.connect()));
	plug.connect_port(2, (void*)(band3.connect()));
	plug.connect_port(3, (void*)(band4.connect()));
	plug.connect_port(4, (void*)(band5.connect()));
	plug.connect_port(5, (void*)(band6.connect()));
	plug.connect_port(6, (void*)(band7.connect()));
	plug.connect_port(7, (void*)(band8.connect()));
	plug.connect_port(8, (void*)(band9.connect()));
	plug.connect_port(9, (void*)inbuff.ptr());
	plug.connect_port(10, (void*)outbuff.ptr());

	plug.init();

	cout << "Loaded plugin, stepping..." << endl;

	// Running
	for (int it = 0; it < iter; ++it) {
		in.get(inbuff, bsize);

		band1.step(bsize);
		band2.step(bsize);
		band3.step(bsize);
		band4.step(bsize);
		band5.step(bsize);
		band6.step(bsize);
		band7.step(bsize);
		band8.step(bsize);
		band9.step(bsize);

		plug.step(bsize);

		out.put(outbuff, bsize);
	}

	//Run a magnitude FFT on the result
	fvec X = plotLogSmoothedDBSpectrum(out.data(), npts,
		20.0f / (Fs / 2.0f), 20000.0f / (Fs / 2.0f), 1 / 48.0f);

	cout << "done." << endl;
}

//==================================================

static void runStressTest(bool automation)
{
	// File I/O setup ---------------------------
	const double dur = 1000.0; //seconds, looped
	const double Fs = 48000.0; //samples/sec
	const int bsize = 128; //samples per block
	const int iter = (int)ceil(dur * Fs / bsize); // to fill the duration

	wavreader in("C:\\Users\\crees\\Music\\Cephean tests\\doo_riffing.wav", dur, 5.0, 25.0);
	assert(in.getSampleRate() == Fs);

	// Plugin setup -----------------------------
	graphiceq plug(Fs);

	automator band1(50.0f); //controls
	automator band2(50.0f);
	automator band3(50.0f);
	automator band4(50.0f);
	automator band5(50.0f);
	automator band6(50.0f);
	automator band7(50.0f);
	automator band8(50.0f);
	automator band9(50.0f);
	if (automation) {
		band1 = automator({ 0,		5 * 48000,	10 * 48000 },
			{ 0.0f,		100.0f,		0.0f });
		band2 = automator({ 0,		1 * 48000,	2 * 48000 },
			{ 0.0f,		100.0f,		0.0f });
		band3 = automator({ 0,		1 * 48000,	2 * 48000 },
			{ 0.0f,		100.0f,		0.0f });
		band4 = automator({ 0,		1 * 48000,	2 * 48000 },
			{ 0.0f,		100.0f,		0.0f });
		band5 = automator({ 0,		1 * 48000,	2 * 48000 },
			{ 0.0f,		100.0f,		0.0f });
		band6 = automator({ 0,		1 * 48000,	2 * 48000 },
			{ 0.0f,		100.0f,		0.0f });
		band7 = automator({ 0,		1 * 48000,	2 * 48000 },
			{ 0.0f,		100.0f,		0.0f });
		band8 = automator({ 0,		1 * 48000,	2 * 48000 },
			{ 0.0f,		100.0f,		0.0f });
		band9 = automator({ 0,		1 * 48000,	2 * 48000 },
			{ 0.0f,		100.0f,		0.0f });
	}
	fvec inbuff(bsize); //I/O
	fvec outbuff(bsize);

	plug.connect_port(0, (void*)(band1.connect())); //ports
	plug.connect_port(1, (void*)(band2.connect()));
	plug.connect_port(2, (void*)(band3.connect()));
	plug.connect_port(3, (void*)(band4.connect()));
	plug.connect_port(4, (void*)(band5.connect()));
	plug.connect_port(5, (void*)(band6.connect()));
	plug.connect_port(6, (void*)(band7.connect()));
	plug.connect_port(7, (void*)(band8.connect()));
	plug.connect_port(8, (void*)(band9.connect()));
	plug.connect_port(9, (void*)inbuff.ptr());
	plug.connect_port(10, (void*)outbuff.ptr());

	plug.init();

	cout << "Loaded plugin, stepping..." << endl;

	// Running
	for (int it = 0; it < iter; ++it) {
		in.get(inbuff, bsize);

		band1.step(bsize);
		band2.step(bsize);
		band3.step(bsize);
		band4.step(bsize);
		band5.step(bsize);
		band6.step(bsize);
		band7.step(bsize);
		band8.step(bsize);
		band9.step(bsize);

		plug.step(bsize);
	}

	cout << "done." << endl;
}

//==================================================

int main()
{
	//runSoundTest(true);
	//runStressTest(true);
	runImpulseTest();

	return 0;
}