
//--------------------------------------------------
// Runs tests on the raw plugin object without an
// LV2 wrapper
//--------------------------------------------------

#include <iostream>
#include "cephean-lv2.h"
#include "cephean-test.h"

#include "reverb2.h"

using namespace std;
using namespace cephean;

//==================================================

static void runSoundTest(bool automation)
{
	// File I/O setup ---------------------------
	const double dur = 10.0; //seconds
	const double Fs = 48000.0; //samples/sec
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
	reverb2 plug(Fs);

	automator tail(50.0f); //controls
	automator color(50.0f);
	automator rolloff(0.0f);
	automator energy(100.0f);
	automator wet(50.0f);
	automator dry(50.0f);
	if (automation) {
		//tail = automator(	{ 0,		1 * 48000,	2 * 48000 },
		//						{ 0.0f,		100.0f,		0.0f });
		//color = automator(		{ 0,		1 * 48000,	2 * 48000 },
		//						{ 100.0f,	0.0f,		100.0f });
		//wet = automator(		{ 0,		1 * 48000,	2 * 48000 },
		//						{ 50.0f,	50.0f,		50.0f });
		//dry = automator(		{ 0,		1 * 48000,	2 * 48000 },
		//						{ 50.0f,	50.0f,		50.0f });
	}

	fvec inbuff(bsize); //I/O
	fvec outbuff(bsize);

	plug.connect_port(0, (void*)(tail.connect())); //ports
	plug.connect_port(1, (void*)(color.connect()));
	plug.connect_port(2, (void*)(rolloff.connect()));
	plug.connect_port(3, (void*)(energy.connect()));
	plug.connect_port(4, (void*)(wet.connect()));
	plug.connect_port(5, (void*)(dry.connect()));
	plug.connect_port(6, (void*)inbuff.ptr());
	plug.connect_port(7, (void*)outbuff.ptr());

	plug.init();

	cout << "Loaded plugin, stepping..." << endl;

	// Running
	for (int it = 0; it < iter; ++it) {
		in.get(inbuff, bsize);

		tail.step(bsize);
		color.step(bsize);
		rolloff.step(bsize);
		energy.step(bsize);
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
	reverb2 plug((double)Fs);

	automator tail(50.0f); //controls
	automator color(50.0f);
	automator rolloff(25.0f);
	automator energy(0.0f);
	automator wet(50.0f);
	automator dry(0.0f);

	fvec inbuff(bsize); //I/O
	fvec outbuff(bsize);

	plug.connect_port(0, (void*)(tail.connect())); //ports
	plug.connect_port(1, (void*)(color.connect()));
	plug.connect_port(2, (void*)(rolloff.connect()));
	plug.connect_port(3, (void*)(energy.connect()));
	plug.connect_port(4, (void*)(wet.connect()));
	plug.connect_port(5, (void*)(dry.connect()));
	plug.connect_port(6, (void*)inbuff.ptr());
	plug.connect_port(7, (void*)outbuff.ptr());

	plug.init();

	cout << "Loaded plugin, stepping..." << endl;

	// Running
	for (int it = 0; it < iter; ++it) {
		in.get(inbuff, bsize);

		tail.step(bsize);
		color.step(bsize);
		rolloff.step(bsize);
		energy.step(bsize);
		wet.step(bsize);
		dry.step(bsize);

		plug.step(bsize);

		out.put(outbuff, bsize);
	}

	//monitor DC offset and total power
	float DC = sum(out.data(), out.size()) / out.size();
	float POW = ssq(out.data(), out.size());

	//Run a smoothed power detector on the result
	vdiv(out.data(), rms(out.data(), npts), out.data(), npts); //power normalize
	fvec xpow = plotSmoothedPower(out.data(), npts, 0.01f * Fs);

	//Run a magnitude FFT on the result
	fvec window(2 * npts); //apply Hann window for better spectral resolution
	hann(window.ptr(), window.size());
	vmult(out.data(), window.ptr() + npts, out.data(), npts);
	fvec X = plotLogSmoothedDBSpectrum(out.data(), npts,
		20.0f / (Fs / 2.0f), 20000.0f / (Fs / 2.0f), 1 / 12.0f);

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
	reverb2 plug(Fs);

	automator tail(50.0f); //controls
	automator color(50.0f);
	automator rolloff(0.0f);
	automator energy(50.0f);
	automator wet(50.0f);
	automator dry(50.0f);
	if (automation) {
		//tail = automator(	{ 0,		1 * 48000,	2 * 48000 },
		//						{ 0.0f,		100.0f,		0.0f });
		//color = automator(		{ 0,		1 * 48000,	2 * 48000 },
		//						{ 100.0f,	0.0f,		100.0f });
		//wet = automator(		{ 0,		1 * 48000,	2 * 48000 },
		//						{ 50.0f,	50.0f,		50.0f });
		//dry = automator(		{ 0,		1 * 48000,	2 * 48000 },
		//						{ 50.0f,	50.0f,		50.0f });
	}

	fvec inbuff(bsize); //I/O
	fvec outbuff(bsize);

	plug.connect_port(0, (void*)(tail.connect())); //ports
	plug.connect_port(1, (void*)(color.connect()));
	plug.connect_port(2, (void*)(rolloff.connect()));
	plug.connect_port(3, (void*)(energy.connect()));
	plug.connect_port(4, (void*)(wet.connect()));
	plug.connect_port(5, (void*)(dry.connect()));
	plug.connect_port(6, (void*)inbuff.ptr());
	plug.connect_port(7, (void*)outbuff.ptr());

	plug.init();

	cout << "Loaded plugin, stepping..." << endl;

	// Running
	for (int it = 0; it < iter; ++it) {
		in.get(inbuff, bsize);

		tail.step(bsize);
		color.step(bsize);
		rolloff.step(bsize);
		energy.step(bsize);
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
	//runImpulseTest();
	runStressTest(false);

	return 0;
}