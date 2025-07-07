
//--------------------------------------------------
// Runs tests on the raw plugin object without an
// LV2 wrapper
//--------------------------------------------------

#include <iostream>
#include "cephean-lv2.h"
#include "cephean-test.h"

#include "reverb.h"

//LEFT OFF HERE <--- setting up my project environment
//
//	It seems easy to copy/paste an entire project folder, then rename things.
//	This will make it easy to write Python script to copy from a templateplug
//	project and find/replace all file names, and inside all source and project
//	files! Next step is to write that utility.

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

	wavreader in("C:\\Users\\crees\\Music\\Cephean tests\\doo_riffing.wav", dur, 5.0);
	assert(in.getSampleRate() == Fs);
	wavwriter out("C:\\Users\\crees\\Music\\Cephean tests\\test_output.wav", 1, Fs, dur);

	gain outlevel;
	outlevel.setLevel(1.0f, true);

	// Plugin setup -----------------------------
	//all seeds: { 488, 845, 268, 390, 943, 262, 849, 724, 531, 899, 652, 651, 403, 586, 792, 329, 969, 829, 636, 320, 326, 314, 236, 667 }
	//final: { 652, 403, 329, 488, 845, 268, 390, 943, 262, 849, 724, 531, 651, 792, 969, 829 }
	reverb plug((double)Fs);

	iautomator seed(0); //controls
	automator duration(50.0f);
	automator color(100.0f);
	automator wet(50.0f);
	automator dry(0.0f);
	if (automation) {
		//seed = iautomator({ 0,		1 * 48000,	2 * 48000 },
		//	{ 0, 1, 2 });
		//duration = automator(	{ 0,		1 * 48000,	2 * 48000 },
		//						{ 0.0f,		100.0f,		0.0f });
		color = automator(		{ 0,		1 * 48000,	2 * 48000 },
								{ 100.0f,	0.0f,		100.0f });
		//wet = automator(		{ 0,		1 * 48000,	2 * 48000 },
		//						{ 50.0f,	50.0f,		50.0f });
		//dry = automator(		{ 0,		1 * 48000,	2 * 48000 },
		//						{ 50.0f,	50.0f,		50.0f });
	}

	fvec inbuff(bsize); //I/O
	fvec outbuff(bsize);

	plug.connect_port(0, (void*)(seed.connect())); //ports
	plug.connect_port(1, (void*)(duration.connect()));
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
		duration.step(bsize);
		color.step(bsize);
		wet.step(bsize);
		dry.step(bsize);

		plug.step(bsize);

		outlevel.stepBlock(outbuff, outbuff, bsize);
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
	//seeds: { 488, 845, 268, 390, 943, 262, 849, 724, 531, 899, 652, 651, 403, 586, 792, 329, 969, 829, 636, 320, 326, 314, 236, 667 }
	reverb plug((double)Fs);

	iautomator seed(0); //controls
	automator duration(100.0f);
	automator color(100.0f);
	automator wet(50.0f);
	automator dry(0.0f);

	fvec inbuff(bsize); //I/O
	fvec outbuff(bsize);

	plug.connect_port(0, (void*)(seed.connect())); //ports
	plug.connect_port(1, (void*)(duration.connect()));
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

		duration.step(bsize);
		color.step(bsize);
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
	reverb plug(Fs);

	iautomator seed(0); //controls
	automator duration(50.0f);
	automator color(50.0f);
	automator wet(50.0f);
	automator dry(50.0f);
	if (automation) {
		//seed = iautomator({ 0,		1 * 48000,	2 * 48000 },
		//	{ 0, 1, 2 });
		duration = automator(	{ 0,		1 * 48000,	2 * 48000 },
								{ 0.0f,		100.0f,		0.0f });
		color = automator(		{ 0,		1 * 48000,	2 * 48000 },
								{ 0.0f,		100.0f,		0.0f });
		wet = automator(		{ 0,		1 * 48000,	2 * 48000 },
								{ 0.0f,		100.0f,		0.0f });
		dry = automator(		{ 0,		1 * 48000,	2 * 48000 },
								{ 0.0f,		100.0f,		0.0f });
	}

	fvec inbuff(bsize); //I/O
	fvec outbuff(bsize);

	plug.connect_port(0, (void*)(seed.connect())); //ports
	plug.connect_port(1, (void*)(duration.connect()));
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

		duration.step(bsize);
		color.step(bsize);
		wet.step(bsize);
		dry.step(bsize);

		plug.step(bsize);
	}

	cout << "done." << endl;
}

//==================================================

static void determineBestSeed()
{
	//Global setup -----------------------------
	const float Fs = 48000.0;
	const int bsize = 128;
	const int npts = 65536;
	const int iter = npts / bsize;

	fvec window(2 * npts); //apply Hann window for smoother results
	hann(window.ptr(), window.size());
	const float Noct = 1 / 12.0f; //resolution of peak finder

	impulse in;
	testbuff out(npts);

	iautomator seed(0); //controls
	automator duration(100.0f);
	automator color(100.0f);
	automator wet(50.0f);
	automator dry(0.0f);

	fvec inbuff(bsize); //I/O
	fvec outbuff(bsize);

	cout << "Entering seed test loop..." << endl;

	int testlen = 1000; //how many seeds to check
	fvec tpeaks(testlen);
	fvec fpeaks(testlen);
	for (int testi = 0; testi < testlen; ++testi) {

		// Plugin setup -----------------------------
		reverb plug((double)Fs, testi);

		plug.connect_port(0, (void*)(seed.connect())); //ports
		plug.connect_port(1, (void*)(duration.connect()));
		plug.connect_port(2, (void*)(color.connect()));
		plug.connect_port(3, (void*)(wet.connect()));
		plug.connect_port(4, (void*)(dry.connect()));
		plug.connect_port(5, (void*)inbuff.ptr());
		plug.connect_port(6, (void*)outbuff.ptr());

		plug.init();
		in.reset();
		out.clear();

		// Running
		for (int it = 0; it < iter; ++it) {
			in.get(inbuff, bsize);

			seed.step(bsize);
			duration.step(bsize);
			color.step(bsize);
			wet.step(bsize);
			dry.step(bsize);

			plug.step(bsize);

			out.put(outbuff, bsize);
		}

		//Run a smoothed power detector on the result
		vdiv(out.data(), rms(out.data(), npts), out.data(), npts); //power normalize
		fvec xpow = plotSmoothedPower(out.data(), npts, 0.01f * Fs);

		//Run a windowed magnitude FFT on the result
		vmult(out.data(), window.ptr() + npts, out.data(), npts);
		fvec X = plotLogSmoothedDBSpectrum(out.data(), npts,
			20.0f / (Fs / 2.0f), 20000.0f / (Fs / 2.0f), Noct);

		//find the largest subband and time peak in this verb tail
		tpeaks[testi] = vfindmax(xpow.ptr(), xpow.size());
		fpeaks[testi] = vfindmax(X.ptr(), X.size());
	}

	//Post process to find the best combinations of time/frequency peakiness
	fvec best_tpeaks(testlen);
	ivec ind_tpeaks(testlen);
	ivec tscores(testlen);
	vsortmin(tpeaks.ptr(), best_tpeaks.ptr(), ind_tpeaks.ptr(), testlen);
	unsort(ind_tpeaks.ptr(), tscores.ptr(), testlen);

	fvec best_fpeaks(testlen);
	ivec ind_fpeaks(testlen);
	ivec fscores(testlen);
	vsortmin(fpeaks.ptr(), best_fpeaks.ptr(), ind_fpeaks.ptr(), testlen);
	unsort(ind_fpeaks.ptr(), fscores.ptr(), testlen);

	ivec score(testlen);
	ivec best_score(testlen);
	ivec ind_score(testlen);
	vadd(tscores.ptr(), fscores.ptr(), score.ptr(), testlen);
	vsortmin(score.ptr(), best_score.ptr(), ind_score.ptr(), testlen);

	int best_seeds[24];
	vcopy(ind_score.ptr(), best_seeds, 24);

	cout << "done." << endl;
}

//==================================================

int main()
{
	runSoundTest(true);
	//runImpulseTest();
	//runStressTest(false);
	//determineBestSeed();

	return 0;
}