
//--------------------------------------------------
// Runs tests on the raw plugin object without an
// LV2 wrapper
//--------------------------------------------------

#include <iostream>
#include "cephean-lv2.h"
#include "cephean-test.h"

#include "fbdelay.h"

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
	//wavreader in("C:\\Users\\crees\\Music\\Cephean tests\\guitar_chugs.wav", dur, 5.0);
	//wavreader in("C:\\Users\\crees\\Music\\Cephean tests\\dry_guitar_lead.wav", dur, 48.0);
	//wavreader in("C:\\Users\\crees\\Music\\Cephean tests\\long_guitar_chords.wav", dur, 17.0);
	//wavreader in("C:\\Users\\crees\\Music\\Cephean tests\\finger_bass_riff.wav", dur, 25.0);
	assert(in.getSampleRate() == Fs);
	wavwriter out("C:\\Users\\crees\\Music\\Cephean tests\\test_output.wav", 1, Fs, dur);

	// Plugin setup -----------------------------
	fbdelay plug(Fs);

	automator time(50.0f); //controls
	automator feedback(70.0f);
	automator wet(100.0f);
	automator spread(0.0f);
	automator low(0.0f);
	automator high(50.0f);
	if (automation) {
		//time = automator({ 0,		5 * 48000,	10 * 48000 },
		//	{ 0.0f,		100.0f,		0.0f });
		//feedback = automator({ 0,		5 * 48000,	10 * 48000 },
		//	{ 0.0f,		100.0f,		0.0f });
		//wet = automator({ 0,		5 * 48000,	10 * 48000 },
		//	{ 0.0f,		100.0f,		0.0f });
		spread = automator({ 0,		5 * 48000,	10 * 48000 },
			{ 0.0f,		100.0f,		0.0f });
		//low = automator({ 0,		5 * 48000,	10 * 48000 },
		//	{ 100.0f,		0.0f,		100.0f });
		//high = automator({ 0,		5 * 48000,	10 * 48000 },
		//	{ 0.0f,		100.0f,		0.0f });
	}
	fvec inbuff(bsize); //I/O
	fvec outbuff(bsize);

	plug.connect_port(0, (void*)(time.connect())); //ports
	plug.connect_port(1, (void*)(feedback.connect()));
	plug.connect_port(2, (void*)(wet.connect()));
	plug.connect_port(3, (void*)(spread.connect()));
	plug.connect_port(4, (void*)(low.connect()));
	plug.connect_port(5, (void*)(high.connect()));
	plug.connect_port(6, (void*)inbuff.ptr());
	plug.connect_port(7, (void*)outbuff.ptr());

	plug.init();

	cout << "Loaded plugin, stepping..." << endl;

	// Running
	for (int it = 0; it < iter; ++it) {
		in.get(inbuff, bsize);

		time.step(bsize);
		feedback.step(bsize);
		wet.step(bsize);
		spread.step(bsize);
		low.step(bsize);
		high.step(bsize);

		if (it == (iter / 2)) {
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
	//{ 689, 353, 90, 387, 24, 668, 126, 135, 516, 938, 341, 202, 293, 277, 99, 709 };
	fbdelay plug((double)Fs, 709);

	automator time(50.0f); //controls
	automator feedback(50.0f);
	automator wet(50.0f);
	automator spread(100.0f);
	automator low(50.0f);
	automator high(50.0f);

	fvec inbuff(bsize); //I/O
	fvec outbuff(bsize);

	plug.connect_port(0, (void*)(time.connect())); //ports
	plug.connect_port(1, (void*)(feedback.connect()));
	plug.connect_port(2, (void*)(wet.connect()));
	plug.connect_port(3, (void*)(spread.connect()));
	plug.connect_port(4, (void*)(low.connect()));
	plug.connect_port(5, (void*)(high.connect()));
	plug.connect_port(6, (void*)inbuff.ptr());
	plug.connect_port(7, (void*)outbuff.ptr());

	plug.init();

	cout << "Loaded plugin, stepping..." << endl;

	// Running
	for (int it = 0; it < iter; ++it) {
		in.get(inbuff, bsize);

		time.step(bsize);
		feedback.step(bsize);
		wet.step(bsize);
		spread.step(bsize);
		low.step(bsize);
		high.step(bsize);

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
	fbdelay plug(Fs);

	automator time(50.0f); //controls
	automator feedback(70.0f);
	automator wet(70.0f);
	automator spread(50.0f);
	automator low(0.0f);
	automator high(50.0f);
	if (automation) {
		//time = automator({ 0,		5 * 48000,	10 * 48000 },
		//	{ 0.0f,		100.0f,		0.0f });
		//feedback = automator({ 0,		5 * 48000,	10 * 48000 },
		//	{ 0.0f,		100.0f,		0.0f });
		//wet = automator({ 0,		5 * 48000,	10 * 48000 },
		//	{ 0.0f,		100.0f,		0.0f });
		spread = automator({ 0,		5 * 48000,	10 * 48000 },
			{ 0.0f,		100.0f,		0.0f });
		//low = automator({ 0,		5 * 48000,	10 * 48000 },
		//	{ 100.0f,		0.0f,		100.0f });
		//high = automator({ 0,		5 * 48000,	10 * 48000 },
		//	{ 0.0f,		100.0f,		0.0f });
	}
	fvec inbuff(bsize); //I/O
	fvec outbuff(bsize);

	plug.connect_port(0, (void*)(time.connect())); //ports
	plug.connect_port(1, (void*)(feedback.connect()));
	plug.connect_port(2, (void*)(wet.connect()));
	plug.connect_port(3, (void*)(spread.connect()));
	plug.connect_port(4, (void*)(low.connect()));
	plug.connect_port(5, (void*)(high.connect()));
	plug.connect_port(6, (void*)inbuff.ptr());
	plug.connect_port(7, (void*)outbuff.ptr());

	plug.init();

	cout << "Loaded plugin, stepping..." << endl;

	// Running
	for (int it = 0; it < iter; ++it) {
		in.get(inbuff, bsize);

		time.step(bsize);
		feedback.step(bsize);
		wet.step(bsize);
		spread.step(bsize);
		low.step(bsize);
		high.step(bsize);

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

	automator time(50.0f); //controls -- unused when testing echo
	automator feedback(50.0f);
	automator wet(50.0f);
	automator spread(100.0f);
	automator low(50.0f);
	automator high(50.0f);

	fvec inbuff(bsize); //I/O
	fvec outbuff(bsize);

	cout << "Entering seed test loop..." << endl;

	int testlen = 1000; //how many seeds to check
	fvec tpeaks(testlen);
	fvec fpeaks(testlen);
	for (int testi = 0; testi < testlen; ++testi) {

		// Plugin setup -----------------------------
		fbdelay plug((double)Fs, testi);

		plug.connect_port(0, (void*)(time.connect())); //ports
		plug.connect_port(1, (void*)(feedback.connect()));
		plug.connect_port(2, (void*)(wet.connect()));
		plug.connect_port(3, (void*)(spread.connect()));
		plug.connect_port(4, (void*)(low.connect()));
		plug.connect_port(5, (void*)(high.connect()));
		plug.connect_port(6, (void*)inbuff.ptr());
		plug.connect_port(7, (void*)outbuff.ptr());

		plug.init();
		in.reset();
		out.clear();

		// Running
		for (int it = 0; it < iter; ++it) {
			in.get(inbuff, bsize);

			time.step(bsize);
			feedback.step(bsize);
			wet.step(bsize);
			spread.step(bsize);
			low.step(bsize);
			high.step(bsize);

			plug.step(bsize);

			out.put(outbuff, bsize);
		}

		//Run a smoothed power detector on the result
		vdiv(out.data(), rms(out.data(), npts), out.data(), npts); //power normalize
		fvec xpow = plotSmoothedPower(out.data(), npts, 0.01f * Fs);

		//Run a windowed magnitude FFT on the result
		vmult(out.data(), window.ptr() + npts, out.data(), npts);
		fvec X = plotLogSmoothedDBSpectrum(out.data(), npts,
			20.0f / (Fs / 2.0f), 20000.0f / (Fs / 2.0f), 1 / 12.0f);

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

	cout << "done." << endl;
}


//==================================================

static void fitPowerCurves()
{
	cout << "Manually fitting echo power curve to data..." << endl;

	//Power values gathered from total SSQ of impulse responses
	//at each echo duration with the default seed before applying
	//power normalization to the wet gain.
	const float durs[11] = { 0.0f, 0.025f, 0.05f, 0.075f, 0.1f, 0.125f, 0.15f, 0.175f, 0.2f, 0.225f, 0.25f };
	float pows[11] = { 0.0f, 0.00353f, 0.01414f, 0.02662f, 0.04087f, 0.05735f,
		0.07610f, 0.09690f, 0.11951f, 0.14367f, 0.16917f };
	const float pow0 = 0.01f; //fixed dry mix with 0.1 gain multiplier

	vadd(pows, pow0, pows, 11); //add power from dry mix

	float approx[11];
	float err[11];
	const float r = 1.6f;
	const float g = (pows[10] - pow0) / powf(durs[10], r);
	for (int i = 0; i < 11; ++i) {
		approx[i] = powf(durs[i], r) * g + pow0; //super-linear approx
		err[i] = 10.0f * log10f(pows[i] / approx[i]); //dB error
	}

	//SHIP IT -- <0.25 dB error over all values

	cout << "done." << endl;
}


//==================================================

static void fitAveragePowerCurve()
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

	float time = 50.0f; //controls -- only spread used when testing echo
	float feedback = 50.0f;
	float wet = 50.0f;
	float spread = 50.0f;
	float low = 50.0f;
	float high = 50.0f;

	fvec inbuff(bsize); //I/O
	fvec outbuff(bsize);

	cout << "Collecting power curves for each echo seed..." << endl;

	int seeds[16] = { 689, 353, 90, 387, 24, 668, 126, 135, 516, 938, 341, 202, 293, 277, 99, 709 };
	float spreads[11] = { 0.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.0f };
	float powers[16][11]; //detected total powers of each IR

	for (int ind1 = 0; ind1 < 16; ++ind1) {

		// Plugin setup -----------------------------
		fbdelay plug((double)Fs, seeds[ind1]);

		plug.connect_port(0, (void*)(&time)); //ports
		plug.connect_port(1, (void*)(&feedback));
		plug.connect_port(2, (void*)(&wet));
		plug.connect_port(3, (void*)(&spread));
		plug.connect_port(4, (void*)(&low));
		plug.connect_port(5, (void*)(&high));
		plug.connect_port(6, (void*)inbuff.ptr());
		plug.connect_port(7, (void*)outbuff.ptr());

		for (int ind2 = 0; ind2 < 11; ++ind2) {

			spread = spreads[ind2];

			plug.init();
			in.reset();
			out.clear();

			// Running
			for (int it = 0; it < iter; ++it) {
				in.get(inbuff, bsize);

				plug.step(bsize);

				out.put(outbuff, bsize);
			}

			// Compute the power of the result and assign it
			powers[ind1][ind2] = ssq(out.data(), out.size());
		}

	}

	cout << "done." << endl;

	cout << "Manually fitting echo power curve to data..." << endl;

	const float durs[11] = { 0.0f, 0.025f, 0.05f, 0.075f, 0.1f, 0.125f, 0.15f, 0.175f, 0.2f, 0.225f, 0.25f };
	const float pow0 = 0.01f; //fixed dry mix with 0.1 gain multiplier

	float approx[16][11];
	float err[11][16];
	float worsterr[11];
	const float r = 1.7f;
	float g[16];
	for (int j = 0; j < 11; ++j) {
		for (int i = 0; i < 16; ++i) {
			g[i] = powers[i][10] / powf(durs[10], r);
			approx[i][j] = powf(durs[j], r) * g[i] + pow0; //super-linear approx
			powers[i][j] += pow0; //include bias dry power
			err[j][i] = 10.0f * fabsf(log10f(powers[i][j] / approx[i][j])); //dB error
		}
		worsterr[j] = vfindmax(err[j], 16); //consider worst case errors
	}

	cout << "done." << endl;
}


//==================================================

int main()
{
	runSoundTest(false);
	//runImpulseTest();
	//runStressTest(false);
	//determineBestSeed();
	//fitPowerCurves();
	//fitAveragePowerCurve();

	//LEFT OFF HERE <-- everything is tested and working, time to export!

	return 0;
}