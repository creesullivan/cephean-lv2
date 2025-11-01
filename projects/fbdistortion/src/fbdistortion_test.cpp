
//--------------------------------------------------
// Runs tests on the raw plugin object without an
// LV2 wrapper
//--------------------------------------------------

#include <iostream>
#include "cephean-lv2.h"
#include "cephean-test.h"

#include "fbdistortion.h"

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
	wavreader in("C:\\Users\\crees\\Music\\Cephean tests\\long_guitar_chords.wav", dur, 17.0);
	//wavreader in("C:\\Users\\crees\\Music\\Cephean tests\\finger_bass_riff.wav", dur, 25.0);
	assert(in.getSampleRate() == Fs);
	wavwriter out("C:\\Users\\crees\\Music\\Cephean tests\\test_output.wav", 1, Fs, dur);

	// Plugin setup -----------------------------
	fbdistortion plug(Fs);

	automator drive(80.0f); //controls
	automator shape(80.0f);
	automator feedback(50.0f);
	automator dirty(50.0f);
	automator color(50.0f);
	automator level(70.0f);
	if (automation) {
		drive = automator({ 0,		1 * 48000,	2 * 48000 },
			{ 80.0f,		80.0f,		80.0f });
		shape = automator({ 0,		1 * 48000,	2 * 48000 },
			{ 80.0f,		80.0f,		80.0f });
		feedback = automator({ 0,		1 * 48000,	2 * 48000 },
			{ 0.0f,		0.0f,		0.0f });
		dirty = automator({ 0,		1 * 48000,	2 * 48000 },
			{ 50.0f,		50.0f,		50.0f });
		color = automator({ 0,		1 * 48000,	2 * 48000 },
			{ 0.0f,		100.0f,		0.0f });
		level = automator({ 0,		1 * 48000,	2 * 48000 },
			{ 70.0f,		70.0f,		70.0f });
	}
	fvec inbuff(bsize); //I/O
	fvec scbuff(bsize);
	fvec outbuff(bsize);

	plug.connect_port(0, (void*)(drive.connect())); //ports
	plug.connect_port(1, (void*)(shape.connect()));
	plug.connect_port(2, (void*)(feedback.connect()));
	plug.connect_port(3, (void*)(dirty.connect()));
	plug.connect_port(4, (void*)(color.connect()));
	plug.connect_port(5, (void*)(level.connect()));
	plug.connect_port(6, (void*)inbuff.ptr());
	plug.connect_port(7, (void*)scbuff.ptr());
	plug.connect_port(8, (void*)outbuff.ptr());

	plug.init();

	cout << "Loaded plugin, stepping..." << endl;

	// Running
	float maxval = 0.0f;
	for (int it = 0; it < iter; ++it) {
		in.get(inbuff, bsize);
		vcopy(inbuff.ptr(), scbuff.ptr(), bsize); //trigger on the same input signal to test

		drive.step(bsize);
		feedback.step(bsize);
		shape.step(bsize);
		dirty.step(bsize);
		color.step(bsize);
		level.step(bsize);

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

static void runStressTest(bool automation)
{
	// File I/O setup ---------------------------
	const double dur = 1000.0; //seconds
	const float Fs = 48000.0f; //samples/sec
	const int bsize = 128; //samples per block
	const int iter = (int)ceil(dur * Fs / bsize); // to fill the duration

	//wavreader in("C:\\Users\\crees\\Music\\Cephean tests\\doo_riffing.wav", dur, 5.0);
	//wavreader in("C:\\Users\\crees\\Music\\Cephean tests\\guitar_chugs.wav", dur, 5.0);
	//wavreader in("C:\\Users\\crees\\Music\\Cephean tests\\dry_guitar_lead.wav", dur, 48.0);
	wavreader in("C:\\Users\\crees\\Music\\Cephean tests\\long_guitar_chords.wav", dur, 17.0);
	//wavreader in("C:\\Users\\crees\\Music\\Cephean tests\\finger_bass_riff.wav", dur, 25.0);
	assert(in.getSampleRate() == Fs);

	// Plugin setup -----------------------------
	fbdistortion plug(Fs);

	automator drive(80.0f); //controls
	automator shape(80.0f);
	automator feedback(0.0f);
	automator dirty(50.0f);
	automator color(70.0f);
	automator level(70.0f);
	if (automation) {
		/* TODO
		level = automator(		{ 0,		1 * 48000,	2 * 48000 },
								{ 0.0f,		100.0f,		0.0f });
								*/
	}
	fvec inbuff(bsize); //I/O
	fvec scbuff(bsize);
	fvec outbuff(bsize);

	plug.connect_port(0, (void*)(drive.connect())); //ports
	plug.connect_port(1, (void*)(shape.connect()));
	plug.connect_port(2, (void*)(feedback.connect()));
	plug.connect_port(3, (void*)(dirty.connect()));
	plug.connect_port(4, (void*)(color.connect()));
	plug.connect_port(5, (void*)(level.connect()));
	plug.connect_port(6, (void*)inbuff.ptr());
	plug.connect_port(7, (void*)scbuff.ptr());
	plug.connect_port(8, (void*)outbuff.ptr());

	plug.init();

	cout << "Loaded plugin, stepping..." << endl;

	// Running
	float maxval = 0.0f;
	for (int it = 0; it < iter; ++it) {
		in.get(inbuff, bsize);
		vcopy(inbuff.ptr(), scbuff.ptr(), bsize); //trigger on the same input signal to test

		drive.step(bsize);
		feedback.step(bsize);
		shape.step(bsize);
		dirty.step(bsize);
		color.step(bsize);
		level.step(bsize);

		plug.step(bsize);
	}

	cout << "done." << endl;
}


//==================================================

int main()
{
	runSoundTest(false);
	//runStressTest(false);

	//NOTES WHILE PLAYING LIVE:
	//	-feedback is good enough now as is =)
	//  -the distortion is too! sounds utterly amazing on light, medium, and heavy tones =)
	//	-take a hard look at the best range for the feedback parameter, drive dependent?
	//	-the buffer edge has a hard time ringing out, maybe investigate, but this is one note and isn't a big issue
	//		-around F#, kinda annoying...
	//
	//TODO ITEMS AFTER THE BUGS ARE SQUASHED:
	//	-I want to have a sidechain detect input for the string model

	return 0;
}