
//--------------------------------------------------
// Runs tests on the raw plugin object without an
// LV2 wrapper
//--------------------------------------------------

#include <iostream>
#include "cephean-lv2.h"
#include "cephean-test.h"

#include "gain.h"

using namespace std;
using namespace cephean;

int main()
{
	//TODO -- make wavreader capable of looping for stress test profiling
	//	-and add a bunch of test features like plotting helpers, control automation, etc.
	
	// I/O setup
	const double dur = 10.0; //seconds
	wavreader fr("C:\\Users\\crees\\Music\\Cephean tests\\doo_riffing.wav", dur, 3.0);
	double Fs = fr.getSampleRate();

	wavwriter fw("C:\\Users\\crees\\Music\\Cephean tests\\test_output.wav", 1, Fs, dur);

	const int bsize = 128; //samples per block
	const int iter = (int)ceil(10.0 * Fs / bsize); // to fill the duration

	// Plugin setup
	gain plug(Fs);
	float ctrl_level = -20.0f; //dB
	fvec inbuff(bsize);
	fvec outbuff(bsize);
	plug.connect_port(0, (void*)(&ctrl_level));
	plug.connect_port(1, (void*)inbuff.ptr());
	plug.connect_port(2, (void*)outbuff.ptr());
	plug.init();
	cout << "Loaded plugin, stepping..." << endl;
	
	// Running
	for (int it = 0; it < iter; ++it) {
		fr.get(inbuff, bsize);
		plug.step((uint32_t)bsize);
		fw.put(outbuff, bsize);
	}
	cout << "done." << endl;

	//Export the audio file
	fw.write();
	
}