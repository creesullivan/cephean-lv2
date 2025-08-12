
#pragma once

//--------------------------------------------------
// Defines the module's audio processing behavior
//--------------------------------------------------

#include "cephean-lv2.h"
#include "cephean-dsp.h"

using namespace cephean;

//======================================================

//config definitions
#define PLUG_URI "https://github.com/creesullivan/cephean/octaver"
#define PLUG_CLASS octaver

#define PLUG_CONTROLS (6)
#define PLUG_INPUTS (1)
#define PLUG_OUTPUTS (1)

//======================================================

//class definition
class octaver : public plugin
{
public:
	octaver(float setFs);
	~octaver();

	void clear();
	void init() override;
	void step(int len) override;
	void deactivate() override;

private:
	//ports -----------------
	fcontrol downone{ this, 0 };
	fcontrol dry{ this, 1 };
	fcontrol uphalf{ this, 2 };
	fcontrol upone{ this, 3 };
	fcontrol uptwo{ this, 4 };
	fcontrol color{ this, 5 };

	input in{ this, 6 };

	output out{ this, 7 };

	//constants -------------
	const float slewdur = 0.1f;				//seconds

	const float thresh = 1e-8f;				//linear power thresh for correlation safety (-80 dB)

	const int bufflen = 2048;				//length of buffer including chunk at nominal rate
	const int chunklen = 448;				//length of chunk at nominal rate
	const int blendlen = 256;				//length of crossfading at nominal rate
	const int minper = 96;					//minimum period length at nominal rate
	const int DSR = 8;						//downsampling rate of correlation at nominal rate
	const int Rlen = (2048 - 448) / 8 + 1;	//length of correlation sequence

	const int maxlen = 128;					//maximum block size
	const float Tblend = 0.95f;				//blending allowance
	const float Tsnap = 0.5f;				//correlation threshold indicating a transient

	//state etc.-------------
	rebuffer rebuff;

	corrbuffer buff; //primary buffer with time locking for correlation, used for down octave
	downshift downoct; //drop octave control
	upshift upfifth; //up fifth control
	upshift upoct; //up octave control
	upshift up2oct; //up 2x octave control

	fvec R; //scratch storage for autocorrelation sequence

	gain gainDownOne; //gains for each voice
	gain gainDry;
	gain gainUpHalf;
	gain gainUpOne;
	gain gainUpTwo;

	fvec inscr, vscr; //scratch memory

	slewed<float> slewcolor; //slew for color filter design
	sofcasc<2> filtColor; //4th order lowpass on wet path for color

	//helper functions ------
	void update(int len = 0);
};
