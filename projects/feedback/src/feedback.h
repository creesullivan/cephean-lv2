
#pragma once

//--------------------------------------------------
// Defines the module's audio processing behavior
//--------------------------------------------------

#include "cephean-lv2.h"
#include "cephean-dsp.h"

using namespace cephean;

//======================================================

//config definitions
#define PLUG_URI "https://github.com/creesullivan/cephean/feedback"
#define PLUG_CLASS feedback

#define PLUG_CONTROLS (2)
#define PLUG_INPUTS (3)
#define PLUG_OUTPUTS (1)

//======================================================

//class definition
class feedback : public plugin
{
public:
	feedback(float setFs);
	~feedback();

	void init() override;
	void step(int len) override;
	void deactivate() override;

private:
	//ports -----------------
	fcontrol amount{ this, 0 };
	fcontrol dry{ this, 1 };

	input in{ this, 2 };
	input sc{ this, 3 };
	input fb{ this, 4 };

	output out{ this, 5 };

	//constants -------------
	const float slewdur = 0.1f;				//seconds duration for parameter slews
	const float stringdur1 = 0.1f;			//seconds string state fade duration during sustain
	const float stringdur2 = 0.02f;			//seconds string state fade duration during transients/mute

	const int maxlen = 128;					//maximum permitted device block size for good operation

	const float thresh = 1e-6f;				//linear power thresh for correlation safety (-60 dB)
	const int bufflen = 2048;				//length of correlation buffer including chunk at nominal rate
	const int chunklen = 800;				//length of chunk at nominal rate
	const int maxper = 2048 - 800;			//maximum period detect length at nominal rate
	const int minper = 129;					//minimum period detect length at nominal rate
	const int DSR = 8;						//downsampling rate of correlation at nominal rate
	const int Rlen = (2048 - 800) / 8 + 1;	//length of downsampled correlation sequence

	const float Rmax = 0.9f;				//R threshold to fully open the string model
	const float Rmin = 0.5f;				//R threshold to mute the string model

	const float fstring = 5000.0f;			//string model lowpass frequency in Hz
	const float tstring = 1.0f;				//string model approximate build time in sec
	const float fDC = 400.0f;				//string model DC rejection post-filter in Hz

	const float decay0 = 0.0f;				//decay periods when string model is fully open
	const float sbeta0 = 1.0f;				//fixed string model reflection numerator coef, set on construct (1 - expf(-1.0f / decay0))

	//state etc.-------------
	rebuffer rebuff;
	fvec inscr;
	fvec scscr;
	fvec fbscr;

	corrbuffer detect; //fundamental detection buffer

	//string state struct is faded every block based on detected properties
	struct stringParam
	{
		stringParam(float sustain = 0.0f, float setPeriod = 0.0f, float setConfidence = 0.0f);
		~stringParam();

		float sustain = 0.0f;
		float period = 0.0f;
		float confidence = 0.0f;

		bool operator ==(stringParam other) const;
		bool operator !=(stringParam other) const;
	};
	faded<stringParam> stateString; //string model state
	circbuffer buffString; //string model core
	flatbuffer stmem; //string output block for delay alignment
	fof filtString; //fixed string model feedback lowpass filter
	fof DCreject; //fixed string model DC rejection filter
	sof filtFeedback; //fixed feedback path lowpass filter
	gain fbGain; //string model feedback drive gain

	gain dryGain; //input signal dry mix gain

	fvec R; //scratch storage for correlation vector
	fvec stscr, scr1, scr2, scr3; //general purpose scratch vectors

	//helper functions ------
	void update(int len = 0);
	void clear();
};
