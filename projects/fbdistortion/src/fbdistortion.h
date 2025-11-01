
#pragma once

//--------------------------------------------------
// Defines the module's audio processing behavior
//--------------------------------------------------

#include "cephean-lv2.h"
#include "cephean-dsp.h"

using namespace cephean;

//======================================================

//config definitions
#define PLUG_URI "https://github.com/creesullivan/cephean/fbdistortion"
#define PLUG_CLASS fbdistortion

#define PLUG_CONTROLS (6)
#define PLUG_INPUTS (2)
#define PLUG_OUTPUTS (1)

//======================================================

//class definition
class fbdistortion : public plugin
{
public:
	fbdistortion(float setFs);
	~fbdistortion();

	void init() override;
	void step(int len) override;
	void deactivate() override;

private:
	//ports -----------------
	fcontrol drive{ this, 0 };
	fcontrol shape{ this, 1 };
	fcontrol feedback{ this, 2 };
	fcontrol dirty{ this, 3 };
	fcontrol color{ this, 4 };
	fcontrol level{ this, 5 };

	input in{ this, 6 };
	input sc{ this, 7 };

	output out{ this, 8 };

	//constants -------------
	const float slewdur = 0.1f;				//seconds duration for parameter slews
	const float stringdur1 = 0.1f;			//seconds string state fade duration during sustain
	const float stringdur2 = 0.02f;			//seconds string state fade duration during transients/mute

	const int maxlen = 128;					//maximum block size

	const float thresh = 1e-6f;				//linear power thresh for correlation safety (-60 dB)
	const int bufflen = 2048;				//length of correlation buffer including chunk at nominal rate
	const int chunklen = 800;				//length of chunk at nominal rate
	const int maxper = 2048 - 800;			//maximum period detect length at nominal rate
	const int minper = 129;					//minimum period detect length at nominal rate
	const int DSR = 8;						//downsampling rate of correlation at nominal rate
	const int Rlen = (2048 - 800) / 8 + 1;	//length of downsampled correlation sequence

	const float Rmax = 0.9f;				//R threshold to fully open the string model
	const float Rmin = 0.5f;				//R threshold to mute the string model
	const float permult = 1.0f;			//period multiplier for string based on detect period (careful not to overrun!)

	const float fstring = 5000.0f;			//string model lowpass frequency in Hz
	const float tstring = 1.0f;				//string model approximate build time in sec
	const float fDC = 400.0f;				//string model DC rejection post-filter in Hz

	const float decay0 = 0.0f;				//decay periods when string model is fully open
	const float sbeta0 = 1.0f;				//fixed string model reflection numerator coef, set on construct (1 - expf(-1.0f / decay0))

	const float attdur = 0.002f;			//attack envelope time in seconds
	const float holddur = 0.02f;			//hold envelope time in seconds
	const float reldur = 0.02f;				//release envelope time in seconds

	//state etc.-------------
	rebuffer rebuff;
	fvec inscr;
	fvec scscr;

	slewed<float> slewColor;

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
	faded<stringParam> stateString; //fractional period length
	circbuffer buffString; //string model core
	fof filtString; //fixed string model feedback lowpass filter
	fof DCreject; //fixed string model DC rejection filter
	gain fbGain; //string model feedback drive gain

	gain driveGain; //drive gain into the nonlinearity

	attrel envAttack; //envelope detection for input compressor
	hold4 envHold;
	attrel envRelease;
	
	fastratio compress; //compression curve
	simshaper distort; //distortion curve with lowpass filtering
	gain makeupGain; //makeup gain after the nonlinearity
	sof postFilter; //high freq makeup boost

	sof filtFeedback; //feedback color resonant LPF
	gain outGain; //final output level gain

	fvec R; //scratch storage for correlation vector
	fvec stscr, scr1, scr2, scr3; //general purpose scratch vectors

	//helper functions ------
	void update(int len = 0);
	void clear();
};
