
#pragma once

//--------------------------------------------------
// Defines the module's audio processing behavior
//--------------------------------------------------

#include "cephean-lv2.h"
#include "cephean-dsp.h"

using namespace cephean;

//======================================================

//config definitions
#define PLUG_URI "https://github.com/creesullivan/cephean/plugintemplate"
#define PLUG_CLASS plugintemplate

#define PLUG_CONTROLS (1)
#define PLUG_INPUTS (1)
#define PLUG_OUTPUTS (1)

//======================================================

//class definition
class plugintemplate : public plugin
{
public:
	plugintemplate(float setFs);
	~plugintemplate();

	void init() override;
	void step(int len) override;
	void deactivate() override;

private:
	//ports -----------------
	fcontrol level{ this, 0 };

	input in{ this, 1 };

	output out{ this, 2 };

	//constants -------------
	const float slewdur = 0.1f;			//seconds

	const int maxlen = 256;				//maximum block size

	//state etc.-------------
	rebuffer rebuff;

	gain slewedGain;

	//helper functions ------
	void update(int len = 0);
};
