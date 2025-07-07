
#pragma once

//--------------------------------------------------
// Defines the module's audio processing behavior
//--------------------------------------------------

#include "cephean-lv2.h"

using namespace cephean;

//======================================================

//config definitions
#define PLUG_URI "https://github.com/creesullivan/cephean/gain"
#define PLUG_CLASS gain

//======================================================

//class definition
class gain : public plugin
{
public:
	gain(double setFs);
	~gain();

	void connect_port(uint32_t port, void* data) override;
	void init() override;
	void step(uint32_t len) override;
	void deactivate() override;

private:
	//controls---------------
	const float* level = nullptr;

	//inputs-----------------
	const float* in = nullptr;

	//outputs----------------
	float* out = nullptr;
};
