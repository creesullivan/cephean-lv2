
//--------------------------------------------------
// Defines the module's audio processing behavior
//--------------------------------------------------

#include "gain.h"

gain::gain(double setFs) : plugin(setFs) {}
gain::~gain() {}

void gain::connect_port(uint32_t port, void* data)
{
	switch (port)
	{
	case 0:
		connect_control(level, data);
		break;
	case 1:
		connect_input(in, data);
		break;
	case 2:
		connect_output(out, data);
		break;
	}
}

void gain::init() {}

void gain::step(uint32_t len)
{
	//fix stack
	const float G = *level;
	const float* const pin = in;
	float* const pout = out;

	//unslewed for now, super basic test
	const float g = (G > -100.0f) ? powf(10.0f, G / 20.0f) : 0.0f;

	for (uint32_t l = 0; l < len; l++) {
		pout[l] = pin[l] * g;
	}
}

void gain::deactivate() {}