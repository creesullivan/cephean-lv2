
//--------------------------------------------------
// Compiles the DSP .dll of the LV2 plugin
//--------------------------------------------------

#include "fbdistortion.h"

typedef PLUG_CLASS plug;

static LV2_Handle instantiate(const LV2_Descriptor* descriptor, double rate,
	const char* bundle_path, const LV2_Feature* const* features)
{
	plug* H = new plug((float)rate);
	return (LV2_Handle)H;
}

static void connect_port(LV2_Handle instance, uint32_t port, void* data)
{
	plug* H = (plug*)instance;
	H->connect_port(port, data);
}

static void activate(LV2_Handle instance)
{
	plug* H = (plug*)instance;
	H->init();
}

static void run(LV2_Handle instance, uint32_t n_samples)
{
	plug* H = (plug*)instance;
	H->step((int)n_samples);
}

static void deactivate(LV2_Handle instance)
{
	plug* H = (plug*)instance;
	H->deactivate();
}

static void cleanup(LV2_Handle instance)
{
	plug* H = (plug*)instance;
	delete H;
}

static const void* extension_data(const char* uri)
{
	return NULL; //unused
}

static const LV2_Descriptor descriptor = { PLUG_URI,
										  instantiate,
										  connect_port,
										  activate,
										  run,
										  deactivate,
										  cleanup,
										  extension_data };

extern "C" LV2_SYMBOL_EXPORT const LV2_Descriptor* lv2_descriptor(uint32_t index)
{
	return index == 0 ? &descriptor : NULL;
}