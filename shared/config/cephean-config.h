
#pragma once

//--------------------------------------------------
// Defines the build environment
//--------------------------------------------------

#define ISWINDOWS (true)
#define ISARM (false)

#if ISWINDOWS
#include <mmintrin.h>
#endif

#if ISARM
#include <fenv.h>
#endif
