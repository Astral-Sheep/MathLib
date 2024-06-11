#pragma once

#ifdef MATH_EXPORT
	#define MATHLIB __declspec(dllexport)
#elif defined MATH_IMPORT
	#define MATHLIB __declspec(dllimport)
#else
	#define MATHLIB
#endif

