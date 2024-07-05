#pragma once

#ifdef MATH_EXPORT
	#define MATHLIB __declspec(dllexport)
#elif defined MATH_IMPORT
	#define MATHLIB __declspec(dllimport)
#else
	#define MATHLIB
#endif

#ifdef HIGH_PRECISION
typedef double float_type;
#else
typedef float float_type;
#endif

