#pragma once

#include "Core.h"
#include "Vector.hpp"

namespace Math
{
#define DECLARE_VEC4(Type, Precision)\
	template<>\
	struct MATHLIB Vector<4, Type, Precision>\
	{\
		VECTOR_BODY(4, Type, Precision)\
	\
	public:\
		union { Type x, r, s; };\
		union { Type y, g, t; };\
		union { Type z, b, p; };\
		union { Type w, a, q; };\
		\
		/* -- Constructors -- */\
		Vector<4, Type, Precision>(const Type x, const Type y, const Type z, const Type w);\
		\
		/* -- Getters -- */\
		static inline Vec4 Zero()\
		{\
			return Vec4(0, 0, 0, 0);\
		}\
		\
		static inline Vec4 One()\
		{\
			return Vec4(1, 1, 1, 1);\
		}\
	}

	DECLARE_VEC4(int, float);
	DECLARE_VEC4(float, float);
	DECLARE_VEC4(double, double);

#undef DECLARE_VEC4
}

