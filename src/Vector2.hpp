#pragma once

#include "Core.h"
#include "Vector.hpp"

namespace Math
{
#define DECLARE_VEC2(Type, Precision)\
	template<>\
	struct MATHLIB Vector<2, Type, Precision>\
	{\
		VECTOR_BODY(2, Type, Precision)\
	\
	public:\
		union { Type x, r, s; };\
		union { Type y, g, t; };\
		\
		/* -- Constructors -- */\
		Vector<2, Type, Precision>(const Type x, const Type y);\
		\
		/* -- Getters --*/\
		Precision Cross(const Vec2 &vector) const;\
		Vec2 Rotated(const Precision angle) const;\
		\
		static inline Vec2 Zero()\
		{\
			return Vec2(0, 0);\
		}\
		\
		static inline Vec2 One()\
		{\
			return Vec2(1, 1);\
		}\
		\
		static inline Vec2 Left()\
		{\
			return Vec2(-1, 0);\
		}\
		\
		static inline Vec2 Right()\
		{\
			return Vec2(1, 0);\
		}\
		\
		static inline Vec2 Down()\
		{\
			return Vec2(0, -1);\
		}\
		\
		static inline Vec2 Up()\
		{\
			return Vec2(0, 1);\
		}\
		\
		/* -- Transformations -- */\
		Vec2 &Rotate(const Precision angle);\
	}

	DECLARE_VEC2(int, float);
	DECLARE_VEC2(float, float);
	DECLARE_VEC2(double, double);

#undef DECLARE_VEC2
}

