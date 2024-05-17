#pragma once

#include "Vector.hpp"

namespace Math
{
	enum class EulerAngleOrder
	{
		XYZ,
		XZY,
		YXZ,
		YZX,
		ZXY,
		ZYX,
	};

#define DECLARE_VEC3(Type, Precision)\
	template<>\
	struct Vector<3, Type, Precision>\
	{\
		VECTOR_BODY(3, Type, Precision)\
	\
	private:\
		Vec3 &RotateAroundX(const Precision cosAngle, const Precision sinAngle);\
		Vec3 &RotateAroundY(const Precision cosAngle, const Precision sinAngle);\
		Vec3 &RotateAroundZ(const Precision cosAngle, const Precision sinAngle);\
	\
	public:\
		union { Type x, r, s; };\
		union { Type y, g, t; };\
		union { Type z, b, p; };\
		\
		/* -- Constructors -- */\
		Vector<3, Type, Precision>(const Type x, const Type y, const Type z);\
		\
		/* -- Getters -- */\
		Vec3 Cross(const Vec3 &vector) const;\
		Vec3 Rotated(const Vector<3, Precision, Precision> &rotation, const EulerAngleOrder order = EulerAngleOrder::XYZ) const;\
		\
		static inline Vec3 Zero()\
		{\
			return Vec3(0, 0, 0);\
		}\
		\
		static inline Vec3 One()\
		{\
			return Vec3(1, 1, 1);\
		}\
		\
		static inline Vec3 Left()\
		{\
			return Vec3(-1, 0, 0);\
		}\
		\
		static inline Vec3 Right()\
		{\
			return Vec3(1, 0, 0);\
		}\
		\
		static inline Vec3 Down()\
		{\
			return Vec3(0, -1, 0);\
		}\
		\
		static inline Vec3 Up()\
		{\
			return Vec3(0, 1, 0);\
		}\
		\
		static inline Vec3 Forward()\
		{\
			return Vec3(0, 0, 1);\
		}\
		\
		static inline Vec3 Backward()\
		{\
			return Vec3(0, 0, -1);\
		}\
		\
		/* -- Transformations -- */\
		Vec3 &RotateAroundX(const Precision angle);\
		Vec3 &RotateAroundY(const Precision angle);\
		Vec3 &RotateAroundZ(const Precision angle);\
		Vec3 &Rotate(const Vector<3, Precision, Precision> &rotation, const EulerAngleOrder order = EulerAngleOrder::XYZ);\
	}

	DECLARE_VEC3(int, float);
	DECLARE_VEC3(float, float);
	DECLARE_VEC3(double, double);
}

