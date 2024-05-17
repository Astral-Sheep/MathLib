#pragma once

#include "Vector3.hpp"
#include <ostream>

namespace Math
{
	struct QuaternionD;

	struct Quaternion
	{
		float x;
		float y;
		float z;
		float w;

		// -- Constructors --
		Quaternion();
		Quaternion(const float x, const float y, const float z, const float w);
		Quaternion(const Vector3 &vector, const float scalar);
		Quaternion(const Quaternion &quaternion);

		// -- Unary arithmetic operators --
		Quaternion &operator+=(const Quaternion &quaternion);
		Quaternion &operator-=(const Quaternion &quaternion);
		Quaternion &operator*=(const Quaternion &quaternion);
		Quaternion &operator*=(const float scalar);
		Quaternion &operator/=(const Quaternion &quaternion);
		Quaternion &operator/=(const float scalar);

		// -- Unary operators --
		Quaternion operator-() const;

		// -- Binary operators --
		Quaternion operator+(const Quaternion &quaternion) const;
		Quaternion operator-(const Quaternion &quaternion) const;
		Quaternion operator*(const Quaternion &quaternion) const;
		Vector3 operator*(const Vector3 &vector) const;
		Quaternion operator*(const float scalar) const;
		Quaternion operator/(const Quaternion &quaternion) const;
		Quaternion operator/(const float scalar) const;

		// -- Boolean operators --
		bool operator==(const Quaternion &quaternion) const;
		bool operator!=(const Quaternion &quaternion) const;

		// -- Convertion operators --
		operator QuaternionD() const;

		// -- Stream operators --
		friend std::ostream &operator<<(std::ostream &ostream, const Quaternion &quaternion);

		// -- Static getters --
		static inline Quaternion Identity()
		{
			return Quaternion(0.0f, 0.0f, 0.0f, 1.0f);
		}

		// -- Static methods --
		static Quaternion AngleAxis(Vector3 axis, const float angle);

		// -- Getters --
		Quaternion Inversed() const;
		float Magnitude() const;
		float MagnitudeSquared() const;

		inline Quaternion GetConjugate() const
		{
			return Quaternion(-x, -y, -z, w);
		}

		inline float Scalar() const
		{
			return w;
		}

		inline Vector3 Vector() const
		{
			return Vector3(x, y, z);
		}

		// -- Transformations --
		Quaternion &Conjugate();
		Quaternion &Inverse();
	};
}

