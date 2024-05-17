#pragma once

#include "Vector3.hpp"
#include <ostream>

namespace Math
{
	struct Quaternion;

	struct QuaternionD
	{
		double x;
		double y;
		double z;
		double w;

		// -- Constructors --
		QuaternionD();
		QuaternionD(const double x, const double y, const double z, const double w);
		QuaternionD(const Vector3D &vector, const double scalar);
		QuaternionD(const QuaternionD &quaternion);

		// -- Unary arithmetic operators --
		QuaternionD &operator+=(const QuaternionD &quaternion);
		QuaternionD &operator-=(const QuaternionD &quaternion);
		QuaternionD &operator*=(const QuaternionD &quaternion);
		QuaternionD &operator*=(const double scalar);
		QuaternionD &operator/=(const QuaternionD &quaternion);
		QuaternionD &operator/=(const double scalar);

		// -- Unary operators --
		QuaternionD operator-() const;

		// -- Binary operators --
		QuaternionD operator+(const QuaternionD &quaternion) const;
		QuaternionD operator-(const QuaternionD &quaternion) const;
		QuaternionD operator*(const QuaternionD &quaternion) const;
		Vector3D operator*(const Vector3D &vector) const;
		QuaternionD operator*(const double scalar) const;
		QuaternionD operator/(const QuaternionD &quaternion) const;
		QuaternionD operator/(const double scalar) const;

		// -- Boolean operators --
		bool operator==(const QuaternionD &quaternion) const;
		bool operator!=(const QuaternionD &quaternion) const;

		// -- Convertion operators --
		explicit operator Quaternion() const;

		// -- Stream operators --
		friend std::ostream &operator<<(std::ostream &ostream, const QuaternionD &quaternion);

		// -- Static getters --
		static inline QuaternionD Identity()
		{
			return QuaternionD(0.0, 0.0, 0.0, 1.0);
		}

		// -- Static methods --
		static QuaternionD AngleAxis(Vector3D axis, const double angle);

		// -- Getters --
		QuaternionD Inversed() const;
		double Magnitude() const;
		double MagnitudeSquared() const;

		inline QuaternionD GetConjugate() const
		{
			return QuaternionD(-x, -y, -z, w);
		}

		inline double Scalar() const
		{
			return w;
		}

		inline Vector3D Vector() const
		{
			return Vector3D(x, y, z);
		}

		// -- Transformations --
		QuaternionD &Conjugate();
		QuaternionD &Inverse();
	};
}

