#pragma once

#include "Vector3.hpp"
#include <ostream>

namespace Math
{
	template<typename T, typename P>
	struct QuaternionT
	{
		T w;
		T x;
		T y;
		T z;

		// -- Constructors --

		QuaternionT() : w(T(0)), x(T(0)), y(T(0)), z(T(0)) {}

		QuaternionT(const T w, const T x, const T y, const T z)
			: w(w), x(x), y(y), z(z) {}

		QuaternionT(const Vector<3, T, P> &vec, const T scalar)
			: w(scalar), x(vec.x), y(vec.y), z(vec.z) {}

		QuaternionT(const QuaternionT &quat)
			: w(quat.w), x(quat.x), y(quat.y), z(quat.z) {}

		// -- Getters --

		inline QuaternionT GetConjugate() const
		{
			return QuaternionT(w, -x, -y, -z);
		}

		inline T GetScalar() const
		{
			return w;
		}

		inline Vector<3, T, P> GetVector() const
		{
			Vector<3, T, P> res;
			res[0] = x;
			res[1] = y;
			res[2] = z;
			return res;
		}

		inline P Magnitude() const
		{
			return std::sqrt(w * w + x * x + y * y + z * z);
		}

		inline P MagnitudeSquared() const
		{
			return w * w + x * x + y * y + z * z;
		}

		QuaternionT Inversed() const
		{
			const T sqrlength = MagnitudeSquared();
			return Quaternion(w / sqrlength, -x / sqrlength, -y / sqrlength, -z / sqrlength);
		}

		Vector<3, T, P> ToEulerAngles() const
		{
			const T pole = x * y + z * w;
			Vector<3, T, P> res;

			if (pole == T(0.5))
			{
				res[0] = T(2) * std::atan2(x, w);
				res[1] = std::asin(T(2) * (w * y - x * z));
				res[2] = T(0);
				return res;
			}

			if (pole == T(-0.5))
			{
				res[0] = T(-2) * std::atan2(x, w);
				res[1] = std::asin(T(2) * (w * y - x * z));
				res[2] = T(0);
				return res;
			}

			res[0] = std::atan2(T(2) * (w * x + y * z), T(1) - T(2) * (x * x + y * y));
			res[1] = std::asin(T(2) * (w * y - x * z));
			res[2] = std::atan2(T(2) * (w * z + x * y), T(1) - T(2) * (y * y + z * z));
			return res;
		}

		// -- Transformations --

		inline QuaternionT &Conjugate()
		{
			x = -x;
			y = -y;
			z = -z;
			return *this;
		}

		QuaternionT &Inverse()
		{
			const T sqrlength = MagnitudeSquared();
			w /= sqrlength;
			x = -x / sqrlength;
			y = -y / sqrlength;
			z = -z / sqrlength;
			return *this;
		}

		// -- Unary arithmetic operators --

		inline QuaternionT &operator+=(const QuaternionT &quat)
		{
			w += quat.w;
			x += quat.x;
			y += quat.y;
			z += quat.z;
			return *this;
		}

		inline QuaternionT &operator-=(const QuaternionT &quat)
		{
			w -= quat.w;
			x -= quat.x;
			y -= quat.y;
			z -= quat.z;
			return *this;
		}

		QuaternionT &operator*=(const QuaternionT &quat)
		{
			const T temp_x = x;
			const T temp_y = y;
			const T temp_z = z;
			const T temp_w = w;

			w = temp_w * quat.w - temp_x * quat.x - temp_y * quat.y - temp_z * quat.z;
			x = temp_w * quat.x + temp_x * quat.w + temp_y * quat.z - temp_z * quat.y;
			y = temp_w * quat.y - temp_x * quat.z + temp_y * quat.w + temp_z * quat.x;
			z = temp_w * quat.z + temp_x * quat.y - temp_y * quat.x + temp_z * quat.w;
			return *this;
		}

		inline QuaternionT &operator*=(const T val)
		{
			w *= val;
			x *= val;
			y *= val;
			z *= val;
			return *this;
		}

		inline QuaternionT &operator/=(const QuaternionT &quat)
		{
			return *this *= quat.Inversed();
		}

		inline QuaternionT &operator/=(const T val)
		{
			w /= val;
			x /= val;
			y /= val;
			z /= val;
			return *this;
		}

		// -- Unary operators --

		inline QuaternionT operator+() const
		{
			return *this;
		}

		inline QuaternionT operator-() const
		{
			return QuaternionT(-w, -x, -y, -z);
		}

		// -- Binary operators --

		inline QuaternionT operator+(const QuaternionT &quat) const
		{
			return QuaternionT(
				w + quat.w,
				x + quat.x,
				y + quat.y,
				z + quat.z
			);
		}

		inline QuaternionT operator-(const QuaternionT &quat) const
		{
			return QuaternionT(
				w - quat.w,
				x - quat.x,
				y - quat.y,
				z - quat.z
			);
		}

		inline QuaternionT operator*(const QuaternionT &quat) const
		{
			return QuaternionT(
				w * quat.w - x * quat.x - y * quat.y - z * quat.z,
				w * quat.x + x * quat.w + y * quat.z - z * quat.y,
				w * quat.y - x * quat.z + y * quat.w + z * quat.x,
				w * quat.z + x * quat.y - y * quat.x + z * quat.w
			);
		}

		Vector<3, T, P> operator*(const Vector<3, T, P> &vec) const
		{
			return Vector<3, T, P>(
				(T(1) - T(2) * y * y - T(2) * z * z) * vec[0] + T(2) * (x * y + w * z) * vec[1] + T(2) * (x * z - w * y) * vec[2],
				T(2) * (x * y - w * z) * vec[0] + (T(1) - T(2) * x * x - T(2) * z * z) * vec[1] + T(2) * (x * z + w * y) * vec[2],
				T(2) * (x * z + w * y) * vec[0] + T(2) * (y * z - w * x) * vec[1] + (T(1) - T(2) * x * x - T(2) * y * y) * vec[2]
			);
		}

		inline QuaternionT operator*(const T val) const
		{
			return QuaternionT(
				w * val,
				x * val,
				y * val,
				z * val
			);
		}

		inline QuaternionT operator/(const QuaternionT &quat) const
		{
			return *this * quat.Inversed();
		}

		inline QuaternionT operator/(const T val) const
		{
			return QuaternionT(
				w / val,
				x / val,
				y / val,
				z / val
			);
		}

		// -- Boolean operators --

		inline bool operator!=(const QuaternionT &quat) const
		{
			return w != quat.w || x != quat.x || y != quat.y || z != quat.z;
		}

		inline bool operator==(const QuaternionT &quat) const
		{
			return !(*this != quat);
		}

		// -- Convertion operator --

		template<typename U, typename Q>
		operator QuaternionT<U, Q>() const
		{
			return QuaternionT<U, Q>((U)w, (U)x, (U)y, (U)z);
		}

		// -- Stream operators --

		friend std::ostream &operator<<(std::ostream &ostream, const QuaternionT &quat)
		{
			ostream << "(w: " << quat.w << ", x: " << quat.x << ", y: " << quat.y << ", z: " << quat.z << ")";
		}

		// -- Static getters --

		static inline QuaternionT Identity()
		{
			return QuaternionT(T(0), T(0), T(0), T(1));
		}

		// -- Static methods --

		static QuaternionT FromAngleAxis(const Vector<3, T, P> &axis, const T angle)
		{
			if (!axis.LengthSquared())
			{
				return Identity();
			}

			if (!axis.IsNormalized())
			{
				axis.Normalize();
			}

			T sin = std::sin(angle * T(0.5));
			return QuaternionT(std::cos(angle * T(0.5)), axis.x * sin, axis.y * sin, axis.z * sin);
		}

		static QuaternionT FromEulerAngles(const T roll, const T pitch, const T yaw)
		{
			T cr = std::cos(roll * T(0.5));
			T sr = std::sin(roll * T(0.5));

			T cp = std::cos(pitch * T(0.5));
			T sp = std::sin(pitch * T(0.5));

			T cy = std::cos(yaw * T(0.5));
			T sy = std::sin(yaw * T(0.5));

			return QuaternionT(
				sr * cp * cy - cr * sp * sy,
				cr * sp * cy + sr * cp * sy,
				cr * cp * sy - sr * sp * cy,
				cr * cp * cy + sr * sp * sy
			);
		}
	};

	typedef QuaternionT<float, float> Quaternion;
	typedef QuaternionT<double, double> QuaternionD;
}

