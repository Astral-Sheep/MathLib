#pragma once

#include "Vector.hpp"
#include <ostream>

namespace Math
{
	template<typename T, typename P>
	struct QuaternionT
	{
	private:
		typedef QuaternionT<T, P> Quat;
	public:
		T w;
		T x;
		T y;
		T z;

		// -- Constructors --

		QuaternionT<T, P>() : w(T(0)), x(T(0)), y(T(0)), z(T(0)) {}

		QuaternionT<T, P>(const T pW, const T pX, const T pY, const T pZ)
			: w(pW), x(pX), y(pY), z(pZ) {}

		QuaternionT<T, P>(const Vector<3, T, P> &pVec, const T pScalar)
			: w(pScalar), x(pVec.x), y(pVec.y), z(pVec.z) {}

		QuaternionT<T, P>(const Quat &pQuat)
			: w(pQuat.w), x(pQuat.x), y(pQuat.y), z(pQuat.z) {}

		// -- Getters --

		inline Quat GetConjugate() const
		{
			return Quat(w, -x, -y, -z);
		}

		inline T GetScalar() const
		{
			return w;
		}

		inline Vector<3, T, P> GetVector() const
		{
			Vector<3, T, P> lRes;
			lRes[0] = x;
			lRes[1] = y;
			lRes[2] = z;
			return lRes;
		}

		inline P Magnitude() const
		{
			return std::sqrt(w * w + x * x + y * y + z * z);
		}

		inline P MagnitudeSquared() const
		{
			return w * w + x * x + y * y + z * z;
		}

		Quat Inversed() const
		{
			const T lSqrlength = MagnitudeSquared();
			return Quat(w / lSqrlength, -x / lSqrlength, -y / lSqrlength, -z / lSqrlength);
		}

		Vector<3, T, P> ToEulerAngles() const
		{
			const T lPole = x * y + z * w;
			Vector<3, T, P> lRes;

			if (lPole == T(0.5))
			{
				lRes[0] = T(2) * std::atan2(x, w);
				lRes[1] = std::asin(T(2) * (w * y - x * z));
				lRes[2] = T(0);
				return lRes;
			}

			if (lPole == T(-0.5))
			{
				lRes[0] = T(-2) * std::atan2(x, w);
				lRes[1] = std::asin(T(2) * (w * y - x * z));
				lRes[2] = T(0);
				return lRes;
			}

			lRes[0] = std::atan2(T(2) * (w * x + y * z), T(1) - T(2) * (x * x + y * y));
			lRes[1] = std::asin(T(2) * (w * y - x * z));
			lRes[2] = std::atan2(T(2) * (w * z + x * y), T(1) - T(2) * (y * y + z * z));
			return lRes;
		}

		// -- Transformations --

		inline Quat &Conjugate()
		{
			x = -x;
			y = -y;
			z = -z;
			return *this;
		}

		Quat &Inverse()
		{
			const T lSqrlength = MagnitudeSquared();
			w /= lSqrlength;
			x = -x / lSqrlength;
			y = -y / lSqrlength;
			z = -z / lSqrlength;
			return *this;
		}

		// -- Unary arithmetic operators --

		inline Quat &operator+=(const Quat &pQuat)
		{
			w += pQuat.w;
			x += pQuat.x;
			y += pQuat.y;
			z += pQuat.z;
			return *this;
		}

		inline Quat &operator-=(const Quat &pQuat)
		{
			w -= pQuat.w;
			x -= pQuat.x;
			y -= pQuat.y;
			z -= pQuat.z;
			return *this;
		}

		Quat &operator*=(const Quat &pQuat)
		{
			const T lTempX = x;
			const T lTempY = y;
			const T lTempZ = z;
			const T lTempW = w;

			w = lTempW * pQuat.w - lTempX * pQuat.x - lTempY * pQuat.y - lTempZ * pQuat.z;
			x = lTempW * pQuat.x + lTempX * pQuat.w + lTempY * pQuat.z - lTempZ * pQuat.y;
			y = lTempW * pQuat.y - lTempX * pQuat.z + lTempY * pQuat.w + lTempZ * pQuat.x;
			z = lTempW * pQuat.z + lTempX * pQuat.y - lTempY * pQuat.x + lTempZ * pQuat.w;
			return *this;
		}

		inline Quat &operator*=(const T pVal)
		{
			w *= pVal;
			x *= pVal;
			y *= pVal;
			z *= pVal;
			return *this;
		}

		inline Quat &operator/=(const Quat &pQuat)
		{
			return *this *= pQuat.Inversed();
		}

		inline Quat &operator/=(const T pVal)
		{
			w /= pVal;
			x /= pVal;
			y /= pVal;
			z /= pVal;
			return *this;
		}

		// -- Unary operators --

		inline Quat operator+() const
		{
			return *this;
		}

		inline Quat operator-() const
		{
			return Quat(-w, -x, -y, -z);
		}

		// -- Binary operators --

		inline Quat operator+(const Quat &pQuat) const
		{
			return Quat(
				w + pQuat.w,
				x + pQuat.x,
				y + pQuat.y,
				z + pQuat.z
			);
		}

		inline Quat operator-(const Quat &pQuat) const
		{
			return Quat(
				w - pQuat.w,
				x - pQuat.x,
				y - pQuat.y,
				z - pQuat.z
			);
		}

		inline Quat operator*(const Quat &pQuat) const
		{
			return Quat(
				w * pQuat.w - x * pQuat.x - y * pQuat.y - z * pQuat.z,
				w * pQuat.x + x * pQuat.w + y * pQuat.z - z * pQuat.y,
				w * pQuat.y - x * pQuat.z + y * pQuat.w + z * pQuat.x,
				w * pQuat.z + x * pQuat.y - y * pQuat.x + z * pQuat.w
			);
		}

		Vector<3, T, P> operator*(const Vector<3, T, P> &pVec) const
		{
			return Vector<3, T, P>(
				(T(1) - T(2) * y * y - T(2) * z * z) * pVec[0] + T(2) * (x * y + w * z) * pVec[1] + T(2) * (x * z - w * y) * pVec[2],
				T(2) * (x * y - w * z) * pVec[0] + (T(1) - T(2) * x * x - T(2) * z * z) * pVec[1] + T(2) * (x * z + w * y) * pVec[2],
				T(2) * (x * z + w * y) * pVec[0] + T(2) * (y * z - w * x) * pVec[1] + (T(1) - T(2) * x * x - T(2) * y * y) * pVec[2]
			);
		}

		inline Quat operator*(const T pVal) const
		{
			return Quat(
				w * pVal,
				x * pVal,
				y * pVal,
				z * pVal
			);
		}

		inline Quat operator/(const Quat &pQuat) const
		{
			return *this * pQuat.Inversed();
		}

		inline Quat operator/(const T pVal) const
		{
			return Quat(
				w / pVal,
				x / pVal,
				y / pVal,
				z / pVal
			);
		}

		// -- Boolean operators --

		inline bool operator!=(const Quat &pQuat) const
		{
			return w != pQuat.w || x != pQuat.x || y != pQuat.y || z != pQuat.z;
		}

		inline bool operator==(const Quat &pQuat) const
		{
			return !(*this != pQuat);
		}

		// -- Convertion operator --

		template<typename U, typename Q>
		operator QuaternionT<U, Q>() const
		{
			return QuaternionT<U, Q>((U)w, (U)x, (U)y, (U)z);
		}

		// -- Stream operators --

		friend std::ostream &operator<<(std::ostream &pOStream, const Quat &pQuat)
		{
			pOStream << "(w: " << pQuat.w << ", x: " << pQuat.x << ", y: " << pQuat.y << ", z: " << pQuat.z << ")";
		}

		// -- Static getters --

		static inline Quat Identity()
		{
			return Quat(T(0), T(0), T(0), T(1));
		}

		// -- Static methods --

		static Quat FromAngleAxis(const Vector<3, T, P> &pAxis, const T pAngle)
		{
			if (!pAxis.LengthSquared())
			{
				return Identity();
			}

			if (!pAxis.IsNormalized())
			{
				pAxis.Normalize();
			}

			const T lSin = std::sin(pAngle * T(0.5));
			return Quat(std::cos(pAngle * T(0.5)), pAxis.x * lSin, pAxis.y * lSin, pAxis.z * lSin);
		}

		static Quat FromEulerAngles(const T pRoll, const T pRitch, const T pYaw)
		{
			const T lCR = std::cos(pRoll * T(0.5));
			const T lSR = std::sin(pRoll * T(0.5));

			const T lCP = std::cos(pRitch * T(0.5));
			const T lSP = std::sin(pRitch * T(0.5));

			const T lCY = std::cos(pYaw * T(0.5));
			const T lSY = std::sin(pYaw * T(0.5));

			return Quat(
				lSR * lCP * lCY - lCR * lSP * lSY,
				lCR * lSP * lCY + lSR * lCP * lSY,
				lCR * lCP * lSY - lSR * lSP * lCY,
				lCR * lCP * lCY + lSR * lSP * lSY
			);
		}
	};

	typedef QuaternionT<float, float> Quaternion;
	typedef QuaternionT<double, double> QuaternionD;
}

