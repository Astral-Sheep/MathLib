#pragma once

#include "Math.hpp"
#include "MatrixTransform.hpp"
#include "Vector.hpp"

namespace Math
{
	template<typename T, typename P>
	struct Vector<3, T, P>
	{
	private:
		typedef Vector<3, T, P> Vec;

		Vec &RotateAroundX(const P pCos, const P pSin)
		{
			const T lTempY = y;
			const T lTempZ = z;
			y = pCos * lTempY - pSin * lTempZ;
			z = pSin * lTempY + pCos * lTempZ;
			return *this;
		}

		Vec &RotateAroundY(const P pCos, const P pSin)
		{
			T lTempX = x;
			T lTempZ = z;
			x = pCos * lTempX + pSin * lTempZ;
			z = -pSin * lTempX + pCos * lTempZ;
			return *this;
		}

		Vec &RotateAroundZ(const P pCos, const P pSin)
		{
			T lTempX = x;
			T lTempY = y;
			x = pCos * lTempX - pSin * lTempY;
			y = pSin * lTempX + pCos * lTempY;
			return *this;
		}

	public:
		union { T x, r, s; };
		union { T y, g, t; };
		union { T z, b, p; };

		// -- Constructors --

		Vector<3, T, P>()
			: x(T(0)), y(T(0)), z(T(0)) {}

		explicit Vector<3, T, P>(const T pVal)
			: x(pVal), y(pVal), z(pVal) {}

		Vector<3, T, P>(const T pX, const T pY, const T pZ)
			: x(pX), y(pY), z(pZ) {}

		Vector<3, T, P>(const T pVals[3])
			: x(pVals[0]), y(pVals[1]), z(pVals[2]) {}

		Vector<3, T, P>(const Vec &pVec)
			: x(pVec.x), y(pVec.y), z(pVec.z) {}

		// -- Accesses --

		inline const T &operator[](const int pIndex) const noexcept
		{
			return *(&x + pIndex);
		}

		inline T &operator[](const int pIndex) noexcept
		{
			return *(&x + pIndex);
		}

		// -- Unary arithmetic operators --

		inline Vec &operator+=(const Vec &pVec)
		{
			x += pVec.x;
			y += pVec.y;
			z += pVec.z;
			return *this;
		}

		inline Vec &operator-=(const Vec &pVec)
		{
			x -= pVec.x;
			y -= pVec.y;
			z -= pVec.z;
			return *this;
		}

		inline Vec &operator*=(const Vec &pVec)
		{
			x *= pVec.x;
			y *= pVec.y;
			z *= pVec.z;
			return *this;
		}

		inline Vec &operator*=(const T pVal)
		{
			x *= pVal;
			y *= pVal;
			z *= pVal;
			return *this;
		}

		inline Vec &operator/=(const Vec &pVec)
		{
			x /= pVec.x;
			y /= pVec.y;
			z /= pVec.z;
			return *this;
		}

		inline Vec &operator/=(const T pVal)
		{
			x /= pVal;
			y /= pVal;
			z /= pVal;
			return *this;
		}

		// -- Unary operators --

		inline Vec operator+() const
		{
			return *this;
		}

		inline Vec operator-() const
		{
			return Vector<3, T, P>(-x, -y, -z);
		}

		// -- Binary operators --

		inline Vec operator+(const Vec &pVec) const
		{
			return Vector<3, T, P>(x + pVec.x, y + pVec.y, z + pVec.z);
		}

		inline Vec operator-(const Vec &pVec) const
		{
			return Vector<3, T, P>(x - pVec.x, y - pVec.y, z - pVec.z);
		}

		inline Vec operator*(const Vec &pVec) const
		{
			return Vector<3, T, P>(x * pVec.x, y * pVec.y, z * pVec.z);
		}

		inline Vec operator*(const T pVal) const
		{
			return Vector<3, T, P>(x * pVal, y * pVal, z * pVal);
		}

		inline Vec operator/(const Vec &pVec) const
		{
			return Vector<3, T, P>(x / pVec.x, y / pVec.y, z / pVec.z);
		}

		inline Vec operator/(const T pVal) const
		{
			return Vector<3, T, P>(x / pVal, y / pVal, z / pVal);
		}

		// -- Boolean operators --

		inline bool operator==(const Vec &pVec) const
		{
			return x == pVec.x && y == pVec.y && z == pVec.z;
		}

		inline bool operator!=(const Vec &pVec) const
		{
			return x != pVec.x || y != pVec.y || z != pVec.z;
		}

		// -- Convertion operators --

		template<typename U, typename Q>
		inline explicit operator Vector<3, U, Q>() const
		{
			return Vector<3, U, Q>((U)x, (U)y, (U)z);
		}

		// -- Stream operators --

		friend std::ostream &operator<<(std::ostream &pOStream, const Vec &pVec)
		{
			return pOStream << "(" << pVec.x << ", " << pVec.y << ", " << pVec.z << ")";
		}

		// -- Getters --

		Vec Abs() const
		{
			return Vector<3, T, P>(Math::Abs(x), Math::Abs(y), Math::Abs(z));
		}

		P AngleTo(const Vec &pVec) const
		{
			return std::acos(Math::Clamp(Dot(pVec), T(-1), T(1)));
		}

		Vec Cross(const Vec &pVec) const
		{
			return Vector<3, T, P>(
				y * pVec.z - z * pVec.y,
				z * pVec.x - x * pVec.z,
				x * pVec.y - y * pVec.x
			);
		}

		P Distance(const Vec &pVec) const
		{
			const P to_x = pVec.x - x;
			const P to_y = pVec.y - y;
			const P to_z = pVec.z - z;
			return std::sqrt(to_x * to_x + to_y * to_y + to_z * to_z);
		}

		P DistanceSquared(const Vec &pVec) const
		{
			const P to_x = pVec.x - x;
			const P to_y = pVec.y - y;
			const P to_z = pVec.z - z;
			return to_x * to_x + to_y * to_y + to_z * to_z;
		}

		T Dot(const Vec &pVec) const
		{
			return x * pVec.x + y * pVec.y + z * pVec.z;
		}

		bool IsNormalized() const
		{
			return x * x + y * y + z * z == T(1);
		}

		P Length() const
		{
			return std::sqrt(x * x + y * y + z * z);
		}

		T LengthSquared() const
		{
			return x * x + y * y + z * z;
		}

		Vec Normalized() const
		{
			P lLength = LengthSquared();

			if (lLength == 0 || lLength == 1)
			{
				return *this;
			}

			lLength = std::sqrt(lLength);
			return Vector<3, T, P>(x / lLength, y / lLength, z / lLength);
		}

		Vec Rotated(const P pAngle, const Vector<3, P, P> &pNormal)
		{
			Matrix<4, 4, T, P> lRot = Math::Rotate(Math::Matrix<4, 4, T, P>(T(1)), pAngle, pNormal);
			Vector<4, T, P> lRes = lRot * Vector<4, T, P>(x, y, z, T(1));
			return Vec(lRes.x, lRes.y, lRes.z);
		}

		Vec Sign() const
		{
			return Vector<3, T, P>(Math::Sign(x), Math::Sign(y), Math::Sign(z));
		}

		inline size_t SizeofField() const
		{
			return sizeof(T);
		}

		// -- Transformations --

		Vec &Normalize()
		{
			P lLength = LengthSquared();

			if (lLength == 0 || lLength == 1)
			{
				return *this;
			}

			lLength = std::sqrt(lLength);
			x /= lLength;
			y /= lLength;
			z /= lLength;
			return *this;
		}

		Vec &RotateAroundX(const P pAngle)
		{
			return RotateAroundX(std::cos(pAngle), std::sin(pAngle));
		}

		Vec &RotateAroundY(const P pAngle)
		{
			return RotateAroundY(std::cos(pAngle), std::sin(pAngle));
		}

		Vec &RotateAroundZ(const P pAngle)
		{
			return RotateAroundZ(std::cos(pAngle), std::sin(pAngle));
		}

		Vec &Rotate(const T pAngle, const Vector<3, T, P> &pNormal)
		{
			const Matrix<4, 4, T, P> lRot = Math::Rotate(Math::Matrix<4, 4, T, P>(T(1)), pAngle, pNormal);
			Vector<4, T, P> lRes = lRot * Vector<4, T, P>(x, y, z, T(1));
			x = lRes.x;
			y = lRes.y;
			z = lRes.z;
			return *this;
		}

		// -- Static getters --

		static Vec Zero()
		{
			return Vector<3, T, P>(T(0), T(0), T(0));
		}

		static Vec One()
		{
			return Vector<3, T, P>(T(1), T(1), T(1));
		}

		static Vec NegOne()
		{
			return Vector<3, T, P>(T(-1), T(-1), T(-1));
		}

		static Vec Left()
		{
			return Vector<3, T, P>(T(-1), T(0), T(0));
		}

		static Vec Right()
		{
			return Vector<3, T, P>(T(1), T(0), T(0));
		}

		static Vec Up()
		{
			return Vector<3, T, P>(T(0), T(1), T(0));
		}

		static Vec Down()
		{
			return Vector<3, T, P>(T(0), T(-1), T(0));
		}

		static Vec Forward()
		{
			return Vector<3, T, P>(T(0), T(0), T(1));
		}

		static Vec Backward()
		{
			return Vector<3, T, P>(T(0), T(0), T(-1));
		}

		// -- Static methods --

		static Vec Lerp(const Vec &pFrom, const Vec &pTo, const P pTime)
		{
			return Vector<3, T, P>(
				Math::Lerp(pFrom.x, pTo.x, pTime),
				Math::Lerp(pFrom.y, pTo.y, pTime),
				Math::Lerp(pFrom.z, pTo.z, pTime)
			);
		}

		static Vec LerpClamped(const Vec &pFrom, const Vec &pTo, const P pTime)
		{
			return Vector<3, T, P>(
				Math::LerpClamped(pFrom.x, pTo.x, pTime),
				Math::LerpClamped(pFrom.y, pTo.y, pTime),
				Math::LerpClamped(pFrom.z, pTo.z, pTime)
			);
		}
	};

	typedef Vector<3, float, float> Vector3;
	typedef Vector<3, int, float> Vector3I;
	typedef Vector<3, double, double> Vector3D;
}

