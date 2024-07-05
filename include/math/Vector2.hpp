#pragma once

#include "Core.h"
#include "Math.hpp"
#include "Vector.hpp"
#include <cmath>

namespace Math
{
	template<typename T, typename P>
	struct Vector<2, T, P>
	{
	private:
		typedef Vector<2, T, P> Vec;

	public:
		union { T x, r, s; };
		union { T y, g, t; };

		// -- Constructors --

		Vector()
			: x(T(0)), y(T(0)) {}

		explicit Vector(const T pVal)
			: x(pVal), y(pVal) {}

		Vector(const T pX, const T pY)
			: x(pX), y(pY) {}

		Vector(const T pVals[2])
			: x(pVals[0]), y(pVals[1]) {}

		Vector(const Vec &pVec)
			: x(pVec.x), y(pVec.y) {}

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
			return *this;
		}

		inline Vec &operator-=(const Vec &pVec)
		{
			x -= pVec.x;
			y -= pVec.y;
			return *this;
		}

		inline Vec &operator*=(const Vec &pVec)
		{
			x *= pVec.x;
			y *= pVec.y;
			return *this;
		}

		inline Vec &operator*=(const T pVal)
		{
			x *= pVal;
			y *= pVal;
			return *this;
		}

		inline Vec &operator/=(const Vec &pVec)
		{
			x /= pVec.x;
			y /= pVec.y;
			return *this;
		}

		inline Vec &operator/=(const T pVal)
		{
			x /= pVal;
			y /= pVal;
			return *this;
		}

		// -- Unary operators --

		inline Vec operator+() const
		{
			return *this;
		}

		inline Vec operator-() const
		{
			return Vec(-x, -y);
		}

		// -- Binary operators --

		inline Vec operator+(const Vec &pVec) const
		{
			return Vec(x + pVec.x, y + pVec.y);
		}

		inline Vec operator-(const Vec &pVec) const
		{
			return Vec(x - pVec.x, y - pVec.y);
		}

		inline Vec operator*(const Vec &pVec) const
		{
			return Vec(x * pVec.x, y * pVec.y);
		}

		inline Vec operator*(const T pVal) const
		{
			return Vec(x * pVal, y * pVal);
		}

		inline friend Vec operator*(const T pVal, const Vec &pVec)
		{
			return Vec(pVal * pVec.x, pVal * pVec.y);
		}

		inline Vec operator/(const Vec &pVec) const
		{
			return Vec(x / pVec.x, y / pVec.y);
		}

		inline Vec operator/(const T pVal) const
		{
			return Vec(x / pVal, y / pVal);
		}

		// -- Boolean operators --

		inline bool operator==(const Vec &pVec) const
		{
			return x == pVec.x && y == pVec.y;
		}

		inline bool operator!=(const Vec &pVec) const
		{
			return x != pVec.x || y != pVec.y;
		}

		// -- Convertion operators --

		template<typename U, typename Q>
		inline explicit operator Vector<2, U, Q>() const
		{
			return Vector<2, U, Q>((U)x, (U)y);
		}

		// -- Stream operators --

		inline friend std::ostream &operator<<(std::ostream &pOStream, const Vec &pVec)
		{
			return pOStream << "(" << pVec.x << ", " << pVec.y << ")";
		}

		// -- Getters --

		Vec Abs() const
		{
			return Vec(Math::Abs(x), Math::Abs(y));
		}

		P AngleTo(const Vec &pVec) const
		{
			return std::acos(Math::Clamp(Dot(pVec), T(-1), T(1)));
		}

		T Cross(const Vec &pVec) const
		{
			return x * pVec.y - y * pVec.x;
		}

		P Distance(const Vec &pVec) const
		{
			const P lToX = pVec.x - x;
			const P lToY = pVec.y - y;
			const P lSqr = lToX * lToX + lToY * lToY;

			if (lSqr == P(0))
			{
				return P(0);
			}

			return std::sqrt(lSqr);
		}

		P DistanceSquared(const Vec &pVec) const
		{
			const P to_x = pVec.x - x;
			const P to_y = pVec.y - y;
			return to_x * to_x + to_y * to_y;
		}

		T Dot(const Vec &pVec) const
		{
			x * pVec.x + y * pVec.y;
		}

		bool IsNormalized() const
		{
			return x * x + y * y == T(1);
		}

		P Length() const
		{
			return std::sqrt(x * x + y * y);
		}

		T LengthSquared() const
		{
			return x * x + y * y;
		}

		Vec Normalized() const
		{
			P lLength = LengthSquared();

			if (lLength == 0 || lLength == 1)
			{
				return *this;
			}

			lLength = std::sqrt(lLength);
			return Vec(x / lLength, y / lLength);
		}

		Vec Rotated(const P pAngle) const
		{
			const P cos = std::cos(pAngle);
			const P sin = std::sin(pAngle);
			return Vec(cos * x - sin * y, sin * x + cos * y);
		}

		Vec Sign() const
		{
			return Vec(Math::Sign(x), Math::Sign(y));
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
			return *this;
		}

		Vec &Rotate(const P pAngle)
		{
			const P lCos = std::cos(pAngle);
			const P lSin = std::sin(pAngle);
			const T lTempX = x;
			const T lTempY = y;

			x = lCos * lTempX - lSin * lTempY;
			y = lSin * lTempX + lCos * lTempY;
			return *this;
		}

		// -- Static getters --

		inline static Vec Zero()
		{
			return Vec(T(0), T(0));
		}

		inline static Vec One()
		{
			return Vec(T(1), T(1));
		}

		inline static Vec NegOne()
		{
			return Vec(T(-1), T(-1));
		}

		inline static Vec Left()
		{
			return Vec(T(-1), T(0));
		}

		inline static Vec Right()
		{
			return Vec(T(1), T(0));
		}

		inline static Vec Up()
		{
			return Vec(T(0), T(1));
		}

		inline static Vec Down()
		{
			return Vec(T(0), T(-1));
		}

		// -- Static methods --

		static Vec Lerp(const Vec &pFrom, const Vec &pTo, const P pTime)
		{
			return Vec(
				Math::Lerp(pFrom.x, pTo.x, pTime),
				Math::Lerp(pFrom.y, pTo.y, pTime)
			);
		}

		static Vec LerpClamped(const Vec &pFrom, const Vec &pTo, const P pTime)
		{
			return Vec(
				Math::LerpClamped(pFrom.x, pTo.x, pTime),
				Math::LerpClamped(pFrom.y, pTo.y, pTime)
			);
		}
	};

	template<typename T = float_type, typename P = T>
	using Vector2 = Vector<2, T, P>;

	typedef Vector<2, float, float> Vector2F;
	typedef Vector<2, int, float> Vector2I;
	typedef Vector<2, double, double> Vector2D;
}

