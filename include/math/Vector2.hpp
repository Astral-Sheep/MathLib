#pragma once

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

		Vector<2, T, P>()
			: x(T(0)), y(T(0)) {}

		explicit Vector<2, T, P>(const T pVal)
			: x(pVal), y(pVal) {}

		Vector<2, T, P>(const T pX, const T pY)
			: x(pX), y(pY) {}

		Vector<2, T, P>(const T pVals[2])
			: x(pVals[0]), y(pVals[1]) {}

		Vector<2, T, P>(const Vec &pVec)
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
			return Vector<2, T, P>(-x, -y);
		}

		// -- Binary operators --

		inline Vec operator+(const Vec &pVec) const
		{
			return Vector<2, T, P>(x + pVec.x, y + pVec.y);
		}

		inline Vec operator-(const Vec &pVec) const
		{
			return Vector<2, T, P>(x - pVec.x, y - pVec.y);
		}

		inline Vec operator*(const Vec &pVec) const
		{
			return Vector<2, T, P>(x * pVec.x, y * pVec.y);
		}

		inline Vec operator*(const T pVal) const
		{
			return Vector<2, T, P>(x * pVal, y * pVal);
		}

		inline Vec operator/(const Vec &pVec) const
		{
			return Vector<2, T, P>(x / pVec.x, y / pVec.y);
		}

		inline Vec operator/(const T pVal) const
		{
			return Vector<2, T, P>(x / pVal, y / pVal);
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

		friend std::ostream &operator<<(std::ostream &pOStream, const Vec &pVec)
		{
			return pOStream << "(" << pVec.x << ", " << pVec.y << ")";
		}

		// -- Getters --

		Vec Abs() const
		{
			return Vector<2, T, P>(Math::Abs(x), Math::Abs(y));
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
			const P to_x = pVec.x - x;
			const P to_y = pVec.y - y;
			return std::sqrt(to_x * to_x - to_y * to_y);
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
			return Vector<2, T, P>(x / lLength, y / lLength);
		}

		Vec Rotated(const P pAngle) const
		{
			const P cos = std::cos(pAngle);
			const P sin = std::sin(pAngle);
			return Vector<2, T, P>(cos * x - sin * y, sin * x + cos * y);
		}

		Vec Sign() const
		{
			return Vector<2, T, P>(Math::Sign(x), Math::Sign(y));
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

		static Vec Zero()
		{
			return Vector<2, T, P>(T(0), T(0));
		}

		static Vec One()
		{
			return Vector<2, T, P>(T(1), T(1));
		}

		static Vec NegOne()
		{
			return Vector<2, T, P>(T(-1), T(-1));
		}

		static Vec Left()
		{
			return Vector<2, T, P>(T(-1), T(0));
		}

		static Vec Right()
		{
			return Vector<2, T, P>(T(1), T(0));
		}

		static Vec Up()
		{
			return Vector<2, T, P>(T(0), T(1));
		}

		static Vec Down()
		{
			return Vector<2, T, P>(T(0), T(-1));
		}

		// -- Static methods --

		static Vec Lerp(const Vec &pFrom, const Vec &pTo, const P pTime)
		{
			return Vector<2, T, P>(
				Math::Lerp(pFrom.x, pTo.x, pTime),
				Math::Lerp(pFrom.y, pTo.y, pTime)
			);
		}

		static Vec LerpClamped(const Vec &pFrom, const Vec &pTo, const P pTime)
		{
			return Vector<2, T, P>(
				Math::LerpClamped(pFrom.x, pTo.x, pTime),
				Math::LerpClamped(pFrom.y, pTo.y, pTime)
			);
		}
	};

	typedef Vector<2, float, float> Vector2;
	typedef Vector<2, int, float> Vector2I;
	typedef Vector<2, double, double> Vector2D;
}

