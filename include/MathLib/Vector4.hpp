#pragma once

#include "Math.hpp"
#include "Vector.hpp"

namespace Math
{
	template<typename T, typename P>
	struct Vector<4, T, P>
	{
	private:
		typedef Vector<4, T, P> Vec;

	public:
		union { T x, r, s; };
		union { T y, g, t; };
		union { T z, b, p; };
		union { T w, a, q; };

		// -- Constructors --

		Vector()
			: x(T(0)), y(T(0)), z(T(0)), w(T(0)) {}

		explicit Vector(const T pVal)
			: x(pVal), y(pVal), z(pVal), w(pVal) {}

		Vector(const T pX, const T pY, const T pZ, const T pW)
			: x(pX), y(pY), z(pZ), w(pW) {}

		Vector(const T pVals[4])
			: x(pVals[0]), y(pVals[1]), z(pVals[2]), w(pVals[3]) {}

		Vector(const Vec &pVec)
			: x(pVec.x), y(pVec.y), z(pVec.z), w(pVec.w) {}

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
			w += pVec.w;
			return *this;
		}

		inline Vec &operator-=(const Vec &pVec)
		{
			x -= pVec.x;
			y -= pVec.y;
			z -= pVec.z;
			w -= pVec.w;
			return *this;
		}

		inline Vec &operator*=(const Vec &pVec)
		{
			x *= pVec.x;
			y *= pVec.y;
			z *= pVec.z;
			w *= pVec.w;
			return *this;
		}

		inline Vec &operator*=(const T pVal)
		{
			x *= pVal;
			y *= pVal;
			z *= pVal;
			w *= pVal;
			return *this;
		}

		inline Vec &operator/=(const Vec &pVec)
		{
			x /= pVec.x;
			y /= pVec.y;
			z /= pVec.z;
			w /= pVec.w;
			return *this;
		}

		inline Vec &operator/=(const T pVal)
		{
			x /= pVal;
			y /= pVal;
			z /= pVal;
			w /= pVal;
			return *this;
		}

		// -- Unary operators --

		inline Vec operator+() const
		{
			return *this;
		}

		inline Vec operator-() const
		{
			return Vec(-x, -y, -z, -w);
		}

		// -- Binary operators --

		inline Vec operator+(const Vec &pVec) const
		{
			return Vec(x + pVec.x, y + pVec.y, z + pVec.z, w + pVec.w);
		}

		inline Vec operator-(const Vec &pVec) const
		{
			return Vec(x - pVec.x, y - pVec.y, z - pVec.z, w - pVec.w);
		}

		inline Vec operator*(const Vec &pVec) const
		{
			return Vec(x * pVec.x, y * pVec.y, z * pVec.z, w * pVec.w);
		}

		inline Vec operator*(const T pVal) const
		{
			return Vec(x * pVal, y * pVal, z * pVal, w * pVal);
		}

		inline friend Vec operator*(const T pVal, const Vec &pVec)
		{
			return Vec(pVal * pVec.x, pVal * pVec.y, pVal * pVec.z, pVal * pVec.w);
		}

		inline Vec operator/(const Vec &pVec) const
		{
			return Vec(x / pVec.x, y / pVec.y, z / pVec.z, w / pVec.w);
		}

		inline Vec operator/(const T pVal) const
		{
			return Vec(x / pVal, y / pVal, z / pVal, w / pVal);
		}

		// -- Boolean operators --

		inline bool operator==(const Vec &pVec) const
		{
			return x == pVec.x && y == pVec.y && z == pVec.z && w == pVec.w;
		}

		inline bool operator!=(const Vec &pVec) const
		{
			return x != pVec.x || y != pVec.y || z != pVec.z || w != pVec.w;
		}

		// -- Convertion operators --

		template<typename U, typename Q>
		inline explicit operator Vector<4, U, Q>() const
		{
			return Vector<4, U, Q>((U)x, (U)y, (U)z, (U)w);
		}

		// -- Stream operators --

		inline friend std::ostream &operator<<(std::ostream &pOStream, const Vec &pVec)
		{
			return pOStream << "(" << pVec.x << ", " << pVec.y << ", " << pVec.z << ", " << pVec.w << ")";
		}

		// -- Getters --

		Vec Abs() const
		{
			return Vec(Math::Abs(x), Math::Abs(y), Math::Abs(z), Math::Abs(w));
		}

		P AngleTo(const Vec &pVec) const
		{
			return std::acos(Math::Clamp(Dot(pVec), T(-1), T(1)));
		}

		P Distance(const Vec &pVec) const
		{
			const P lToX = pVec.x - x;
			const P lToY = pVec.y - y;
			const P lToZ = pVec.z - z;
			const P lToW = pVec.w - w;
			const P lSqr = lToX * lToX + lToY * lToY + lToZ * lToZ + lToW * lToW;

			if (lSqr == P(0))
			{
				return P(0);
			}

			return std::sqrt(lSqr);
		}

		P DistanceSquared(const Vec &pVec) const
		{
			const P lToX = pVec.x - x;
			const P lToY = pVec.y - y;
			const P lToZ = pVec.z - z;
			const P lToW = pVec.w - w;
			return lToX * lToX + lToY * lToY + lToZ * lToZ + lToW * lToW;
		}

		T Dot(const Vec &pVec) const
		{
			return x * pVec.x + y * pVec.y + z * pVec.z + w * pVec.w;
		}

		bool IsNormalized() const
		{
			return x * x + y * y + z * z + w * w == T(1);
		}

		P Length() const
		{
			return std::sqrt(x * x + y * y + z * z + w * w);
		}

		T LengthSquared() const
		{
			return x * x + y * y + z * z + w * w;
		}

		Vec Normalized() const
		{
			P lLength = LengthSquared();

			if (lLength == 0 || lLength == 1)
			{
				return *this;
			}

			lLength = std::sqrt(lLength);
			return Vec(x / lLength, y / lLength, z / lLength, w / lLength);
		}

		Vec Sign() const
		{
			return Vec(Math::Sign(x), Math::Sign(y), Math::Sign(z), Math::Sign(w));
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
			w /= lLength;
			return *this;
		}

		// -- Static getters --

		inline static Vec Zero()
		{
			return Vec(T(0), T(0), T(0), T(0));
		}

		inline static Vec One()
		{
			return Vec(T(1), T(1), T(1), T(1));
		}

		inline static Vec NegOne()
		{
			return Vec(T(-1), T(-1), T(-1), T(-1));
		}

		inline static Vec PosX()
		{
			return Vec(T(1), T(0), T(0), T(0));
		}

		inline static Vec NegX()
		{
			return Vec(T(-1), T(0), T(0), T(0));
		}

		inline static Vec PosY()
		{
			return Vec(T(0), T(1), T(0), T(0));
		}

		inline static Vec NegY()
		{
			return Vec(T(0), T(-1), T(0), T(0));
		}

		inline static Vec PosZ()
		{
			return Vec(T(0), T(0), T(1), T(0));
		}

		inline static Vec NegZ()
		{
			return Vec(T(0), T(0), T(-1), T(0));
		}

		inline static Vec PosW()
		{
			return Vec(T(0), T(0), T(0), T(1));
		}

		inline static Vec NegW()
		{
			return Vec(T(0), T(0), T(0), T(-1));
		}

		// -- Static methods --

		static Vec Lerp(const Vec &pFrom, const Vec &pTo, const P pTime)
		{
			return Vec(
				Math::Lerp(pFrom.x, pTo.x, pTime),
				Math::Lerp(pFrom.y, pTo.y, pTime),
				Math::Lerp(pFrom.z, pTo.z, pTime),
				Math::Lerp(pFrom.w, pTo.w, pTime)
			);
		}

		static Vec LerpClamped(const Vec &pFrom, const Vec &pTo, const P pTime)
		{
			return Vec(
				Math::LerpClamped(pFrom.x, pTo.x, pTime),
				Math::LerpClamped(pFrom.y, pTo.y, pTime),
				Math::LerpClamped(pFrom.z, pTo.z, pTime),
				Math::LerpClamped(pFrom.w, pTo.w, pTime)
			);
		}
	};

	typedef Vector<4, float, float> Vector4;
	typedef Vector<4, int, float> Vector4I;
	typedef Vector<4, double, double> Vector4D;
}

