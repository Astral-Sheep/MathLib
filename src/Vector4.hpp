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

		Vector<4, T, P>()
			: x(T(0)), y(T(0)), z(T(0)), w(T(0)) {}

		explicit Vector<4, T, P>(const T pVal)
			: x(pVal), y(pVal), z(pVal), w(pVal) {}

		Vector <4, T, P>(const T pX, const T pY, const T pZ, const T pW)
			: x(pX), y(pY), z(pZ), w(pW) {}

		Vector<4, T, P>(const T pVals[4])
			: x(pVals[0]), y(pVals[1]), z(pVals[2]), w(pVals[3]) {}

		Vector<4, T, P>(const Vec &pVec)
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
			return Vector<4, T, P>(-x, -y, -z, -w);
		}

		// -- Binary operators --

		inline Vec operator+(const Vec &pVec) const
		{
			return Vector<4, T, P>(x + pVec.x, y + pVec.y, z + pVec.z, w + pVec.w);
		}

		inline Vec operator-(const Vec &pVec) const
		{
			return Vector<4, T, P>(x - pVec.x, y - pVec.y, z - pVec.z, w - pVec.w);
		}

		inline Vec operator*(const Vec &pVec) const
		{
			return Vector<4, T, P>(x * pVec.x, y * pVec.y, z * pVec.z, w * pVec.w);
		}

		inline Vec operator*(const T pVal) const
		{
			return Vector<4, T, P>(x * pVal, y * pVal, z * pVal, w * pVal);
		}

		inline Vec operator/(const Vec &pVec) const
		{
			return Vector<4, T, P>(x / pVec.x, y / pVec.y, z / pVec.z, w / pVec.w);
		}

		inline Vec operator/(const T pVal) const
		{
			return Vector<4, T, P>(x / pVal, y / pVal, z / pVal, w / pVal);
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

		friend std::ostream &operator<<(std::ostream &pOStream, const Vec &pVec)
		{
			return pOStream << "(" << pVec.x << ", " << pVec.y << ", " << pVec.z << ", " << pVec.w << ")";
		}

		// -- Getters --

		Vec Abs() const
		{
			return Vector<4, T, P>(Math::Abs(x), Math::Abs(y), Math::Abs(z), Math::Abs(w));
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
			return std::sqrt(lToX * lToX + lToY * lToY + lToZ * lToZ + lToW * lToW);
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
			return Vector<4, T, P>(x / lLength, y / lLength, z / lLength, w / lLength);
		}

		Vec Sign() const
		{
			return Vector<4, T, P>(Math::Sign(x), Math::Sign(y), Math::Sign(z), Math::Sign(w));
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

		static inline Vec Zero()
		{
			return Vector<4, T, P>(T(0), T(0), T(0), T(0));
		}

		static inline Vec One()
		{
			return Vector<4, T, P>(T(1), T(1), T(1), T(1));
		}

		static inline Vec NegOne()
		{
			return Vector<4, T, P>(T(-1), T(-1), T(-1), T(-1));
		}

		static inline Vec PosX()
		{
			return Vector<4, T, P>(T(1), T(0), T(0), T(0));
		}

		static inline Vec NegX()
		{
			return Vector<4, T, P>(T(-1), T(0), T(0), T(0));
		}

		static inline Vec PosY()
		{
			return Vector<4, T, P>(T(0), T(1), T(0), T(0));
		}

		static inline Vec NegY()
		{
			return Vector<4, T, P>(T(0), T(-1), T(0), T(0));
		}

		static inline Vec PosZ()
		{
			return Vector<4, T, P>(T(0), T(0), T(1), T(0));
		}

		static inline Vec NegZ()
		{
			return Vector<4, T, P>(T(0), T(0), T(-1), T(0));
		}

		static inline Vec PosW()
		{
			return Vector<4, T, P>(T(0), T(0), T(0), T(1));
		}

		static inline Vec NegW()
		{
			return Vector<4, T, P>(T(0), T(0), T(0), T(-1));
		}

		// -- Static methods --

		static Vec Lerp(const Vec &pFrom, const Vec &pTo, const P pTime)
		{
			return Vector<4, T, P>(
				Math::Lerp(pFrom.x, pTo.x, pTime),
				Math::Lerp(pFrom.y, pTo.y, pTime),
				Math::Lerp(pFrom.z, pTo.z, pTime),
				Math::Lerp(pFrom.w, pTo.w, pTime)
			);
		}

		static Vec LerpClamped(const Vec &pFrom, const Vec &pTo, const P pTime)
		{
			return Vector<4, T, P>(
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

