#pragma once

#include "Math.hpp"
#include <cmath>
#include <ostream>

namespace Math
{
	template<unsigned int S, typename T, typename P>
	struct Vector
	{
	private:
		typedef Vector<S, T, P> Vec;
		T mValues[S];

	public:
		// -- Constructors --

		Vector()
		{
			for (int i = 0; i < S; i++)
			{
				mValues[i] = T(0);
			}
		}

		explicit Vector(const T pVal)
		{
			for (int i = 0; i < S; i++)
			{
				mValues[i] = pVal;
			}
		}

		Vector(const T pVals[S])
		{
			for (int i = 0; i < S; i++)
			{
				mValues[i] = pVals[i];
			}
		}

		Vector(const Vec &pVec)
		{
			for (int i = 0; i < S; i++)
			{
				mValues[i] = pVec[i];
			}
		}

		// -- Accesses --

		inline const T &operator[](const int pIndex) const noexcept
		{
			return mValues[pIndex];
		}

		inline T &operator[](const int pIndex) noexcept
		{
			return mValues[pIndex];
		}

		// -- Unary arithmetic operators --

		Vec &operator+=(const Vec &pVec)
		{
			for (int i = 0; i < S; i++)
			{
				mValues[i] += pVec[i];
			}

			return *this;
		}

		Vec &operator-=(const Vec &pVec)
		{
			for (int i = 0; i < S; i++)
			{
				mValues[i] -= pVec[i];
			}

			return *this;
		}

		Vec &operator*=(const Vec &pVec)
		{
			for (int i = 0; i < S; i++)
			{
				mValues[i] *= pVec[i];
			}

			return *this;
		}

		Vec &operator*=(const T pVal)
		{
			for (int i = 0; i < S; i++)
			{
				mValues[i] *= pVal;
			}

			return *this;
		}

		Vec &operator/=(const Vec &pVec)
		{
			for (int i = 0; i < S; i++)
			{
				mValues[i] /= pVec[i];
			}

			return *this;
		}

		Vec &operator/=(const T pVal)
		{
			for (int i = 0; i < S; i++)
			{
				mValues[i] /= pVal;
			}

			return *this;
		}

		// -- Unary operators --

		Vec operator+() const
		{
			return *this;
		}

		Vec operator-() const
		{
			Vec lRes;

			for (int i = 0; i < S; i++)
			{
				lRes[i] = -mValues[i];
			}

			return lRes;
		}

		// -- Binary operators --

		Vec operator+(const Vec &pVec) const
		{
			Vec lRes;

			for (int i = 0; i < S; i++)
			{
				lRes[i] = mValues[i] + pVec[i];
			}

			return lRes;
		}

		Vec operator-(const Vec &pVec) const
		{
			Vec lRes;

			for (int i = 0; i < S; i++)
			{
				lRes[i] = mValues[i] - pVec[i];
			}

			return lRes;
		}

		Vec operator*(const Vec &pVec) const
		{
			Vec lRes;

			for (int i = 0; i < S; i++)
			{
				lRes[i] = mValues[i] * pVec[i];
			}

			return lRes;
		}

		Vec operator*(const T pVal) const
		{
			Vec lRes;

			for (int i = 0; i < S; i++)
			{
				lRes[i] = mValues[i] * pVal;
			}

			return lRes;
		}

		friend Vec operator*(const T pVal, const Vec &pVec)
		{
			Vec lRes;

			for (int i = 0; i < S; i++)
			{
				lRes[i] = pVal * pVec[i];
			}

			return lRes;
		}

		Vec operator/(const Vec &pVec) const
		{
			Vec lRes;

			for (int i = 0; i < S; i++)
			{
				lRes[i] = mValues[i] / pVec[i];
			}

			return lRes;
		}

		Vec operator/(const T pVal) const
		{
			Vec lRes;

			for (int i = 0; i < S; i++)
			{
				lRes[i] = mValues[i] / pVal;
			}

			return lRes;
		}

		// -- Boolean operators --

		bool operator==(const Vec &pVec) const
		{
			for (int i = 0; i < S; i++)
			{
				if (mValues[i] != pVec[i])
				{
					return false;
				}
			}

			return true;
		}

		bool operator!=(const Vec &pVec) const
		{
			for (int i = 0; i < S; i++)
			{
				if (mValues[i] == pVec[i])
				{
					return false;
				}
			}

			return true;
		}

		// -- Convertion operators --

		template<typename U, typename Q>
		operator Vector<S, U, Q>() const
		{
			Vector<S, U, Q> lRes;

			for (int i = 0; i < S; i++)
			{
				lRes[i] = (U)mValues[i];
			}

			return lRes;
		}

		// -- Stream operators --

		friend std::ostream &operator<<(std::ostream &pOStream, const Vec &pVec)
		{
			for (int i = 0; i < S; i++)
			{
				pOStream << pVec[i];

				if (i < S - 1)
				{
					pOStream << ',';
				}
			}

			return pOStream;
		}

		// -- Static methods --

		static Vec Lerp(const Vec &pFrom, const Vec &pTo, const float pTime)
		{
			return pFrom + (pTo - pFrom) * pTime;
		}

		static Vec LerpClamped(const Vec &pFrom, const Vec &pTo, const float pTime)
		{
			return pFrom + (pTo - pFrom) * Math::Clamp01(pTime);
		}

		// -- Getters --

		Vec Abs() const
		{
			Vec lRes;

			for (int i = 0; i < S; i++)
			{
				lRes[i] = Math::Abs(mValues[i]);
			}

			return lRes;
		}

		P Angle(const Vec &pVec) const
		{
			return std::acos(Math::Clamp(Dot(pVec), T(-1), T(1)));
		}

		P Distance(const Vec &pVec) const
		{
			P lSqrdist = P(0);

			for (int i = 0; i < S; i++)
			{
				lSqrdist += (pVec[i] - mValues[i]) * (pVec[i] - mValues[i]);
			}

			if (lSqrdist == T(0))
			{
				return T(0);
			}

			return std::sqrt(lSqrdist);
		}

		P DistanceSquared(const Vec &pVec) const
		{
			P lSqrdist = P(0);

			for (int i = 0; i < S; i++)
			{
				lSqrdist += (pVec[i] - mValues[i]) * (pVec[i] - mValues[i]);
			}

			return lSqrdist;
		}

		P Dot(const Vec &pVec) const
		{
			P lDot = P(0);

			for (int i = 0; i < S; i++)
			{
				lDot = mValues[i] * pVec[i];
			}

			return lDot;
		}

		inline bool IsNormalized() const
		{
			return LengthSquared() == P(1);
		}

		P Length() const
		{
			P lSqrLength = P(0);

			for (int i = 0; i < S; i++)
			{
				lSqrLength += mValues[i] * mValues[i];
			}

			return std::sqrt(lSqrLength);
		}

		P LengthSquared() const
		{
			P lSqrLength = P(0);

			for (int i = 0; i < S; i++)
			{
				lSqrLength += mValues[i] * mValues[i];
			}

			return lSqrLength;
		}

		Vec Normalized() const
		{
			P lLength = LengthSquared();

			if (lLength == 0 || lLength == 1)
			{
				return *this;
			}

			Vec lRes;
			lLength = std::sqrt(lLength);

			for (int i = 0; i < S; i++)
			{
				lRes[i] = mValues[i] / lLength;
			}

			return lRes;
		}

		Vec Sign() const
		{
			Vec lRes;

			for (int i = 0; i < S; i++)
			{
				lRes[i] = Math::Sign(mValues[i]);
			}

			return lRes;
		}

		inline size_t SizeOfField() const
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

			for (int i = 0; i < S; i++)
			{
				mValues[i] /= lLength;
			}

			return *this;
		}
	};

	template<unsigned int S>
	using VectorI = Vector<S, int, float>;

	template<unsigned int S>
	using VectorF = Vector<S, float, float>;

	template<unsigned int S>
	using VectorD = Vector<S, double, double>;
}

