#pragma once

#include "Core.h"
#include "Vector.hpp"

namespace Math
{
	template<unsigned int C, unsigned int R, typename T, typename P>
	struct Matrix
	{
	private:
		typedef Matrix<C, R, T, P> MatCxR;
		typedef Matrix<R, C, T, P> MatTranspose;
		typedef Vector<R, T, P> Column;
		Column mValues[C];

	public:
		// -- Constructors --

		Matrix()
		{
			for (int i = 0; i < C; i++)
			{
				for (int j = 0; j < R; j++)
				{
					mValues[i][j] = T(0);
				}
			}
		}

		explicit Matrix(const T pVal)
		{
			for (int i = 0; i < C; i++)
			{
				for (int j = 0; j < R; j++)
				{
					mValues[i][j] = i == j ? pVal : T(0);
				}
			}
		}

		Matrix(const Column pVals[C])
		{
			for (int i = 0; i < C; i++)
			{
				mValues[i] = pVals[i];
			}
		}

		Matrix(const T pVals[C * R])
		{
			for (int i = 0; i < C; i++)
			{
				for (int j = 0; j < R; j++)
				{
					mValues[i][j] = pVals[i * R + C];
				}
			}
		}

		Matrix(const MatCxR &pMat)
		{
			for (int i = 0; i < C; i++)
			{
				mValues[i] = pMat[i];
			}
		}

		// -- Accesses --

		inline const Column &operator[](const int pIndex) const noexcept
		{
			return mValues[pIndex];
		}

		inline Column &operator[](const int pIndex) noexcept
		{
			return mValues[pIndex];
		}

		// -- Unary arithmetic operators --

		MatCxR &operator+=(const MatCxR &pMat)
		{
			for (int i = 0; i < C; i++)
			{
				mValues[i] += pMat[i];
			}

			return *this;
		}

		MatCxR &operator-=(const MatCxR &pMat)
		{
			for (int i = 0; i < C; i++)
			{
				mValues[i] -= pMat[i];
			}

			return *this;
		}

		MatCxR &operator*=(const T pVal)
		{
			for (int i = 0; i < C; i++)
			{
				mValues[i] *= pVal;
			}

			return *this;
		}

		MatCxR &operator/=(const T pVal)
		{
			for (int i = 0; i < C; i++)
			{
				mValues[i] /= pVal;
			}

			return *this;
		}

		// -- Unary operators --

		MatCxR operator+() const
		{
			return *this;
		}

		MatCxR operator-() const
		{
			MatCxR lRes;

			for (int i = 0; i < C; i++)
			{
				lRes[i] = -mValues[i];
			}

			return lRes;
		}

		// -- Binary operators --

		MatCxR operator+(const MatCxR &pMat) const
		{
			MatCxR lRes;

			for (int i = 0; i < C; i++)
			{
				lRes[i] = mValues[i] + pMat[i];
			}

			return lRes;
		}

		MatCxR operator-(const MatCxR &pMat) const
		{
			MatCxR lRes;

			for (int i = 0; i < C; i++)
			{
				lRes[i] = mValues[i] - pMat[i];
			}

			return lRes;
		}

		template<unsigned int S>
		Matrix<S, R, T, P> operator*(const Matrix<S, C, T, P> &lMat) const
		{
			Matrix<S, R, T, P> lRes;

			for (int i = 0; i < S; i++)
			{
				for (int j = 0; j < R; j++)
				{
					for (int k = 0; k < C; k++)
					{
						lRes[i][j] += mValues[k][j] * lMat[i][k];
					}
				}
			}

			return lRes;
		}

		MatCxR operator*(const T pVal) const
		{
			MatCxR lRes;

			for (int i = 0; i < C; i++)
			{
				lRes[i] = mValues[i] * pVal;
			}

			return lRes;
		}

		friend MatCxR operator*(const T pVal, const MatCxR &pMat)
		{
			MatCxR lRes;

			for (int i = 0; i < C; i++)
			{
				lRes[i] = pVal * pMat[i];
			}

			return lRes;
		}

		MatCxR operator/(const T pVal) const
		{
			MatCxR lRes;

			for (int i = 0; i < C; i++)
			{
				lRes[i] = mValues[i] / pVal;
			}

			return lRes;
		}

		// -- Boolean operators --

		bool operator==(const MatCxR &pMat) const
		{
			for (int i = 0; i < C; i++)
			{
				if (mValues[i] != pMat[i])
				{
					return false;
				}
			}

			return true;
		}

		bool operator!=(const MatCxR &pMat) const
		{
			for (int i = 0; i < C; i++)
			{
				if (mValues[i] == pMat[i])
				{
					return false;
				}
			}

			return true;
		}

		// -- Stream operators --

		friend std::ostream &operator<<(std::ostream &pOStream, const MatCxR &pMat)
		{
			for (int i = 0; i < C; i++)
			{
				for (int j = 0; j < R; j++)
				{
					pOStream << pMat[i][j] << (j < R - 1 ? ',' : '\n');
				}
			}

			return pOStream;
		}

		// -- Getters --
		MatTranspose Transposed() const
		{
			MatTranspose lRes;

			for (int i = 0; i < C; i++)
			{
				for (int j = 0; j < R; j++)
				{
					lRes[j][i] = mValues[i][j];
				}
			}

			return lRes;
		}
	};

	typedef Matrix<1, 1, float_type, float_type> Matrix1x1;
	typedef Matrix<1, 2, float_type, float_type> Matrix1x2;
	typedef Matrix<2, 1, float_type, float_type> Matrix2x1;
	typedef Matrix<1, 3, float_type, float_type> Matrix1x3;
	typedef Matrix<2, 3, float_type, float_type> Matrix2x3;
	typedef Matrix<3, 1, float_type, float_type> Matrix3x1;
	typedef Matrix<3, 2, float_type, float_type> Matrix3x2;
	typedef Matrix<1, 4, float_type, float_type> Matrix1x4;
	typedef Matrix<2, 4, float_type, float_type> Matrix2x4;
	typedef Matrix<3, 4, float_type, float_type> Matrix3x4;
	typedef Matrix<4, 1, float_type, float_type> Matrix4x1;
	typedef Matrix<4, 2, float_type, float_type> Matrix4x2;
	typedef Matrix<4, 3, float_type, float_type> Matrix4x3;
}

