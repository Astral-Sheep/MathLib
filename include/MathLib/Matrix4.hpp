#pragma once

#include "Vector4.hpp"
#include "Matrix.hpp"

namespace Math
{
	template<typename T, typename P>
	struct Matrix<4, 4, T, P>
	{
	private:
		typedef Matrix<4, 4, T, P> Mat;
		typedef Vector<4, T, P> Col;
		Col mValues[4];

	public:
		// -- Constructors --

		Matrix()
			: mValues{ Col(), Col(), Col(), Col() } {}

		explicit Matrix(const T pVal)
			: mValues{
				Col(pVal, T(0), T(0), T(0)),
				Col(T(0), pVal, T(0), T(0)),
				Col(T(0), T(0), pVal, T(0)),
				Col(T(0), T(0), T(0), pVal)
			} {}

		Matrix(
			const T pX0, const T pX1, const T pX2, const T pX3,
			const T pX4, const T pX5, const T pX6, const T pX7,
			const T pX8, const T pX9, const T pX10, const T pX11,
			const T pX12, const T pX13, const T pX14, const T pX15
		)
			: mValues{
				Col(pX0, pX4, pX8, pX12),
				Col(pX1, pX5, pX9, pX13),
				Col(pX2, pX6, pX10, pX14),
				Col(pX3, pX7, pX11, pX15)
			} {}

		Matrix(const T pVals[16])
			: mValues{
				Col(pVals[0], pVals[4], pVals[8], pVals[12]),
				Col(pVals[1], pVals[5], pVals[9], pVals[13]),
				Col(pVals[2], pVals[6], pVals[10], pVals[14]),
				Col(pVals[3], pVals[7], pVals[11], pVals[15])
			} {}

		Matrix(const Col &pX0, const Col &pX1, const Col &pX2, const Col &pX3)
			: mValues{
				Col(pX0[0], pX1[0], pX2[0], pX3[0]),
				Col(pX0[1], pX1[1], pX2[1], pX3[1]),
				Col(pX0[2], pX1[2], pX2[2], pX3[2]),
				Col(pX0[3], pX1[3], pX2[3], pX3[3]),
			} {}

		Matrix(const Col pVals[4])
			: mValues{
				Col(pVals[0][0], pVals[1][0], pVals[2][0], pVals[3][0]),
				Col(pVals[0][1], pVals[1][1], pVals[2][1], pVals[3][1]),
				Col(pVals[0][2], pVals[1][2], pVals[2][2], pVals[3][2]),
				Col(pVals[0][3], pVals[1][3], pVals[2][3], pVals[3][3])
			} {}

		Matrix (const Mat &pMat)
			: mValues{ pMat[0], pMat[1], pMat[2], pMat[3] } {}

		// -- Accesses --

		inline const Col &operator[](const int pIndex) const noexcept
		{
			return mValues[pIndex];
		}

		inline Col &operator[](const int pIndex) noexcept
		{
			return mValues[pIndex];
		}

		// -- Unary arithmetic operators --

		inline Mat &operator+=(const Mat &pMat)
		{
			mValues[0] += pMat[0];
			mValues[1] += pMat[1];
			mValues[2] += pMat[2];
			mValues[3] += pMat[3];
			return *this;
		}

		inline Mat &operator-=(const Mat &pMat)
		{
			mValues[0] -= pMat[0];
			mValues[1] -= pMat[1];
			mValues[2] -= pMat[2];
			mValues[3] -= pMat[3];
			return *this;
		}

		Mat &operator*=(const Mat &pMat)
		{
			const Mat lTemp(*this);

			mValues[0][0] = lTemp[0][0] * pMat[0][0] + lTemp[1][0] * pMat[0][1] + lTemp[2][0] * pMat[0][2] + lTemp[3][0] * pMat[0][3];
			mValues[0][1] = lTemp[0][1] * pMat[0][0] + lTemp[1][1] * pMat[0][1] + lTemp[2][1] * pMat[0][2] + lTemp[3][1] * pMat[0][3];
			mValues[0][2] = lTemp[0][2] * pMat[0][0] + lTemp[1][2] * pMat[0][1] + lTemp[2][2] * pMat[0][2] + lTemp[3][2] * pMat[0][3];
			mValues[0][3] = lTemp[0][3] * pMat[0][0] + lTemp[1][3] * pMat[0][1] + lTemp[2][3] * pMat[0][2] + lTemp[3][3] * pMat[0][3];

			mValues[1][0] = lTemp[0][0] * pMat[1][0] + lTemp[1][0] * pMat[1][1] + lTemp[2][0] * pMat[1][2] + lTemp[3][0] * pMat[1][3];
			mValues[1][1] = lTemp[0][1] * pMat[1][0] + lTemp[1][1] * pMat[1][1] + lTemp[2][1] * pMat[1][2] + lTemp[3][1] * pMat[1][3];
			mValues[1][2] = lTemp[0][2] * pMat[1][0] + lTemp[1][2] * pMat[1][1] + lTemp[2][2] * pMat[1][2] + lTemp[3][2] * pMat[1][3];
			mValues[1][3] = lTemp[0][3] * pMat[1][0] + lTemp[1][3] * pMat[1][1] + lTemp[2][3] * pMat[1][2] + lTemp[3][3] * pMat[1][3];

			mValues[2][0] = lTemp[0][0] * pMat[2][0] + lTemp[1][0] * pMat[2][1] + lTemp[2][0] * pMat[2][2] + lTemp[3][0] * pMat[2][3];
			mValues[2][1] = lTemp[0][1] * pMat[2][0] + lTemp[1][1] * pMat[2][1] + lTemp[2][1] * pMat[2][2] + lTemp[3][1] * pMat[2][3];
			mValues[2][2] = lTemp[0][2] * pMat[2][0] + lTemp[1][2] * pMat[2][1] + lTemp[2][2] * pMat[2][2] + lTemp[3][2] * pMat[2][3];
			mValues[2][3] = lTemp[0][3] * pMat[2][0] + lTemp[1][3] * pMat[2][1] + lTemp[2][3] * pMat[2][2] + lTemp[3][3] * pMat[2][3];

			mValues[3][0] = lTemp[0][0] * pMat[3][0] + lTemp[1][0] * pMat[3][1] + lTemp[2][0] * pMat[3][2] + lTemp[3][0] * pMat[3][3];
			mValues[3][1] = lTemp[0][1] * pMat[3][0] + lTemp[1][1] * pMat[3][1] + lTemp[2][1] * pMat[3][2] + lTemp[3][1] * pMat[3][3];
			mValues[3][2] = lTemp[0][2] * pMat[3][0] + lTemp[1][2] * pMat[3][1] + lTemp[2][2] * pMat[3][2] + lTemp[3][2] * pMat[3][3];
			mValues[3][3] = lTemp[0][3] * pMat[3][0] + lTemp[1][3] * pMat[3][1] + lTemp[2][3] * pMat[3][2] + lTemp[3][3] * pMat[3][3];

			return *this;
		}

		inline Mat &operator*=(const T pVal)
		{
			mValues[0] *= pVal;
			mValues[1] *= pVal;
			mValues[2] *= pVal;
			mValues[3] *= pVal;
			return *this;
		}

		inline Mat &operator/=(const Mat &pMat)
		{
			return *this *= pMat.Inverted();
		}

		inline Mat &operator/=(const T pVal)
		{
			const T inv = T(1) / pVal;
			mValues[0] *= inv;
			mValues[1] *= inv;
			mValues[2] *= inv;
			mValues[3] *= inv;
			return *this;
		}

		// -- Unary operators --

		inline Mat operator+() const
		{
			return *this;
		}

		inline Mat operator-() const
		{
			return Mat(
				-mValues[0][0], -mValues[1][0], -mValues[2][0], -mValues[3][0],
				-mValues[0][1], -mValues[1][1], -mValues[2][1], -mValues[3][1],
				-mValues[0][2], -mValues[1][2], -mValues[2][2], -mValues[3][2],
				-mValues[0][3], -mValues[1][3], -mValues[2][3], -mValues[3][3]
			);
		}

		// -- Binary operators --

		inline Mat operator+(const Mat &pMat) const
		{
			return Mat(
				mValues[0][0] + pMat[0][0], mValues[1][0] + pMat[1][0], mValues[2][0] + pMat[2][0], mValues[3][0] + pMat[3][0],
				mValues[0][1] + pMat[0][1], mValues[1][1] + pMat[1][1], mValues[2][1] + pMat[2][1], mValues[3][1] + pMat[3][1],
				mValues[0][2] + pMat[0][2], mValues[1][2] + pMat[1][2], mValues[2][2] + pMat[2][2], mValues[3][2] + pMat[3][2],
				mValues[0][3] + pMat[0][3], mValues[1][3] + pMat[1][3], mValues[2][3] + pMat[2][3], mValues[3][3] + pMat[3][3]
			);
		}

		inline Mat operator-(const Mat &pMat) const
		{
			return Mat(
				mValues[0][0] - pMat[0][0], mValues[1][0] - pMat[1][0], mValues[2][0] - pMat[2][0], mValues[3][0] - pMat[3][0],
				mValues[0][1] - pMat[0][1], mValues[1][1] - pMat[1][1], mValues[2][1] - pMat[2][1], mValues[3][1] - pMat[3][1],
				mValues[0][2] - pMat[0][2], mValues[1][2] - pMat[1][2], mValues[2][2] - pMat[2][2], mValues[3][2] - pMat[3][2],
				mValues[0][3] - pMat[0][3], mValues[1][3] - pMat[1][3], mValues[2][3] - pMat[2][3], mValues[3][3] - pMat[3][3]
			);
		}

		inline Mat operator*(const Mat &pMat) const
		{
			return Mat(
				// Row 0
				mValues[0][0] * pMat[0][0] + mValues[1][0] * pMat[0][1] + mValues[2][0] * pMat[0][2] + mValues[3][0] * pMat[0][3],
				mValues[0][0] * pMat[1][0] + mValues[1][0] * pMat[1][1] + mValues[2][0] * pMat[1][2] + mValues[3][0] * pMat[1][3],
				mValues[0][0] * pMat[2][0] + mValues[1][0] * pMat[2][1] + mValues[2][0] * pMat[2][2] + mValues[3][0] * pMat[2][3],
				mValues[0][0] * pMat[3][0] + mValues[1][0] * pMat[3][1] + mValues[2][0] * pMat[3][2] + mValues[3][0] * pMat[3][3],

				// Row 1
				mValues[0][1] * pMat[0][0] + mValues[1][1] * pMat[0][1] + mValues[2][1] * pMat[0][2] + mValues[3][1] * pMat[0][3],
				mValues[0][1] * pMat[1][0] + mValues[1][1] * pMat[1][1] + mValues[2][1] * pMat[1][2] + mValues[3][1] * pMat[1][3],
				mValues[0][1] * pMat[2][0] + mValues[1][1] * pMat[2][1] + mValues[2][1] * pMat[2][2] + mValues[3][1] * pMat[2][3],
				mValues[0][1] * pMat[3][0] + mValues[1][1] * pMat[3][1] + mValues[2][1] * pMat[3][2] + mValues[3][1] * pMat[3][3],

				// Row 2
				mValues[0][2] * pMat[0][0] + mValues[1][2] * pMat[0][1] + mValues[2][2] * pMat[0][2] + mValues[3][2] * pMat[0][3],
				mValues[0][2] * pMat[1][0] + mValues[1][2] * pMat[1][1] + mValues[2][2] * pMat[1][2] + mValues[3][2] * pMat[1][3],
				mValues[0][2] * pMat[2][0] + mValues[1][2] * pMat[2][1] + mValues[2][2] * pMat[2][2] + mValues[3][2] * pMat[2][3],
				mValues[0][2] * pMat[3][0] + mValues[1][2] * pMat[3][1] + mValues[2][2] * pMat[3][2] + mValues[3][2] * pMat[3][3],

				// Row 3
				mValues[0][3] * pMat[0][0] + mValues[1][3] * pMat[0][1] + mValues[2][3] * pMat[0][2] + mValues[3][3] * pMat[0][3],
				mValues[0][3] * pMat[1][0] + mValues[1][3] * pMat[1][1] + mValues[2][3] * pMat[1][2] + mValues[3][3] * pMat[1][3],
				mValues[0][3] * pMat[2][0] + mValues[1][3] * pMat[2][1] + mValues[2][3] * pMat[2][2] + mValues[3][3] * pMat[2][3],
				mValues[0][3] * pMat[3][0] + mValues[1][3] * pMat[3][1] + mValues[2][3] * pMat[3][2] + mValues[3][3] * pMat[3][3]
			);
		}

		inline Col operator*(const Col &pVec) const
		{
			return Col(
				mValues[0][0] * pVec.x + mValues[1][0] * pVec.y + mValues[2][0] * pVec.z + mValues[3][0] * pVec.w,
				mValues[0][1] * pVec.x + mValues[1][1] * pVec.y + mValues[2][1] * pVec.z + mValues[3][1] * pVec.w,
				mValues[0][2] * pVec.x + mValues[1][2] * pVec.y + mValues[2][2] * pVec.z + mValues[3][2] * pVec.w,
				mValues[0][3] * pVec.x + mValues[1][3] * pVec.y + mValues[2][3] * pVec.z + mValues[3][3] * pVec.w
			);
		}

		inline Mat operator*(const T pVal) const
		{
			return Mat(
				mValues[0][0] * pVal, mValues[1][0] * pVal, mValues[2][0] * pVal, mValues[3][0] * pVal,
				mValues[0][1] * pVal, mValues[1][1] * pVal, mValues[2][1] * pVal, mValues[3][1] * pVal,
				mValues[0][2] * pVal, mValues[1][2] * pVal, mValues[2][2] * pVal, mValues[3][2] * pVal,
				mValues[0][3] * pVal, mValues[1][3] * pVal, mValues[2][3] * pVal, mValues[3][3] * pVal
			);
		}

		inline friend Mat operator*(const T pVal, const Mat &pMat)
		{
			return Mat(
				pVal * pMat[0][0], pVal * pMat[1][0], pVal * pMat[2][0], pVal * pMat[3][0],
				pVal * pMat[0][1], pVal * pMat[1][1], pVal * pMat[2][1], pVal * pMat[3][1],
				pVal * pMat[0][2], pVal * pMat[1][2], pVal * pMat[2][2], pVal * pMat[3][2],
				pVal * pMat[0][3], pVal * pMat[1][3], pVal * pMat[2][3], pVal * pMat[3][3]
			);
		}

		Matrix<1, 4, T, P> operator*(const Matrix<1, 4, T, P> &pMat) const
		{
			Matrix<1, 4, T, P> lRes;
			lRes[0] = Col(
				mValues[0][0] * pMat[0][0] + mValues[1][0] * pMat[0][1] + mValues[2][0] * pMat[0][2] + mValues[3][0] * pMat[0][3],
				mValues[0][1] * pMat[0][0] + mValues[1][1] * pMat[0][1] + mValues[2][1] * pMat[0][2] + mValues[3][1] * pMat[0][3],
				mValues[0][2] * pMat[0][0] + mValues[1][2] * pMat[0][1] + mValues[2][2] * pMat[0][2] + mValues[3][2] * pMat[0][3],
				mValues[0][3] * pMat[0][0] + mValues[1][3] * pMat[0][1] + mValues[2][3] * pMat[0][2] + mValues[3][3] * pMat[0][3]
			);
			return lRes;
		}

		Matrix<2, 4, T, P> operator*(const Matrix<2, 4, T, P> &pMat) const
		{
			Matrix<2, 4, T, P> lRes;
			lRes[0] = Col(
				mValues[0][0] * pMat[0][0] + mValues[1][0] * pMat[0][1] + mValues[2][0] * pMat[0][2] + mValues[3][0] * pMat[0][3],
				mValues[0][1] * pMat[0][0] + mValues[1][1] * pMat[0][1] + mValues[2][1] * pMat[0][2] + mValues[3][1] * pMat[0][3],
				mValues[0][2] * pMat[0][0] + mValues[1][2] * pMat[0][1] + mValues[2][2] * pMat[0][2] + mValues[3][2] * pMat[0][3],
				mValues[0][3] * pMat[0][0] + mValues[1][3] * pMat[0][1] + mValues[2][3] * pMat[0][2] + mValues[3][3] * pMat[0][3]
			);
			lRes[1] = Col(
				mValues[0][0] * pMat[1][0] + mValues[1][0] * pMat[1][1] + mValues[2][0] * pMat[1][2] + mValues[3][0] * pMat[1][3],
				mValues[0][1] * pMat[1][0] + mValues[1][1] * pMat[1][1] + mValues[2][1] * pMat[1][2] + mValues[3][1] * pMat[1][3],
				mValues[0][2] * pMat[1][0] + mValues[1][2] * pMat[1][1] + mValues[2][2] * pMat[1][2] + mValues[3][2] * pMat[1][3],
				mValues[0][3] * pMat[1][0] + mValues[1][3] * pMat[1][1] + mValues[2][3] * pMat[1][2] + mValues[3][3] * pMat[1][3]
			);
			return lRes;
		}

		Matrix<3, 4, T, P> operator*(const Matrix<3, 4, T, P> &pMat) const
		{
			Matrix<3, 4, T, P> lRes;
			lRes[0] = Col(
				mValues[0][0] * pMat[0][0] + mValues[1][0] * pMat[0][1] + mValues[2][0] * pMat[0][2] + mValues[3][0] * pMat[0][3],
				mValues[0][1] * pMat[0][0] + mValues[1][1] * pMat[0][1] + mValues[2][1] * pMat[0][2] + mValues[3][1] * pMat[0][3],
				mValues[0][2] * pMat[0][0] + mValues[1][2] * pMat[0][1] + mValues[2][2] * pMat[0][2] + mValues[3][2] * pMat[0][3],
				mValues[0][3] * pMat[0][0] + mValues[1][3] * pMat[0][1] + mValues[2][3] * pMat[0][2] + mValues[3][3] * pMat[0][3]
			);
			lRes[1] = Col(
				mValues[0][0] * pMat[1][0] + mValues[1][0] * pMat[1][1] + mValues[2][0] * pMat[1][2] + mValues[3][0] * pMat[1][3],
				mValues[0][1] * pMat[1][0] + mValues[1][1] * pMat[1][1] + mValues[2][1] * pMat[1][2] + mValues[3][1] * pMat[1][3],
				mValues[0][2] * pMat[1][0] + mValues[1][2] * pMat[1][1] + mValues[2][2] * pMat[1][2] + mValues[3][2] * pMat[1][3],
				mValues[0][3] * pMat[1][0] + mValues[1][3] * pMat[1][1] + mValues[2][3] * pMat[1][2] + mValues[3][3] * pMat[1][3]
			);
			lRes[2] = Col(
				mValues[0][0] * pMat[2][0] + mValues[1][0] * pMat[2][1] + mValues[2][0] * pMat[2][2] + mValues[3][0] * pMat[2][3],
				mValues[0][1] * pMat[2][0] + mValues[1][1] * pMat[2][1] + mValues[2][1] * pMat[2][2] + mValues[3][1] * pMat[2][3],
				mValues[0][2] * pMat[2][0] + mValues[1][2] * pMat[2][1] + mValues[2][2] * pMat[2][2] + mValues[3][2] * pMat[2][3],
				mValues[0][3] * pMat[2][0] + mValues[1][3] * pMat[2][1] + mValues[2][3] * pMat[2][2] + mValues[3][3] * pMat[2][3]
			);
			return lRes;
		}

		template<unsigned int S>
		Matrix<S, 4, T, P> operator*(const Matrix<S, 4, T, P> &pMat) const
		{
			Matrix<S, 4, T, P> lRes;

			for (int i = 0; i < S; i++)
			{
				lRes[i] = Col(
					mValues[0][0] * pMat[i][0] + mValues[1][0] * pMat[i][1] + mValues[2][0] * pMat[i][2] + mValues[3][0] * pMat[i][3],
					mValues[0][1] * pMat[i][0] + mValues[1][1] * pMat[i][1] + mValues[2][1] * pMat[i][2] + mValues[3][1] * pMat[i][3],
					mValues[0][2] * pMat[i][0] + mValues[1][2] * pMat[i][1] + mValues[2][2] * pMat[i][2] + mValues[3][2] * pMat[i][3],
					mValues[0][3] * pMat[i][0] + mValues[1][3] * pMat[i][1] + mValues[2][3] * pMat[i][2] + mValues[3][3] * pMat[i][3]
				);
			}

			return lRes;
		}

		inline Mat operator/(const Mat &pMat) const
		{
			return *this * pMat.Inverted();
		}

		inline Mat operator/(const T pVal) const
		{
			const T lInv = T(1) / pVal;
			return Mat(
				mValues[0][0] * lInv, mValues[1][0] * lInv, mValues[2][0] * lInv, mValues[3][0] * lInv,
				mValues[0][1] * lInv, mValues[1][1] * lInv, mValues[2][1] * lInv, mValues[3][1] * lInv,
				mValues[0][2] * lInv, mValues[1][2] * lInv, mValues[2][2] * lInv, mValues[3][2] * lInv,
				mValues[0][3] * lInv, mValues[1][3] * lInv, mValues[2][3] * lInv, mValues[3][3] * lInv
			);
		}

		// -- Boolean operators --

		inline bool operator==(const Mat &pMat) const
		{
			return mValues[0] == pMat[0] && mValues[1] == pMat[1] && mValues[2] == pMat[2] && mValues[3] == pMat[3];
		}

		inline bool operator!=(const Mat &pMat) const
		{
			return mValues[0] != pMat[0] || mValues[1] != pMat[1] || mValues[2] != pMat[2] || mValues[3] != pMat[3];
		}

		// -- Convertion operators --

		template<typename U, typename Q>
		inline operator Matrix<4, 4, U, Q>() const
		{
			return Matrix<4, 4, U, Q>(
				(U)mValues[0][0], (U)mValues[1][0], (U)mValues[2][0], (U)mValues[3][0],
				(U)mValues[0][1], (U)mValues[1][1], (U)mValues[2][1], (U)mValues[3][1],
				(U)mValues[0][2], (U)mValues[1][2], (U)mValues[2][2], (U)mValues[3][2],
				(U)mValues[0][3], (U)mValues[1][3], (U)mValues[2][3], (U)mValues[3][3]
			);
		}

		// -- Stream operators --

		inline friend std::ostream &operator<<(std::ostream &pOStream, const Mat &pMat)
		{
			return pOStream	<< pMat[0][0] << ',' << pMat[1][0] << ',' << pMat[2][0] << ',' << pMat[3][0] << '\n'
							<< pMat[0][1] << ',' << pMat[1][1] << ',' << pMat[2][1] << ',' << pMat[3][1] << '\n'
							<< pMat[0][2] << ',' << pMat[1][2] << ',' << pMat[2][2] << ',' << pMat[3][2] << '\n'
							<< pMat[0][3] << ',' << pMat[1][3] << ',' << pMat[2][3] << ',' << pMat[3][3] << '\n';
		}

		// -- Getters --

		Mat GetAdjugate() const
		{
			Mat lRes;

			// Col 1
			lRes[0][0] = mValues[0][0] * (
				  mValues[1][1] * (mValues[2][2] * mValues[3][3] - mValues[3][2] * mValues[2][3])
				- mValues[2][1] * (mValues[1][2] * mValues[3][3] - mValues[3][2] * mValues[1][3])
				+ mValues[3][1] * (mValues[1][2] * mValues[2][3] - mValues[2][2] * mValues[1][3])
			);
			lRes[0][1] = -mValues[1][0] * (
				  mValues[0][1] * (mValues[2][2] * mValues[3][3] - mValues[3][2] * mValues[2][3])
				- mValues[2][1] * (mValues[0][2] * mValues[3][3] - mValues[3][2] * mValues[0][3])
				+ mValues[3][1] * (mValues[0][2] * mValues[2][3] - mValues[2][2] * mValues[0][3])
			);
			lRes[0][2] = mValues[2][0] * (
				  mValues[0][1] * (mValues[1][2] * mValues[3][3] - mValues[3][2] * mValues[1][3])
				- mValues[1][1] * (mValues[0][2] * mValues[3][3] - mValues[3][2] * mValues[0][3])
				+ mValues[3][1] * (mValues[0][2] * mValues[1][3] - mValues[1][2] * mValues[0][3])
			);
			lRes[0][3] = -mValues[3][0] * (
				  mValues[0][1] * (mValues[1][2] * mValues[2][3] - mValues[2][2] * mValues[1][3])
				- mValues[1][1] * (mValues[0][2] * mValues[2][3] - mValues[2][2] * mValues[0][3])
				+ mValues[2][1] * (mValues[0][2] * mValues[1][3] - mValues[1][2] * mValues[0][3])
			);

			// Col 2
			lRes[1][0] = -mValues[0][1] * (
				  mValues[1][0] * (mValues[2][2] * mValues[3][3] - mValues[3][2] * mValues[2][3])
				- mValues[2][0] * (mValues[1][2] * mValues[3][3] - mValues[3][2] * mValues[1][3])
				+ mValues[3][0] * (mValues[1][2] * mValues[2][3] - mValues[2][2] * mValues[1][3])
			);
			lRes[1][1] = mValues[1][1] * (
				  mValues[0][0] * (mValues[2][2] * mValues[3][3] - mValues[3][2] * mValues[2][3])
				- mValues[2][0] * (mValues[0][2] * mValues[3][3] - mValues[3][2] * mValues[0][3])
				+ mValues[3][0] * (mValues[0][2] * mValues[2][3] - mValues[2][2] * mValues[0][3])
			);
			lRes[1][2] = -mValues[2][1] * (
				  mValues[0][0] * (mValues[1][2] * mValues[3][3] - mValues[3][2] * mValues[1][3])
				- mValues[1][0] * (mValues[0][2] * mValues[3][3] - mValues[3][2] * mValues[0][3])
				+ mValues[3][0] * (mValues[0][2] * mValues[1][3] - mValues[1][2] * mValues[0][3])
			);
			lRes[1][3] = mValues[3][1] * (
				  mValues[0][0] * (mValues[1][2] * mValues[2][3] - mValues[2][2] * mValues[1][3])
				- mValues[1][0] * (mValues[0][2] * mValues[2][3] - mValues[2][2] * mValues[0][3])
				+ mValues[2][0] * (mValues[0][2] * mValues[1][3] - mValues[1][2] * mValues[0][3])
			);

			// Col 3
			lRes[2][0] = mValues[0][2] * (
				  mValues[1][0] * (mValues[2][1] * mValues[3][3] - mValues[3][1] * mValues[2][3])
				- mValues[2][0] * (mValues[1][1] * mValues[3][3] - mValues[3][1] * mValues[1][3])
				+ mValues[3][0] * (mValues[1][1] * mValues[2][3] - mValues[2][1] * mValues[1][3])
			);
			lRes[2][1] = -mValues[1][2] * (
				  mValues[0][0] * (mValues[2][1] * mValues[3][3] - mValues[3][1] * mValues[2][3])
				- mValues[2][0] * (mValues[0][1] * mValues[3][3] - mValues[3][1] * mValues[0][3])
				+ mValues[3][0] * (mValues[0][1] * mValues[2][3] - mValues[2][1] * mValues[0][3])
			);
			lRes[2][2] = mValues[2][2] * (
				  mValues[0][0] * (mValues[1][1] * mValues[3][3] - mValues[3][1] * mValues[1][3])
				- mValues[1][0] * (mValues[0][1] * mValues[3][3] - mValues[3][1] * mValues[0][3])
				+ mValues[3][0] * (mValues[0][1] * mValues[1][3] - mValues[1][1] * mValues[0][3])
			);
			lRes[2][3] = -mValues[3][2] * (
				  mValues[0][0] * (mValues[2][1] * mValues[3][3] - mValues[3][1] * mValues[2][3])
				- mValues[1][0] * (mValues[1][1] * mValues[3][3] - mValues[3][1] * mValues[1][3])
				+ mValues[2][0] * (mValues[1][1] * mValues[2][3] - mValues[2][1] * mValues[1][3])
			);

			// Col 4
			lRes[3][0] = -mValues[0][3] * (
				  mValues[1][0] * (mValues[2][1] * mValues[3][2] - mValues[3][1] * mValues[2][2])
				- mValues[2][0] * (mValues[1][1] * mValues[3][2] - mValues[3][1] * mValues[1][2])
				+ mValues[3][0] * (mValues[1][1] * mValues[2][2] - mValues[2][1] * mValues[1][2])
			);
			lRes[3][1] = mValues[1][3] * (
				  mValues[0][0] * (mValues[2][1] * mValues[3][2] - mValues[3][1] * mValues[2][2])
				- mValues[2][0] * (mValues[0][1] * mValues[3][2] - mValues[3][1] * mValues[0][2])
				+ mValues[3][0] * (mValues[0][1] * mValues[2][2] - mValues[2][1] * mValues[0][2])
			);
			lRes[3][2] = -mValues[2][3] * (
				  mValues[0][0] * (mValues[1][1] * mValues[3][2] - mValues[3][1] * mValues[1][2])
				- mValues[1][0] * (mValues[0][1] * mValues[3][2] - mValues[3][1] * mValues[0][2])
				+ mValues[3][0] * (mValues[0][1] * mValues[1][2] - mValues[1][1] * mValues[0][2])
			);
			lRes[3][3] = mValues[3][3] * (
				  mValues[0][0] * (mValues[1][1] * mValues[2][2] - mValues[2][1] * mValues[1][2])
				- mValues[1][0] * (mValues[0][1] * mValues[2][2] - mValues[2][1] * mValues[0][2])
				+ mValues[2][0] * (mValues[0][1] * mValues[1][2] - mValues[1][1] * mValues[0][2])
			);

			return lRes;
		}

		Mat GetCofactor() const
		{
			Mat lRes;

			// Col 1
			lRes[0][0] = mValues[0][0] * (
				  mValues[1][1] * (mValues[2][2] * mValues[3][3] - mValues[3][2] * mValues[2][3])
				- mValues[2][1] * (mValues[1][2] * mValues[3][3] - mValues[3][2] * mValues[1][3])
				+ mValues[3][1] * (mValues[1][2] * mValues[2][3] - mValues[2][2] * mValues[1][3])
			);
			lRes[0][1] = -mValues[0][1] * (
				  mValues[1][0] * (mValues[2][2] * mValues[3][3] - mValues[3][2] * mValues[2][3])
				- mValues[2][0] * (mValues[1][2] * mValues[3][3] - mValues[3][2] * mValues[1][3])
				+ mValues[3][0] * (mValues[1][2] * mValues[2][3] - mValues[2][2] * mValues[1][3])
			);
			lRes[0][2] = mValues[0][2] * (
				  mValues[1][0] * (mValues[2][1] * mValues[3][3] - mValues[3][1] * mValues[2][3])
				- mValues[2][0] * (mValues[1][1] * mValues[3][3] - mValues[3][1] * mValues[1][3])
				+ mValues[3][0] * (mValues[1][1] * mValues[2][3] - mValues[2][1] * mValues[1][3])
			);
			lRes[0][3] = -mValues[0][3] * (
				  mValues[1][0] * (mValues[2][1] * mValues[3][2] - mValues[3][1] * mValues[2][2])
				- mValues[2][0] * (mValues[1][1] * mValues[3][2] - mValues[3][1] * mValues[1][2])
				+ mValues[3][0] * (mValues[1][1] * mValues[2][2] - mValues[2][1] * mValues[1][2])
			);

			// Col 2
			lRes[1][0] = -mValues[1][0] * (
				  mValues[0][1] * (mValues[2][2] * mValues[3][3] - mValues[3][2] * mValues[2][3])
				- mValues[2][1] * (mValues[0][2] * mValues[3][3] - mValues[3][2] * mValues[0][3])
				+ mValues[3][1] * (mValues[0][2] * mValues[2][3] - mValues[2][2] * mValues[0][3])
			);
			lRes[1][1] = mValues[1][1] * (
				  mValues[0][0] * (mValues[2][2] * mValues[3][3] - mValues[3][2] * mValues[2][3])
				- mValues[2][0] * (mValues[0][2] * mValues[3][3] - mValues[3][2] * mValues[0][3])
				+ mValues[3][0] * (mValues[0][2] * mValues[2][3] - mValues[2][2] * mValues[0][3])
			);
			lRes[1][2] = -mValues[1][2] * (
				  mValues[0][0] * (mValues[2][1] * mValues[3][3] - mValues[3][1] * mValues[2][3])
				- mValues[2][0] * (mValues[0][1] * mValues[3][3] - mValues[3][1] * mValues[0][3])
				+ mValues[3][0] * (mValues[0][1] * mValues[2][3] - mValues[2][1] * mValues[0][3])
			);
			lRes[1][3] = mValues[1][3] * (
				  mValues[0][0] * (mValues[2][1] * mValues[3][2] - mValues[3][1] * mValues[2][2])
				- mValues[2][0] * (mValues[0][1] * mValues[3][2] - mValues[3][1] * mValues[0][2])
				+ mValues[3][0] * (mValues[0][1] * mValues[2][2] - mValues[2][1] * mValues[0][2])
			);

			// Col 3
			lRes[2][0] = mValues[2][0] * (
				  mValues[0][1] * (mValues[1][2] * mValues[3][3] - mValues[3][2] * mValues[1][3])
				- mValues[1][1] * (mValues[0][2] * mValues[3][3] - mValues[3][2] * mValues[0][3])
				+ mValues[3][1] * (mValues[0][2] * mValues[1][3] - mValues[1][2] * mValues[0][3])
			);
			lRes[2][1] = -mValues[2][1] * (
				  mValues[0][0] * (mValues[1][2] * mValues[3][3] - mValues[3][2] * mValues[1][3])
				- mValues[1][0] * (mValues[0][2] * mValues[3][3] - mValues[3][2] * mValues[0][3])
				+ mValues[3][0] * (mValues[0][2] * mValues[1][3] - mValues[1][2] * mValues[0][3])
			);
			lRes[2][2] = mValues[2][2] * (
				  mValues[0][0] * (mValues[1][1] * mValues[3][3] - mValues[3][1] * mValues[1][3])
				- mValues[1][0] * (mValues[0][1] * mValues[3][3] - mValues[3][1] * mValues[0][3])
				+ mValues[3][0] * (mValues[0][1] * mValues[1][3] - mValues[1][1] * mValues[0][3])
			);
			lRes[2][3] = -mValues[2][3] * (
				  mValues[0][0] * (mValues[1][1] * mValues[3][2] - mValues[3][1] * mValues[1][2])
				- mValues[1][0] * (mValues[0][1] * mValues[3][2] - mValues[3][1] * mValues[0][2])
				+ mValues[3][0] * (mValues[0][1] * mValues[1][2] - mValues[1][1] * mValues[0][2])
			);

			// Col 4
			lRes[3][0] = -mValues[3][0] * (
				  mValues[0][1] * (mValues[1][2] * mValues[2][3] - mValues[2][2] * mValues[1][3])
				- mValues[1][1] * (mValues[0][2] * mValues[2][3] - mValues[2][2] * mValues[0][3])
				+ mValues[2][1] * (mValues[0][2] * mValues[1][3] - mValues[1][2] * mValues[0][3])
			);
			lRes[3][1] = mValues[3][1] * (
				  mValues[0][0] * (mValues[1][2] * mValues[2][3] - mValues[2][2] * mValues[1][3])
				- mValues[1][0] * (mValues[0][2] * mValues[2][3] - mValues[2][2] * mValues[0][3])
				+ mValues[2][0] * (mValues[0][2] * mValues[1][3] - mValues[1][2] * mValues[0][3])
			);
			lRes[3][2] = -mValues[3][2] * (
				  mValues[0][0] * (mValues[2][1] * mValues[3][3] - mValues[3][1] * mValues[2][3])
				- mValues[1][0] * (mValues[1][1] * mValues[3][3] - mValues[3][1] * mValues[1][3])
				+ mValues[2][0] * (mValues[1][1] * mValues[2][3] - mValues[2][1] * mValues[1][3])
			);
			lRes[3][3] = mValues[3][3] * (
				  mValues[0][0] * (mValues[1][1] * mValues[2][2] - mValues[2][1] * mValues[1][2])
				- mValues[1][0] * (mValues[0][1] * mValues[2][2] - mValues[2][1] * mValues[0][2])
				+ mValues[2][0] * (mValues[0][1] * mValues[1][2] - mValues[1][1] * mValues[0][2])
			);

			return lRes;
		}

		T GetDeterminant() const
		{
			return mValues[0][0] * (
				  mValues[1][1] * (mValues[2][2] * mValues[3][3] - mValues[3][2] * mValues[2][3])
				- mValues[2][1] * (mValues[1][2] * mValues[3][3] - mValues[3][2] * mValues[1][3])
				+ mValues[3][1] * (mValues[1][2] * mValues[2][3] - mValues[2][2] * mValues[1][3])
			)
			- mValues[1][0] * (
				  mValues[0][1] * (mValues[2][2] * mValues[3][3] - mValues[3][2] * mValues[2][3])
				- mValues[2][1] * (mValues[0][2] * mValues[3][3] - mValues[3][2] * mValues[0][3])
				+ mValues[3][1] * (mValues[0][2] * mValues[2][3] - mValues[2][2] * mValues[0][3])
			)
			+ mValues[2][0] * (
				  mValues[0][1] * (mValues[1][2] * mValues[3][3] - mValues[3][2] * mValues[1][3])
				- mValues[1][1] * (mValues[0][2] * mValues[3][3] - mValues[3][2] * mValues[0][3])
				+ mValues[3][1] * (mValues[0][2] * mValues[1][3] - mValues[1][2] * mValues[0][3])
			)
			- mValues[3][0] * (
				  mValues[0][1] * (mValues[1][2] * mValues[2][3] - mValues[2][2] * mValues[1][3])
				- mValues[1][1] * (mValues[0][2] * mValues[2][3] - mValues[2][2] * mValues[0][3])
				+ mValues[2][1] * (mValues[0][2] * mValues[1][3] - mValues[1][2] * mValues[0][3])
			);
		}

		Mat Inverted() const
		{
			return GetAdjugate() * (T(1) / GetDeterminant());
		}

		Mat Transposed() const
		{
			return Mat(
				mValues[0][0], mValues[0][1], mValues[0][2], mValues[0][3],
				mValues[1][0], mValues[1][1], mValues[1][2], mValues[1][3],
				mValues[2][0], mValues[2][1], mValues[2][2], mValues[2][3],
				mValues[3][0], mValues[3][1], mValues[3][2], mValues[3][3]
			);
		}

		// -- Transformations --

		Mat &Invert()
		{
			*this = GetAdjugate() * (T(1) / GetDeterminant());
			return *this;
		}

		Mat &Transpose()
		{
			mValues[0][1] += mValues[1][0];
			mValues[1][0] = mValues[0][1] - mValues[1][0];
			mValues[0][1] -= mValues[1][0];

			mValues[0][2] += mValues[2][0];
			mValues[2][0] = mValues[0][2] - mValues[2][0];
			mValues[0][2] -= mValues[2][0];

			mValues[0][3] += mValues[3][0];
			mValues[3][0] = mValues[0][3] - mValues[3][0];
			mValues[0][3] -= mValues[3][0];

			mValues[2][1] += mValues[1][2];
			mValues[1][2] = mValues[2][1] - mValues[1][2];
			mValues[2][1] -= mValues[1][2];

			mValues[3][1] += mValues[1][3];
			mValues[1][3] = mValues[3][1] - mValues[1][3];
			mValues[3][1] -= mValues[1][3];

			mValues[3][2] += mValues[2][3];
			mValues[2][3] = mValues[3][2] - mValues[2][3];
			mValues[3][2] -= mValues[2][3];
			return *this;
		}

		// -- Static getters --

		inline static Mat Identity()
		{
			return Mat(T(1));
		}
	};

	typedef Matrix<4, 4, float, float> Matrix4;
	typedef Matrix<4, 4, float, float> Matrix4x4;
}

