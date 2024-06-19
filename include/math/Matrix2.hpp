#pragma once

#include "Vector2.hpp"
#include "Matrix.hpp"

namespace Math
{
	template<typename T, typename P>
	struct Matrix<2, 2, T, P>
	{
	private:
		typedef Matrix<2, 2, T, P> Mat;
		typedef Vector<2, T, P> Col;
		Col mValues[2];

	public:
		// -- Constructors --

		Matrix()
			: mValues{ Col(), Col() }
		{

		}

		explicit Matrix(const T pVal)
			: mValues{
				Col(pVal, T(0)),
				Col(T(0), pVal)
			}
		{

		}

		Matrix(
			const T pX0, const T pX1,
			const T pX2, const T pX3
		)
			: mValues{
				Col(pX0, pX2),
				Col(pX1, pX3)
			}
		{

		}

		Matrix(const T pVals[4])
			: mValues{
				Col(pVals[0], pVals[2]),
				Col(pVals[1], pVals[3])
			}
		{

		}

		Matrix(const Col &pX0, const Col &pX1)
			: mValues{
				Col(pX0[0], pX1[0]),
				Col(pX0[1], pX1[1])
			}
		{

		}

		Matrix(const Col pVals[2])
			: mValues{
				Col(pVals[0][0], pVals[1][0]),
				Col(pVals[0][1], pVals[1][1])
			}
		{

		}

		Matrix(const Mat &pMat)
			: mValues{ pMat[0], pMat[1] }
		{

		}

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

		Mat &operator+=(const Mat &pMat)
		{
			mValues[0] += pMat[0];
			mValues[1] += pMat[1];
			return *this;
		}

		Mat &operator-=(const Mat &pMat)
		{
			mValues[0] -= pMat[0];
			mValues[1] -= pMat[1];
			return *this;
		}

		Mat &operator*=(const Mat &pMat)
		{
			const Mat lTemp(*this);
			mValues[0][0] = lTemp[0][0] * pMat[0][0] + lTemp[1][0] * pMat[0][1];
			mValues[0][1] = lTemp[0][1] * pMat[0][0] + lTemp[1][1] * pMat[0][1];
			mValues[1][0] = lTemp[0][0] * pMat[1][0] + lTemp[1][0] * pMat[1][1];
			mValues[1][1] = lTemp[0][1] * pMat[1][0] + lTemp[1][1] * pMat[1][1];
			return *this;
		}

		Mat &operator*=(const T pVal)
		{
			mValues[0] *= pVal;
			mValues[1] *= pVal;
			return *this;
		}

		Mat &operator/=(const Mat &pMat)
		{
			return *this *= pMat.Inverted();
		}

		Mat &operator/=(const T pVal)
		{
			const T lInv = T(1) / pVal;
			mValues[0] *= lInv;
			mValues[1] *= lInv;
			return *this;
		}

		// -- Unary operators --

		Mat operator+() const
		{
			return *this;
		}

		Mat operator-() const
		{
			return Mat(
				-mValues[0][0], -mValues[1][0],
				-mValues[0][1], -mValues[1][1]
			);
		}

		// -- Binary operators --

		Mat operator+(const Mat &pMat) const
		{
			return Mat(
				mValues[0][0] + pMat[0][0], mValues[1][0] + pMat[1][0],
				mValues[0][1] + pMat[0][1], mValues[1][1] + pMat[1][1]
			);
		}

		Mat operator-(const Mat &pMat) const
		{
			return Mat(
				mValues[0][0] - pMat[0][0], mValues[1][0] - pMat[1][0],
				mValues[0][1] - pMat[0][1], mValues[1][1] - pMat[1][1]
			);
		}

		Mat operator*(const Mat &pMat) const
		{
			return Mat(
				mValues[0][0] * pMat[0][0] + mValues[1][0] * pMat[0][1],
				mValues[0][0] * pMat[1][0] + mValues[1][0] * pMat[1][1],
				mValues[0][1] * pMat[0][0] + mValues[1][1] * pMat[0][1],
				mValues[0][1] * pMat[1][0] + mValues[1][1] * pMat[1][1]
			);
		}

		Col operator*(const Col &pVec) const
		{
			return Col(
				mValues[0][0] * pVec[0] + mValues[1][0] * pVec[1],
				mValues[0][1] * pVec[0] + mValues[1][1] * pVec[1]
			);
		}

		Mat operator*(const T pVal) const
		{
			return Mat(
				mValues[0][0] * pVal, mValues[1][0] * pVal,
				mValues[0][1] * pVal, mValues[1][1] * pVal
			);
		}

		Matrix<1, 2, T, P> operator*(const Matrix<1, 2, T, P> &pMat) const
		{
			Matrix<1, 2, T, P> lRes;
			lRes[0] = Col(
				mValues[0][0] * pMat[0][0] + mValues[1][0] * pMat[0][1],
				mValues[0][1] * pMat[0][0] + mValues[1][1] * pMat[0][1]
			);
			return lRes;
		}

		Matrix<3, 2, T, P> operator*(const Matrix<3, 2, T, P> &pMat) const
		{
			Matrix<3, 2, T, P> lRes;
			lRes[0] = Col(
				mValues[0][0] * pMat[0][0] + mValues[1][0] * pMat[0][1],
				mValues[0][1] * pMat[0][0] + mValues[1][1] * pMat[0][1]
			);
			lRes[1] = Col(
				mValues[0][0] * pMat[1][0] + mValues[1][0] * pMat[1][1],
				mValues[0][1] * pMat[1][0] + mValues[1][1] * pMat[1][1]
			);
			lRes[2] = Col(
				mValues[0][0] * pMat[2][0] + mValues[1][0] * pMat[2][1],
				mValues[0][1] * pMat[2][0] + mValues[1][1] * pMat[2][1]
			);
			return lRes;
		}

		Matrix<4, 2, T, P> operator*(const Matrix<4, 2, T, P> &pMat) const
		{
			Matrix<4, 2, T, P> lRes;
			lRes[0] = Col(
				mValues[0][0] * pMat[0][0] + mValues[1][0] * pMat[0][1],
				mValues[0][1] * pMat[0][0] + mValues[1][1] * pMat[0][1]
			);
			lRes[1] = Col(
				mValues[0][0] * pMat[1][0] + mValues[1][0] * pMat[1][1],
				mValues[0][1] * pMat[1][0] + mValues[1][1] * pMat[1][1]
			);
			lRes[2] = Col(
				mValues[0][0] * pMat[2][0] + mValues[1][0] * pMat[2][1],
				mValues[0][1] * pMat[2][0] + mValues[1][1] * pMat[2][1]
			);
			lRes[3] = Col(
				mValues[0][0] * pMat[3][0] + mValues[1][0] * pMat[3][1],
				mValues[0][1] * pMat[3][0] + mValues[1][1] * pMat[3][1]
			);
			return lRes;
		}

		template<unsigned int S>
		Matrix<S, 2, T, P> operator*(const Matrix<S, 2, T, P> &pMat) const
		{
			Matrix<S, 2, T, P> lRes;

			for (int i = 0; i < S; i++)
			{
				lRes[i] = Col(
					mValues[0][0] * pMat[i][0] + mValues[1][0] * pMat[i][1],
					mValues[0][1] * pMat[i][0] + mValues[1][1] * pMat[i][1]
				);
			}

			return lRes;
		}

		Mat operator/(const Mat &pMat) const
		{
			return *this * pMat.Inverted();
		}

		Mat operator/(const T pVal) const
		{
			const T lInv = T(1) / pVal;
			return Mat(
				mValues[0][0] * lInv, mValues[1][0] * lInv,
				mValues[0][1] * lInv, mValues[1][1] * lInv
			);
		}

		// -- Boolean operators --

		bool operator==(const Mat &pMat) const
		{
			return mValues[0] == pMat[0] && mValues[1] == pMat[1];
		}

		bool operator!=(const Mat &mat) const
		{
			return mValues[0] != mat[0] || mValues[1] != mat[1];
		}

		// -- Convertion operators --

		template<typename U, typename Q>
		operator Matrix<2, 2, U, Q>() const
		{
			return Matrix<2, 2, U, Q>(
				(U)mValues[0][0], (U)mValues[1][0],
				(U)mValues[0][1], (U)mValues[1][1]
			);
		}

		// -- Stream operators --

		friend std::ostream &operator<<(std::ostream &lOStream, const Mat &pMat)
		{
			return lOStream	<< pMat[0][0] << ',' << pMat[1][0] << '\n'
							<< pMat[0][1] << ',' << pMat[1][1] << '\n';
		}

		// -- Getters --

		Mat GetAdjugate() const
		{
			return Mat(
				mValues[1][1], -mValues[1][0],
				-mValues[0][1], mValues[0][0]
			);
		}

		Mat GetCofactor() const
		{
			return Mat(
				mValues[1][1], -mValues[0][1],
				-mValues[1][0], mValues[0][0]
			);
		}

		T GetDeterminant() const
		{
			return mValues[0][0] * mValues[1][1] - mValues[1][0] * mValues[0][1];
		}

		Mat Inverted() const
		{
			return GetAdjugate() * (T(1) / GetDeterminant());
		}

		Mat Transposed() const
		{
			return Mat(
				mValues[0][0], mValues[0][1],
				mValues[1][0], mValues[1][1]
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
			return *this;
		}

		// -- Static getters --

		static inline Mat Identity()
		{
			return Mat(T(1));
		}
	};

	typedef Matrix<2, 2, float, float> Matrix2;
	typedef Matrix<2, 2, float, float> Matrix2x2;
}

