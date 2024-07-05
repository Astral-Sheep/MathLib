#pragma once

#include "Vector3.hpp"
#include "Matrix.hpp"

namespace Math
{
	template<typename T, typename P>
	struct Matrix<3, 3, T, P>
	{
	private:
		typedef Matrix<3, 3, T, P> Mat;
		typedef Vector<3, T, P> Col;
		Col mValues[3];

	public:
		// -- Constructors --

		Matrix()
			: mValues{ Col(), Col(), Col() } {}

		explicit Matrix(const T pVal)
			: mValues{
				Col(pVal, T(0), T(0)),
				Col(T(0), pVal, T(0)),
				Col(T(0), T(0), pVal)
			} {}

		Matrix(
			const T pX0, const T pX1, const T pX2,
			const T pX3, const T pX4, const T pX5,
			const T pX6, const T pX7, const T pX8
		)
			: mValues{
				Col(pX0, pX3, pX6),
				Col(pX1, pX4, pX7),
				Col(pX2, pX5, pX8)
			} {}

		Matrix(const T pVals[9])
			: mValues{
				Col(pVals[0], pVals[3], pVals[6]),
				Col(pVals[1], pVals[4], pVals[7]),
				Col(pVals[2], pVals[5], pVals[8])
			} {}

		Matrix(const Col &pX0, const Col &pX1, const Col &pX2)
			: mValues{
				Col(pX0[0], pX1[0], pX2[0]),
				Col(pX0[1], pX1[1], pX2[1]),
				Col(pX0[2], pX1[2], pX2[2])
			} {}

		Matrix(const Col pVals[3])
			: mValues{
				Col(pVals[0][0], pVals[1][0], pVals[2][0]),
				Col(pVals[0][1], pVals[1][1], pVals[2][1]),
				Col(pVals[0][2], pVals[1][2], pVals[2][2])
			} {}

		Matrix(const Mat &pMat)
			: mValues{ pMat[0], pMat[1], pMat[2] } {}

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
			return *this;
		}

		inline Mat &operator-=(const Mat &pMat)
		{
			mValues[0] -= pMat[0];
			mValues[1] -= pMat[1];
			mValues[2] -= pMat[2];
			return *this;
		}

		inline Mat &operator*=(const Mat &pMat)
		{
			const Mat lTemp(*this);
			mValues[0][0] = lTemp[0][0] * pMat[0][0] + lTemp[1][0] * pMat[0][1] + lTemp[2][0] * pMat[0][2];
			mValues[0][1] = lTemp[0][1] * pMat[0][0] + lTemp[1][1] * pMat[0][1] + lTemp[2][1] * pMat[0][2];
			mValues[0][2] = lTemp[0][2] * pMat[0][0] + lTemp[1][2] * pMat[0][1] + lTemp[2][2] * pMat[0][2];

			mValues[1][0] = lTemp[0][0] * pMat[1][0] + lTemp[1][0] * pMat[1][1] + lTemp[2][0] * pMat[1][2];
			mValues[1][1] = lTemp[0][1] * pMat[1][0] + lTemp[1][1] * pMat[1][1] + lTemp[2][1] * pMat[1][2];
			mValues[1][2] = lTemp[0][2] * pMat[1][0] + lTemp[1][2] * pMat[1][1] + lTemp[2][2] * pMat[1][2];

			mValues[2][0] = lTemp[0][0] * pMat[2][0] + lTemp[1][0] * pMat[2][1] + lTemp[2][0] * pMat[2][2];
			mValues[2][1] = lTemp[0][1] * pMat[2][0] + lTemp[1][1] * pMat[2][1] + lTemp[2][1] * pMat[2][2];
			mValues[2][2] = lTemp[0][2] * pMat[2][0] + lTemp[1][2] * pMat[2][1] + lTemp[2][2] * pMat[2][2];
			return *this;
		}

		inline Mat &operator*=(const T pVal)
		{
			mValues[0] *= pVal;
			mValues[1] *= pVal;
			mValues[2] *= pVal;
			return *this;
		}

		inline Mat &operator/=(const Mat &pMat)
		{
			return *this *= pMat.Inverted();
		}

		inline Mat &operator/=(const T pVal)
		{
			const T lInv = T(1) / pVal;
			mValues[0] *= lInv;
			mValues[1] *= lInv;
			mValues[2] *= lInv;
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
				-mValues[0][0], -mValues[1][0], -mValues[2][0],
				-mValues[0][1], -mValues[1][1], -mValues[2][1],
				-mValues[0][2], -mValues[1][2], -mValues[2][2]
			);
		}

		// -- Binary operators --

		inline Mat operator+(const Mat &pMat) const
		{
			return Mat(
				mValues[0][0] + pMat[0][0], mValues[1][0] + pMat[1][0], mValues[2][0] + pMat[2][0],
				mValues[0][1] + pMat[0][1], mValues[1][1] + pMat[1][1], mValues[2][1] + pMat[2][1],
				mValues[0][2] + pMat[0][2], mValues[1][2] + pMat[1][2], mValues[2][2] + pMat[2][2]
			);
		}

		inline Mat operator-(const Mat &pMat) const
		{
			return Mat(
				mValues[0][0] - pMat[0][0], mValues[1][0] - pMat[1][0], mValues[2][0] - pMat[2][0],
				mValues[0][1] - pMat[0][1], mValues[1][1] - pMat[1][1], mValues[2][1] - pMat[2][1],
				mValues[0][2] - pMat[0][2], mValues[1][2] - pMat[1][2], mValues[2][2] - pMat[2][2]
			);
		}

		inline Mat operator*(const Mat &pMat) const
		{
			return Mat(
				mValues[0][0] * pMat[0][0] + mValues[1][0] * pMat[0][1] + mValues[2][0] * pMat[0][2],
				mValues[0][0] * pMat[1][0] + mValues[1][0] * pMat[1][1] + mValues[2][0] * pMat[1][2],
				mValues[0][0] * pMat[2][0] + mValues[1][0] * pMat[2][1] + mValues[2][0] * pMat[2][2],

				mValues[0][1] * pMat[0][0] + mValues[1][1] * pMat[0][1] + mValues[2][1] * pMat[0][2],
				mValues[0][1] * pMat[1][0] + mValues[1][1] * pMat[1][1] + mValues[2][1] * pMat[1][2],
				mValues[0][1] * pMat[2][0] + mValues[1][1] * pMat[2][1] + mValues[2][1] * pMat[2][2],

				mValues[0][2] * pMat[0][0] + mValues[1][2] * pMat[0][1] + mValues[2][2] * pMat[0][2],
				mValues[0][2] * pMat[1][0] + mValues[1][2] * pMat[1][1] + mValues[2][2] * pMat[1][2],
				mValues[0][2] * pMat[2][0] + mValues[1][2] * pMat[2][1] + mValues[2][2] * pMat[2][2]
			);
		}

		inline Col operator*(const Col &pVec) const
		{
			return Col(
				mValues[0][0] * pVec.x + mValues[1][0] * pVec.y + mValues[2][0] * pVec.z,
				mValues[0][1] * pVec.x + mValues[1][1] * pVec.y + mValues[2][1] * pVec.z,
				mValues[0][2] * pVec.x + mValues[1][2] * pVec.y + mValues[2][2] * pVec.z
			);
		}

		inline Mat operator*(const T pVal) const
		{
			return Mat(
				mValues[0][0] * pVal, mValues[1][0] * pVal, mValues[2][0] * pVal,
				mValues[0][1] * pVal, mValues[1][1] * pVal, mValues[2][1] * pVal,
				mValues[0][2] * pVal, mValues[1][2] * pVal, mValues[2][2] * pVal
			);
		}

		inline friend Mat operator*(const T pVal, const Mat &pMat)
		{
			return Mat(
				pVal * pMat[0][0], pVal * pMat[1][0], pVal * pMat[2][0],
				pVal * pMat[0][1], pVal * pMat[1][1], pVal * pMat[2][1],
				pVal * pMat[0][2], pVal * pMat[1][2], pVal * pMat[2][2]
			);
		}

		Matrix<1, 3, T, P> operator*(const Matrix<1, 3, T, P> &pMat) const
		{
			Matrix<1, 3, T, P> lRes;
			lRes[0] = Col(
				mValues[0][0] * pMat[0][0] + mValues[1][0] * pMat[0][1] + mValues[2][0] * pMat[0][2],
				mValues[0][1] * pMat[0][0] + mValues[1][1] * pMat[0][1] + mValues[2][1] * pMat[0][2],
				mValues[0][2] * pMat[0][0] + mValues[1][2] * pMat[0][1] + mValues[2][2] * pMat[0][2]
			);
			return lRes;
		}

		Matrix<2, 3, T, P> operator*(const Matrix<2, 3, T, P> &pMat) const
		{
			Matrix<2, 3, T, P> lRes;
			lRes[0] = Col(
				mValues[0][0] * pMat[0][0] + mValues[1][0] * pMat[0][1] + mValues[2][0] * pMat[0][2],
				mValues[0][1] * pMat[0][0] + mValues[1][1] * pMat[0][1] + mValues[2][1] * pMat[0][2],
				mValues[0][2] * pMat[0][0] + mValues[1][2] * pMat[0][1] + mValues[2][2] * pMat[0][2]
			);
			lRes[1] = Col(
				mValues[0][0] * pMat[1][0] + mValues[1][0] * pMat[1][1] + mValues[2][0] * pMat[1][2],
				mValues[0][1] * pMat[1][0] + mValues[1][1] * pMat[1][1] + mValues[2][1] * pMat[1][2],
				mValues[0][2] * pMat[1][0] + mValues[1][2] * pMat[1][1] + mValues[2][2] * pMat[1][2]
			);
			return lRes;
		}

		Matrix<4, 3, T, P> operator*(const Matrix<4, 3, T, P> &pMat) const
		{
			Matrix<4, 3, T, P> lRes;
			lRes[0] = Col(
				mValues[0][0] * pMat[0][0] + mValues[1][0] * pMat[0][1] + mValues[2][0] * pMat[0][2],
				mValues[0][1] * pMat[0][0] + mValues[1][1] * pMat[0][1] + mValues[2][1] * pMat[0][2],
				mValues[0][2] * pMat[0][0] + mValues[1][2] * pMat[0][1] + mValues[2][2] * pMat[0][2]
			);
			lRes[1] = Col(
				mValues[0][0] * pMat[1][0] + mValues[1][0] * pMat[1][1] + mValues[2][0] * pMat[1][2],
				mValues[0][1] * pMat[1][0] + mValues[1][1] * pMat[1][1] + mValues[2][1] * pMat[1][2],
				mValues[0][2] * pMat[1][0] + mValues[1][2] * pMat[1][1] + mValues[2][2] * pMat[1][2]
			);
			lRes[2] = Col(
				mValues[0][0] * pMat[2][0] + mValues[1][0] * pMat[2][1] + mValues[2][0] * pMat[2][2],
				mValues[0][1] * pMat[2][0] + mValues[1][1] * pMat[2][1] + mValues[2][1] * pMat[2][2],
				mValues[0][2] * pMat[2][0] + mValues[1][2] * pMat[2][1] + mValues[2][2] * pMat[2][2]
			);
			lRes[3] = Col(
				mValues[0][0] * pMat[3][0] + mValues[1][0] * pMat[3][1] + mValues[2][0] * pMat[3][2],
				mValues[0][1] * pMat[3][0] + mValues[1][1] * pMat[3][1] + mValues[2][1] * pMat[3][2],
				mValues[0][2] * pMat[3][0] + mValues[1][2] * pMat[3][1] + mValues[2][2] * pMat[3][2]
			);
			return lRes;
		}

		template<unsigned int S>
		Matrix<S, 3, T, P> operator*(const Matrix<S, 3, T, P> &pMat) const
		{
			Matrix<S, 3, T, P> lRes;

			for (int i = 0; i < S; i++)
			{
				lRes[i] = Col(
					mValues[0][0] * pMat[i][0] + mValues[1][0] * pMat[i][1] + mValues[2][0] * pMat[i][2],
					mValues[0][1] * pMat[i][0] + mValues[1][1] * pMat[i][1] + mValues[2][1] * pMat[i][2],
					mValues[0][2] * pMat[i][0] + mValues[1][2] * pMat[i][1] + mValues[2][2] * pMat[i][2]
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
				mValues[0][0] * lInv, mValues[1][0] * lInv, mValues[2][0] * lInv,
				mValues[0][1] * lInv, mValues[1][1] * lInv, mValues[2][1] * lInv,
				mValues[0][2] * lInv, mValues[1][2] * lInv, mValues[2][2] * lInv
			);
		}

		// -- Boolean operators --

		inline bool operator==(const Mat &pMat) const
		{
			return mValues[0] == pMat[0] && mValues[1] == pMat[1] && mValues[2] == pMat[2];
		}

		inline bool operator!=(const Mat &pMat) const
		{
			return mValues[0] != pMat[0] || mValues[1] != pMat[1] || mValues[2] != pMat[2];
		}

		// -- Convertion operators --

		template<typename U, typename Q>
		inline operator Matrix<3, 3, U, Q>() const
		{
			return Matrix<3, 3, U, Q>(
				(U)mValues[0][0], (U)mValues[1][0], (U)mValues[2][0],
				(U)mValues[0][1], (U)mValues[1][1], (U)mValues[2][1],
				(U)mValues[0][2], (U)mValues[1][2], (U)mValues[2][2]
			);
		}

		// -- Stream operators --

		inline friend std::ostream &operator<<(std::ostream &pOStream, const Mat &pMat)
		{
			return pOStream	<< pMat[0][0] << ',' << pMat[1][0] << ',' << pMat[2][0] << '\n'
							<< pMat[0][1] << ',' << pMat[1][1] << ',' << pMat[2][1] << '\n'
							<< pMat[0][2] << ',' << pMat[1][2] << ',' << pMat[2][2] << '\n';
		}

		// -- Getters --

		Mat GetAdjugate() const
		{
			return Mat(
				mValues[1][1] * mValues[2][2] - mValues[2][1] * mValues[1][2],
				mValues[2][0] * mValues[1][2] - mValues[1][0] * mValues[2][2],
				mValues[1][0] * mValues[2][1] - mValues[2][0] * mValues[1][1],

				mValues[2][1] * mValues[0][2] - mValues[0][1] * mValues[2][2],
				mValues[0][0] * mValues[2][2] - mValues[2][0] * mValues[0][2],
				mValues[2][0] * mValues[0][1] - mValues[0][0] * mValues[2][1],

				mValues[0][1] * mValues[1][2] - mValues[1][1] * mValues[0][2],
				mValues[1][0] * mValues[0][2] - mValues[0][0] * mValues[1][2],
				mValues[0][0] * mValues[1][1] - mValues[1][0] * mValues[0][1]
			);
		}

		Mat GetCofactor() const
		{
			return Mat(
				mValues[1][1] * mValues[2][2] - mValues[2][1] * mValues[1][2],
				mValues[2][1] * mValues[0][2] - mValues[0][1] * mValues[2][2],
				mValues[0][1] * mValues[1][2] - mValues[1][1] * mValues[0][2],

				mValues[2][0] * mValues[1][2] - mValues[1][0] * mValues[2][2],
				mValues[0][0] * mValues[2][2] - mValues[2][0] * mValues[0][2],
				mValues[1][0] * mValues[0][2] - mValues[0][0] * mValues[1][2],

				mValues[1][0] * mValues[2][1] - mValues[2][0] * mValues[1][1],
				mValues[2][0] * mValues[0][1] - mValues[0][0] * mValues[2][1],
				mValues[0][0] * mValues[1][1] - mValues[1][0] * mValues[0][1]
			);
		}

		Mat GetDeterminant() const
		{
			return mValues[0][0] * (mValues[1][1] * mValues[2][2] - mValues[2][1] * mValues[1][2])
				-  mValues[1][0] * (mValues[0][1] * mValues[2][2] - mValues[2][1] * mValues[0][2])
				+  mValues[2][0] * (mValues[0][1] * mValues[1][2] - mValues[1][1] * mValues[0][2]);
		}

		Mat Inverted() const
		{
			return GetAdjugate() * (T(1) / GetDeterminant());
		}

		Mat Transposed() const
		{
			return Mat(
				mValues[0][0], mValues[0][1], mValues[0][2],
				mValues[1][0], mValues[1][1], mValues[1][2],
				mValues[2][0], mValues[2][1], mValues[2][2]
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

			mValues[1][2] += mValues[2][1];
			mValues[2][1] = mValues[1][2] - mValues[2][1];
			mValues[1][2] -= mValues[2][1];

			return *this;
		}

		// -- Static getters --

		inline static Mat Identity()
		{
			return Mat(T(1));
		}
	};

	typedef Matrix<3, 3, float, float> Matrix3;
	typedef Matrix<3, 3, float, float> Matrix3x3;
}

