#pragma once

#include "Matrix.hpp"
#include "Vector.hpp"

namespace Math
{
	template<unsigned int N, typename T, typename P>
	struct Matrix<N, N, T, P>
	{
	private:
		typedef Matrix<N, N, T, P> Mat;
		typedef Vector<N, T, P> Col;
		Col mValues[N];

	public:
		// -- Constructors --

		Matrix()
		{
			for (int i = 0; i < N; i++)
			{
				mValues[i] = Col();
			}
		}

		explicit Matrix(const T pVal)
		{
			for (int i = 0; i < N; i++)
			{
				mValues[i] = Col();
				mValues[i][i] = pVal;
			}
		}

		Matrix(const T pVals[N * N])
		{
			for (int i = 0; i < N; i++)
			{
				mValues[i] = Col();

				for (int j = 0; j < N; j++)
				{
					mValues[i][j] = pVals[i + j * N];
				}
			}
		}

		Matrix(const Col pVals[N])
		{
			for (int i = 0; i < N; i++)
			{
				mValues[i] = pVals[i];
			}
		}

		Matrix(const Mat &pMat)
		{
			for (int i = 0; i < N; i++)
			{
				mValues[i] = pMat[i];
			}
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

#define OPERATOR(lhs, op, rhs)\
	for (int i = 0; i < N; i++)\
	{\
		for (int j = 0; j < N; j++)\
		{\
			lhs op rhs;\
		}\
	}

		Mat &operator+=(const Mat &pMat)
		{
			OPERATOR(mValues[i][j], +=, pMat[i][j])
			return *this;
		}

		Mat &operator-=(const Mat &pMat)
		{
			OPERATOR(mValues[i][j], -=, pMat[i][j])
			return *this;
		}

		Mat &operator*=(const Mat &pMat)
		{
			const Mat lTemp(*this);
			Col lColumn;

			for (int i = 0; i < N; i++)
			{
				for (int j = 0; j < N; j++)
				{
					for (int k = 0; k < N; k++)
					{
						lColumn[j] += lTemp[k][i] * pMat[j][k];
					}
				}

				mValues[i] = lColumn;
			}

			return *this;
		}

		Mat &operator*=(const T pVal)
		{
			OPERATOR(mValues[i][j], *=, pVal)
			return *this;
		}

		Mat &operator/=(const Mat &pMat)
		{
			return *this *= pMat.Inverted();
		}

		Mat &operator/=(const T pVal)
		{
			const T lInv = T(1) / pVal;
			OPERATOR(mValues[i][j], *=, lInv)
			return *this;
		}

#define CONST_OPERATOR(op)\
	Mat lRes;\
	\
	for (int i = 0; i < N; i++)\
	{\
		for (int j = 0; j < N; j++)\
		{\
			lRes[i][j] = op;\
		}\
	}\
	\
	return lRes;

		// -- Unary operators --

		Mat operator+() const
		{
			return *this;
		}

		Mat operator-() const
		{
			CONST_OPERATOR(-mValues[i][j]);
		}

		// -- Binary operator --

		Mat operator+(const Mat &pMat) const
		{
			CONST_OPERATOR(mValues[i][j] + pMat[i][j])
		}

		Mat operator-(const Mat &pMat) const
		{
			CONST_OPERATOR(mValues[i][j] - pMat[i][j])
		}

		Mat operator*(const Mat &pMat) const
		{
			Mat lRes;

			for (int i = 0; i < N; i++)
			{
				for (int j = 0; j < N; j++)
				{
					for (int k = 0; k < N; k++)
					{
						lRes[i][j] += mValues[k][j] * pMat[i][k];
					}
				}
			}

			return lRes;
		}

		Col operator*(const Col &pVec) const
		{
			Col lRes;

			for (int i = 0; i < N; i++)
			{
				for (int j = 0; j < N; j++)
				{
					lRes[i] += mValues[j][i] * pVec[j];
				}
			}

			return lRes;
		}

		Mat operator*(const T pVal) const
		{
			CONST_OPERATOR(mValues[i][j] * pVal)
		}

		template<unsigned int S>
		Matrix<S, N, T, P> operator*(const Matrix<S, N, T, P> &pMat) const
		{
			Matrix<S, N, T, P> lRes;

			for (int i = 0; i < S; i++)
			{
				for (int j = 0; j < N; j++)
				{
					for (int k = 0; k < N; k++)
					{
						pMat[i][j] += mValues[k][j] * mValues[i][k];
					}
				}
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
			CONST_OPERATOR(mValues[i][j] * lInv)
		}

		// -- Boolean operators --

#define BOOL_OPERATOR(op)\
	for (int i = 0; i < N; i++)\
	{\
		for (int j = 0; j < N; j++)\
		{\
			if (mValues[i][j] op pMat[i][j])\
			{\
				return false;\
			}\
		}\
	}\
	\
	return true;

		bool operator==(const Mat &pMat) const
		{
			BOOL_OPERATOR(!=)
		}

		bool operator!=(const Mat &pMat) const
		{
			BOOL_OPERATOR(==)
		}

		// -- Convertion operators --

		template<typename U, typename Q>
		operator Matrix<N, N, U, Q>() const
		{
			Matrix<N, N, U, Q> lRes;

			for (int i = 0; i < N; i++)
			{
				for (int j = 0; j < N; j++)
				{
					lRes[i][j] = (U)mValues[i][j];
				}
			}

			return lRes;
		}

		// -- Stream operators --

		friend std::ostream &operator<<(std::ostream &pOStream, const Mat &pMat)
		{
			for (int i = 0; i < N; i++)
			{
				for (int j = 0; j < N; j++)
				{
					pOStream << pMat[j][i] << (j < N - 1 ? ',' : '\n');
				}
			}

			return pOStream;
		}

		// -- Getters --

		Mat GetAdjugate() const
		{
			Mat lRes;

			for (int i = 0; i < N; i++)
			{
				for (int j = 0; j < N; j++)
				{
					Matrix<N - 1, N - 1, T, P> lSubmat;
					int x = 0;
					int y = 0;

					for (int k = 0; k < N; k++)
					{
						if (k == i)
						{
							continue;
						}

						for (int l = 0; l < N; l++)
						{
							if (l == j)
							{
								continue;
							}

							lSubmat[x][y] = mValues[i][j];
							y++;
						}

						x++;
					}

					lRes[j][i] = lSubmat.GetDeterminant() * (((i + j + 1) % 2) * 2 - 1);
				}
			}

			return lRes;
		}

		Mat GetCofactor() const
		{
			Mat lRes;

			for (int i = 0; i < N; i++)
			{
				for (int j = 0; j < N; j++)
				{
					Matrix<N - 1, N - 1, T, P> lSubmat;
					int x = 0;
					int y = 0;

					for (int k = 0; k < N; k++)
					{
						if (k == i)
						{
							continue;
						}

						for (int l = 0; l < N; l++)
						{
							if (l == j)
							{
								continue;
							}

							lSubmat[x][y] = mValues[i][j];
							y++;
						}

						x++;
					}

					lRes[i][j] = lSubmat.GetDeterminant() * (((i + j + 1) % 2) * 2 - 1);
				}
			}

			return lRes;
		}

		T GetDeterminant() const
		{
			if (N == 0)
			{
				return 0;
			}

			if (N == 1)
			{
				return mValues[0][0];
			}

			if (N == 2)
			{
				return mValues[0][0] * mValues[1][1] - mValues[1][0] * mValues[0][1];
			}

			// Matrix3x3 and above

			T lDet = T(0);
			int x;
			int y;

			for (int i = 0; i < N; i++)
			{
				Matrix<N - 1, N - 1, T, P> lSubmat;
				x = 0;

				for (int j = 0; j < N; j++)
				{
					if (j == i)
					{
						continue;
					}

					y = 0;

					for (int k = 0; k < N; k++)
					{
						if (k == 0)
						{
							continue;
						}

						lSubmat[x][y] = mValues[j][k];
						y++;
					}

					x++;
				}

				lDet += mValues[0][i] * lSubmat.GetDeterminant() * ((i % 2) ? -1 : 1);
			}

			return lDet;
		}

		Mat Inverted() const
		{
			return GetAdjugate() * (T(1) / GetDeterminant());
		}

		Mat Transposed() const
		{
			CONST_OPERATOR(mValues[j][i]);
		}

		// -- Transformations --

		Mat &Invert()
		{
			*this = GetAdjugate() * (T(1) / GetDeterminant());
			return *this;
		}

		Mat &Transpose()
		{
			const Mat lTemp(*this);

			for (int i = 0; i < N; i++)
			{
				for (int j = 0; j < N; j++)
				{
					mValues[i][j] = lTemp[j][i];
				}
			}

			return *this;
		}

		// -- Static getters --

		static inline Mat Identity()
		{
			return Mat(T(1));
		}
#undef OPERATOR
#undef CONST_OPERATOR
#undef BOOL_OPERATOR
	};

	template<unsigned int N, typename T, typename P>
	using MatrixNxN = Matrix<N, N, T, P>;

	template<unsigned int N, typename T, typename P>
	using MatrixN = Matrix<N, N, T, P>;
}

