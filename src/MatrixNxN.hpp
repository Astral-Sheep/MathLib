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
		Col values[N];

	public:
		// -- Constructors --

		Matrix()
		{
			for (int i = 0; i < N; i++)
			{
				values[i] = Col();
			}
		}

		Matrix(const T val)
		{
			for (int i = 0; i < N; i++)
			{
				values[i] = Col();
				values[i][i] = val;
			}
		}

		Matrix(const T vals[N * N])
		{
			for (int i = 0; i < N; i++)
			{
				values[i] = Col();

				for (int j = 0; j < N; j++)
				{
					values[i][j] = vals[i + j * N];
				}
			}
		}

		Matrix(const Col vals[N])
		{
			for (int i = 0; i < N; i++)
			{
				values[i] = vals[i];
			}
		}

		Matrix(const Mat &mat)
		{
			for (int i = 0; i < N; i++)
			{
				values[i] = mat[i];
			}
		}

		// -- Accesses --

		inline const Col &operator[](const int index) const noexcept
		{
			return values[index];
		}

		inline Col &operator[](const int index) noexcept
		{
			return values[index];
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

		Mat &operator+=(const Mat &mat)
		{
			OPERATOR(values[i][j], +=, mat[i][j])
			return *this;
		}

		Mat &operator-=(const Mat &mat)
		{
			OPERATOR(values[i][j], -=, mat[i][j])
			return *this;
		}

		Mat &operator*=(const Mat &mat)
		{
			Mat temp(*this);
			Col column;

			for (int i = 0; i < N; i++)
			{
				for (int j = 0; j < N; j++)
				{
					for (int k = 0; k < N; k++)
					{
						column[j] += temp[k][i] * mat[j][k];
					}
				}

				values[i] = column;
			}

			return *this;
		}

		Mat &operator*=(const T val)
		{
			OPERATOR(values[i][j], *=, val)
			return *this;
		}

		Mat &operator/=(const Mat &mat)
		{
			return *this *= mat.Inverted();
		}

		Mat &operator/=(const T val)
		{
			const T inv = T(1) / val;
			OPERATOR(values[i][j], *=, inv)
			return *this;
		}

#define CONST_OPERATOR(op)\
	Mat res;\
	\
	for (int i = 0; i < N; i++)\
	{\
		for (int j = 0; j < N; j++)\
		{\
			res[i][j] = op;\
		}\
	}\
	\
	return res;

		// -- Unary operators --

		Mat operator+() const
		{
			return *this;
		}

		Mat operator-() const
		{
			CONST_OPERATOR(-values[i][j]);
		}

		// -- Binary operator --

		Mat operator+(const Mat &mat) const
		{
			CONST_OPERATOR(values[i][j] + mat[i][j])
		}

		Mat operator-(const Mat &mat) const
		{
			CONST_OPERATOR(values[i][j] - mat[i][j])
		}

		Mat operator*(const Mat &mat) const
		{
			Mat res;

			for (int i = 0; i < N; i++)
			{
				for (int j = 0; j < N; j++)
				{
					for (int k = 0; k < N; k++)
					{
						res[i][j] += values[k][j] * mat[i][k];
					}
				}
			}

			return res;
		}

		Col operator*(const Col &vec) const
		{
			Col res;

			for (int i = 0; i < N; i++)
			{
				for (int j = 0; j < N; j++)
				{
					res[i] += values[j][i] * vec[j];
				}
			}

			return res;
		}

		Mat operator*(const T val) const
		{
			CONST_OPERATOR(values[i][j] * val)
		}

		template<unsigned int S>
		Matrix<S, N, T, P> operator*(const Matrix<S, N, T, P> &mat) const
		{
			Matrix<S, N, T, P> res;

			for (int i = 0; i < S; i++)
			{
				for (int j = 0; j < N; j++)
				{
					for (int k = 0; k < N; k++)
					{
						mat[i][j] += values[k][j] * values[i][k];
					}
				}
			}

			return res;
		}

		Mat operator/(const Mat &mat) const
		{
			return *this * mat.Inverted();
		}

		Mat operator/(const T val) const
		{
			const T inv = T(1) / val;
			CONST_OPERATOR(values[i][j] * inv)
		}

		// -- Boolean operators --

#define BOOL_OPERATOR(op)\
	for (int i = 0; i < N; i++)\
	{\
		for (int j = 0; j < N; j++)\
		{\
			if (values[i][j] op mat[i][j])\
			{\
				return false;\
			}\
		}\
	}\
	\
	return true;

		bool operator==(const Mat &mat) const
		{
			BOOL_OPERATOR(!=)
		}

		bool operator!=(const Mat &mat) const
		{
			BOOL_OPERATOR(==)
		}

		// -- Convertion operators --

		template<typename U, typename Q>
		operator Matrix<N, N, U, Q>() const
		{
			CONST_OPERATOR((U)values[i][j])
		}

		// -- Stream operators --

		friend std::ostream &operator<<(std::ostream &ostream, const Mat &mat)
		{
			for (int i = 0; i < N; i++)
			{
				for (int j = 0; j < N; j++)
				{
					ostream << mat[j][i] << (j < N - 1 ? ',' : '\n');
				}
			}

			return ostream;
		}

		// -- Getters --

		Mat GetAdjugate() const
		{
			Mat res;

			for (int i = 0; i < N; i++)
			{
				for (int j = 0; j < N; j++)
				{
					Matrix<N - 1, N - 1, T, P> submat;
					int x = 0;
					int y = 0;

					for (int k = 0; k < N; k++)
					{
						if (k == i)
						{
							continue;
						}

						for (int l = 0; k < N; k++)
						{
							if (l == j)
							{
								continue;
							}

							submat[x][y] = values[i][j];
							y++;
						}

						x++;
					}

					res[j][i] = submat.GetDeterminant() * (((i + j) % 2) ? -1 : 1);
				}
			}

			return res;
		}

		Mat GetCofactor() const
		{
			Mat res;

			for (int i = 0; i < N; i++)
			{
				for (int j = 0; j < N; j++)
				{
					Matrix<N - 1, N - 1, T, P> submat;
					int x = 0;
					int y = 0;

					for (int k = 0; k < N; k++)
					{
						if (k == i)
						{
							continue;
						}

						for (int l = 0; k < N; k++)
						{
							if (l == j)
							{
								continue;
							}

							submat[x][y] = values[i][j];
							y++;
						}

						x++;
					}

					res[i][j] = submat.GetDeterminant() * (((i + j) % 2) ? -1 : 1);
				}
			}

			return res;
		}

		T GetDeterminant() const
		{
			if (N == 0)
			{
				return 0;
			}

			if (N == 1)
			{
				return values[0][0];
			}

			if (N == 2)
			{
				return values[0][0] * values[1][1] - values[1][0] * values[0][1];
			}

			// Matrix3x3 and above

			T determinant = T(0);
			int x;
			int y;

			for (int i = 0; i < N; i++)
			{
				Matrix<N - 1, N - 1, T, P> subMat;
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

						subMat[x][y] = values[j][k];
						y++;
					}

					x++;
				}

				determinant += values[0][i] * subMat.GetDeterminant() * ((i % 2) ? -1 : 1);
			}

			return determinant;
		}

		Mat Inverted() const
		{
			return GetAdjugate() * (T(1) / GetDeterminant());
		}

		Mat Transposed() const
		{
			Mat res;

			for (int i = 0; i < N; i++)
			{
				for (int j = 0; j < N; j++)
				{
					res[i][j] = values[j][i];
				}
			}

			return res;
		}

		// -- Transformations --

		Mat &Invert()
		{
			*this = GetAdjugate() * (T(1) / GetDeterminant());
			return *this;
		}

		Mat &Transpose()
		{
			const Mat temp(*this);

			for (int i = 0; i < N; i++)
			{
				for (int j = 0; j < N; j++)
				{
					values[i][j] = temp[j][i];
				}
			}

			return *this;
		}

		// -- Static getters --

		static inline Mat Identity()
		{
			return Mat(T(1));
		}
	};

#undef OPERATOR
#undef CONST_OPERATOR
#undef BOOL_OPERATOR
}

