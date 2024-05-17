#pragma once

#include "Vector.hpp"

namespace Math
{
	template<unsigned int C, unsigned int R, typename T>
	struct Matrix
	{
	private:
		typedef Matrix<C, R, T> MatCxR;
		typedef Matrix<R, C, T> MatTranspose;
		typedef Vector<R, T, float> Column;
		Column values[C];

	public:
		// -- Constructors --

		Matrix()
		{
			for (int i = 0; i < C; i++)
			{
				for (int j = 0; j < R; j++)
				{
					values[i][j] = T(0);
				}
			}
		}

		Matrix(const T value)
		{
			for (int i = 0; i < C; i++)
			{
				for (int j = 0; j < R; j++)
				{
					values[i][j] = i == j ? value : T(0);
				}
			}
		}

		Matrix(const Column values[C])
		{
			for (int i = 0; i < C; i++)
			{
				this->values[i] = values[i];
			}
		}

		Matrix(const T values[C * R])
		{
			for (int i = 0; i < C; i++)
			{
				for (int j = 0; j < R; j++)
				{
					this->values[i][j] = values[i * R + C];
				}
			}
		}

		Matrix(const MatCxR &matrix)
		{
			for (int i = 0; i < C; i++)
			{
				values[i] = matrix[i];
			}
		}

		// -- Accesses --

		const Column &operator[](const int index) const noexcept
		{
			return values[index];
		}

		Column &operator[](const int index) noexcept
		{
			return values[index];
		}

		// -- Unary arithmetic operators --

		MatCxR &operator+=(const MatCxR &matrix)
		{
			for (int i = 0; i < C; i++)
			{
				values[i] += matrix[i];
			}

			return *this;
		}

		MatCxR &operator-=(const MatCxR &matrix)
		{
			for (int i = 0; i < C; i++)
			{
				values[i] -= matrix[i];
			}

			return *this;
		}

		MatCxR &operator*=(const T scalar)
		{
			for (int i = 0; i < C; i++)
			{
				values[i] *= scalar;
			}

			return *this;
		}

		MatCxR &operator/=(const T scalar)
		{
			for (int i = 0; i < C; i++)
			{
				values[i] /= scalar;
			}

			return *this;
		}

		// -- Unary operators --

		MatCxR operator-() const
		{
			MatCxR res;

			for (int i = 0; i < C; i++)
			{
				res[i] = -values[i];
			}

			return res;
		}

		// -- Binary operators --

		MatCxR operator+(const MatCxR &matrix) const
		{
			MatCxR res;

			for (int i = 0; i < C; i++)
			{
				res[i] = values[i] + matrix[i];
			}

			return res;
		}

		MatCxR operator-(const MatCxR &matrix) const
		{
			MatCxR res;

			for (int i = 0; i < C; i++)
			{
				res[i] = values[i] - matrix[i];
			}

			return res;
		}

		template<unsigned int S>
		Matrix<S, R, T> operator*(const Matrix<S, C, T> &matrix) const
		{
			Matrix<S, R, T> res;

			for (int i = 0; i < S; i++)
			{
				for (int j = 0; j < R; j++)
				{
					for (int k = 0; k < C; k++)
					{
						res[i][j] += values[k][j] * matrix[i][k];
					}
				}
			}

			return res;
		}

		MatCxR operator*(const T scalar) const
		{
			MatCxR res;

			for (int i = 0; i < C; i++)
			{
				res[i] = values[i] * scalar;
			}

			return res;
		}

		MatCxR operator/(const T scalar) const
		{
			MatCxR res;

			for (int i = 0; i < C; i++)
			{
				res[i] = values[i] / scalar;
			}

			return res;
		}

		// -- Boolean operators --

		bool operator==(const MatCxR &matrix) const
		{
			for (int i = 0; i < C; i++)
			{
				if (values[i] != matrix[i])
				{
					return false;
				}
			}

			return true;
		}

		bool operator!=(const MatCxR &matrix) const
		{
			for (int i = 0; i < C; i++)
			{
				if (values[i] == matrix[i])
				{
					return false;
				}
			}

			return true;
		}

		// -- Stream operators --

		friend std::ostream &operator<<(std::ostream &ostream, const MatCxR &matrix)
		{
			for (int i = 0; i < C; i++)
			{
				for (int j = 0; j < R; j++)
				{
					ostream << matrix[i][j] << (j < R - 1 ? ',' : '\n');
				}
			}

			return ostream;
		}

		// -- Getters --
		MatTranspose Transposed() const
		{
			MatTranspose res;

			for (int i = 0; i < C; i++)
			{
				for (int j = 0; j < R; j++)
				{
					res[j][i] = values[i][j];
				}
			}

			return res;
		}
	};

	typedef Matrix<1, 1, float> Matrix1x1;
	typedef Matrix<2, 2, float> Matrix2;
	typedef Matrix<1, 2, float> Matrix1x2;
	typedef Matrix<2, 1, float> Matrix2x1;
	typedef Matrix<2, 2, float> Matrix2x2;
	typedef Matrix<3, 3, float> Matrix3;
	typedef Matrix<1, 3, float> Matrix1x3;
	typedef Matrix<2, 3, float> Matrix2x3;
	typedef Matrix<3, 1, float> Matrix3x1;
	typedef Matrix<3, 2, float> Matrix3x2;
	typedef Matrix<3, 3, float> Matrix3x3;
	typedef Matrix<4, 4, float> Matrix4;
	typedef Matrix<1, 4, float> Matrix1x4;
	typedef Matrix<2, 4, float> Matrix2x4;
	typedef Matrix<3, 4, float> Matrix3x4;
	typedef Matrix<4, 1, float> Matrix4x1;
	typedef Matrix<4, 2, float> Matrix4x2;
	typedef Matrix<4, 3, float> Matrix4x3;
	typedef Matrix<4, 4, float> Matrix4x4;

	template<unsigned int C, typename T>
	using MatrixNxN = Matrix<C, C, T>;

#define MATRIX_BODY(C, R, Type, Precision)\
	private:\
		typedef Matrix<C, R, Type> Mat##C##x##R;\
		typedef Matrix<R, C, Type> Mat##R##x##C;\
		typedef Matrix<R, R, Type> Mat##R##x##R;\
		typedef Vector<R, Type, Precision> Col;\
		Col values[C];\
	\
	public:\
		/* -- Constructors -- */\
		Matrix();\
		Matrix(const Type value);\
		Matrix(const Col values[C]);\
		Matrix(const Type values[C * R]);\
		Matrix(const Mat##C##x##R &matrix);\
		\
		/* -- Accesses -- */\
		const Col &operator[](const int index) const noexcept;\
		Col &operator[](const int index) noexcept;\
		\
		/* -- Unary arithmetic operators -- */\
		Mat##C##x##R &operator+=(const Mat##C##x##R &matrix);\
		Mat##C##x##R &operator-=(const Mat##C##x##R &matrix);\
		Mat##C##x##R &operator*=(const Type scalar);\
		Mat##C##x##R &operator/=(const Type scalar);\
		\
		/* -- Unary operators -- */\
		Mat##C##x##R operator-() const;\
		\
		/* -- Binary operators -- */\
		Mat##C##x##R operator+(const Mat##C##x##R &matrix) const;\
		Mat##C##x##R operator-(const Mat##C##x##R &matrix) const;\
		Mat##R##x##R operator*(const Mat##R##x##C &matrix) const;\
		template<unsigned int S>\
		Matrix<S, R, Type> operator*(const Matrix<S, C, Type> &matrix) const;\
		Mat##C##x##R operator*(const Type scalar) const;\
		Mat##C##x##R operator/(const Type scalar) const;\
		\
		/* -- Boolean operators -- */\
		bool operator==(const Mat##C##x##R &matrix) const;\
		bool operator!=(const Mat##C##x##R &matrix) const;\
		\
		/* -- Stream operators -- */\
		friend std::ostream &operator<<(std::ostream &ostream, const Mat##C##x##R &matrix);\
		\
		/* -- Getters -- */\
		Mat##R##x##C Transposed() const;

#define SQUARE_MAT_FUNCTIONS(Size, Type)\
		/* -- Unary arithmetic operators -- */\
		Mat##Size##x##Size &operator*=(const Mat##Size##x##Size &matrix);\
		Mat##Size##x##Size &operator/=(const Mat##Size##x##Size &matrix);\
		\
		/* -- Binary operators -- */\
		Mat##Size##x##Size operator/(const Mat##Size##x##Size &matrix) const;\
		\
		/* -- Getters -- */\
		Mat##Size##x##Size GetAdjugate() const;\
		Mat##Size##x##Size GetCofactor() const;\
		Type GetDeterminant() const;\
		Mat##Size##x##Size Inverted() const;\
		\
		/* -- Transformations -- */\
		Mat##Size##x##Size &Invert();\
		Mat##Size##x##Size &Transpose();
}

