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
		Col values[2];

	public:
		// -- Constructors --

		Matrix()
			: values{ Col(), Col() }
		{

		}

		Matrix(const T val)
			: values{
				Col(val, T(0)),
				Col(T(0), val)
			}
		{

		}

		Matrix(
			const T x0, const T x1,
			const T x2, const T x3
		)
			: values{
				Col(x0, x2),
				Col(x1, x3)
			}
		{

		}

		Matrix(const T vals[4])
			: values{
				Col(vals[0], vals[2]),
				Col(vals[1], vals[3])
			}
		{

		}

		Matrix(const Col &x0, const Col &x1)
			: values{
				Col(x0[0], x1[0]),
				Col(x0[1], x1[1])
			}
		{

		}

		Matrix(const Col vals[2])
			: values{
				Col(vals[0][0], vals[1][0]),
				Col(vals[0][1], vals[1][1])
			}
		{

		}

		Matrix(const Mat &mat)
			: values{ mat[0], mat[1] }
		{

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

		Mat &operator+=(const Mat &mat)
		{
			values[0] += mat[0];
			values[1] += mat[1];
			return *this;
		}

		Mat &operator-=(const Mat &mat)
		{
			values[0] -= mat[0];
			values[1] -= mat[1];
			return *this;
		}

		Mat &operator*=(const Mat &mat)
		{
			Mat temp(*this);
			values[0][0] = temp[0][0] * mat[0][0] + temp[1][0] * mat[0][1];
			values[0][1] = temp[0][1] * mat[0][0] + temp[1][1] * mat[0][1];
			values[1][0] = temp[0][0] * mat[1][0] + temp[1][0] * mat[1][1];
			values[1][1] = temp[0][1] * mat[1][0] + temp[1][1] * mat[1][1];
			return *this;
		}

		Mat &operator*=(const T val)
		{
			values[0] *= val;
			values[1] *= val;
			return *this;
		}

		Mat &operator/=(const Mat &mat)
		{
			return *this *= mat.Inverted();
		}

		Mat &operator/=(const T val)
		{
			const T inv = T(1) / val;
			values[0] *= inv;
			values[1] *= inv;
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
				-values[0][0], -values[1][0],
				-values[0][1], -values[1][1]
			);
		}

		// -- Binary operators --

		Mat operator+(const Mat &mat) const
		{
			return Mat(
				values[0][0] + mat[0][0], values[1][0] + mat[1][0],
				values[0][1] + mat[0][1], values[1][1] + mat[1][1]
			);
		}

		Mat operator-(const Mat &mat) const
		{
			return Mat(
				values[0][0] - mat[0][0], values[1][0] - mat[1][0],
				values[0][1] - mat[0][1], values[1][1] - mat[1][1]
			);
		}

		Mat operator*(const Mat &mat) const
		{
			return Mat(
				values[0][0] * mat[0][0] + values[1][0] * mat[0][1],
				values[0][0] * mat[1][0] + values[1][0] * mat[1][1],
				values[0][1] * mat[0][0] + values[1][1] * mat[0][1],
				values[0][1] * mat[1][0] + values[1][1] * mat[1][1]
			);
		}

		Col operator*(const Col &vec) const
		{
			return Col(
				values[0][0] * vec[0] + values[1][0] * vec[1],
				values[0][1] * vec[0] + values[1][1] * vec[1]
			);
		}

		Mat operator*(const T val) const
		{
			return Mat(
				values[0][0] * val, values[1][0] * val,
				values[0][1] * val, values[1][1] * val
			);
		}

		Matrix<1, 2, T, P> operator*(const Matrix<1, 2, T, P> &mat) const
		{
			Matrix<1, 2, T, P> res;
			res[0] = Col(
				values[0][0] * mat[0][0] + values[1][0] * mat[0][1],
				values[0][1] * mat[0][0] + values[1][1] * mat[0][1]
			);
			return res;
		}

		Matrix<3, 2, T, P> operator*(const Matrix<3, 2, T, P> &mat) const
		{
			Matrix<3, 2, T, P> res;
			res[0] = Col(
				values[0][0] * mat[0][0] + values[1][0] * mat[0][1],
				values[0][1] * mat[0][0] + values[1][1] * mat[0][1]
			);
			res[1] = Col(
				values[0][0] * mat[1][0] + values[1][0] * mat[1][1],
				values[0][1] * mat[1][0] + values[1][1] * mat[1][1]
			);
			res[2] = Col(
				values[0][0] * mat[2][0] + values[1][0] * mat[2][1],
				values[0][1] * mat[2][0] + values[1][1] * mat[2][1]
			);
			return res;
		}

		Matrix<4, 2, T, P> operator*(const Matrix<4, 2, T, P> &mat) const
		{
			Matrix<4, 2, T, P> res;
			res[0] = Col(
				values[0][0] * mat[0][0] + values[1][0] * mat[0][1],
				values[0][1] * mat[0][0] + values[1][1] * mat[0][1]
			);
			res[1] = Col(
				values[0][0] * mat[1][0] + values[1][0] * mat[1][1],
				values[0][1] * mat[1][0] + values[1][1] * mat[1][1]
			);
			res[2] = Col(
				values[0][0] * mat[2][0] + values[1][0] * mat[2][1],
				values[0][1] * mat[2][0] + values[1][1] * mat[2][1]
			);
			res[3] = Col(
				values[0][0] * mat[3][0] + values[1][0] * mat[3][1],
				values[0][1] * mat[3][0] + values[1][1] * mat[3][1]
			);
			return res;
		}

		template<unsigned int S>
		Matrix<S, 2, T, P> operator*(const Matrix<S, 2, T, P> &mat) const
		{
			Matrix<S, 2, T, P> res;

			for (int i = 0; i < S; i++)
			{
				res[i] = Col(
					values[0][0] * mat[i][0] + values[1][0] * mat[i][1],
					values[0][1] * mat[i][0] + values[1][1] * mat[i][1]
				);
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
			return Mat(
				values[0][0] * inv, values[1][0] * inv,
				values[0][1] * inv, values[1][1] * inv
			);
		}

		// -- Boolean operators --

		bool operator==(const Mat &mat) const
		{
			return values[0] == mat[0] && values[1] == mat[1];
		}

		bool operator!=(const Mat &mat) const
		{
			return values[0] != mat[0] || values[1] != mat[1];
		}

		// -- Convertion operators --

		template<typename U, typename Q>
		operator Matrix<2, 2, U, Q>() const
		{
			return Matrix<2, 2, U, Q>(
				(U)values[0][0], (U)values[1][0],
				(U)values[0][1], (U)values[1][1]
			);
		}

		// -- Stream operators --

		friend std::ostream &operator<<(std::ostream &ostream, const Mat &mat)
		{
			return ostream	<< mat[0][0] << ',' << mat[1][0] << '\n'
							<< mat[0][1] << ',' << mat[1][1] << '\n';
		}

		// -- Getters --

		Mat GetAdjugate() const
		{
			return Mat(
				values[1][1], -values[1][0],
				-values[0][1], values[0][0]
			);
		}

		Mat GetCofactor() const
		{
			return Mat(
				values[1][1], -values[0][1],
				-values[1][0], values[0][0]
			);
		}

		T GetDeterminant() const
		{
			return values[0][0] * values[1][1] - values[1][0] * values[0][1];
		}

		Mat Inverted() const
		{
			return GetAdjugate() * (T(1) / GetDeterminant());
		}

		Mat Transposed() const
		{
			return Mat(
				values[0][0], values[0][1],
				values[1][0], values[1][1]
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
			values[0][1] += values[1][0];
			values[1][0] = values[0][1] - values[1][0];
			values[0][1] -= values[1][0];
			return *this;
		}

		// -- Static getters --

		static inline Mat Identity()
		{
			return Mat(T(1));
		}
	};
}

