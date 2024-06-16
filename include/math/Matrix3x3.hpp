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
		Col values[3];

	public:
		// -- Constructors --

		Matrix()
			: values{ Col(), Col(), Col() }
		{

		}

		Matrix(const T val)
			: values{
				Col(val, T(0), T(0)),
				Col(T(0), val, T(0)),
				Col(T(0), T(0), val)
			}
		{

		}

		Matrix(
			const T x0, const T x1, const T x2,
			const T x3, const T x4, const T x5,
			const T x6, const T x7, const T x8
		)
			: values{
				Col(x0, x3, x6),
				Col(x1, x4, x7),
				Col(x2, x5, x8)
			}
		{

		}

		Matrix(const T vals[9])
			: values{
				Col(vals[0], vals[3], vals[6]),
				Col(vals[1], vals[4], vals[7]),
				Col(vals[2], vals[5], vals[8])
			}
		{

		}

		Matrix(const Col &x0, const Col &x1, const Col &x2)
			: values{
				Col(x0[0], x1[0], x2[0]),
				Col(x0[1], x1[1], x2[1]),
				Col(x0[2], x1[2], x2[2])
			}
		{

		}

		Matrix(const Col vals[3])
			: values{
				Col(vals[0][0], vals[1][0], vals[2][0]),
				Col(vals[0][1], vals[1][1], vals[2][1]),
				Col(vals[0][2], vals[1][2], vals[2][2])
			}
		{

		}

		Matrix(const Mat &mat)
			: values{ mat[0], mat[1], mat[2] }
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
			values[2] += mat[2];
			return *this;
		}

		Mat &operator-=(const Mat &mat)
		{
			values[0] -= mat[0];
			values[1] -= mat[1];
			values[2] -= mat[2];
			return *this;
		}

		Mat &operator*=(const Mat &mat)
		{
			Mat temp(*this);
			values[0][0] = temp[0][0] * mat[0][0] + temp[1][0] * mat[0][1] + temp[2][0] * mat[0][2];
			values[0][1] = temp[0][1] * mat[0][0] + temp[1][1] * mat[0][1] + temp[2][1] * mat[0][2];
			values[0][2] = temp[0][2] * mat[0][0] + temp[1][2] * mat[0][1] + temp[2][2] * mat[0][2];

			values[1][0] = temp[0][0] * mat[1][0] + temp[1][0] * mat[1][1] + temp[2][0] * mat[1][2];
			values[1][1] = temp[0][1] * mat[1][0] + temp[1][1] * mat[1][1] + temp[2][1] * mat[1][2];
			values[1][2] = temp[0][2] * mat[1][0] + temp[1][2] * mat[1][1] + temp[2][2] * mat[1][2];

			values[2][0] = temp[0][0] * mat[2][0] + temp[1][0] * mat[2][1] + temp[2][0] * mat[2][2];
			values[2][1] = temp[0][1] * mat[2][0] + temp[1][1] * mat[2][1] + temp[2][1] * mat[2][2];
			values[2][2] = temp[0][2] * mat[2][0] + temp[1][2] * mat[2][1] + temp[2][2] * mat[2][2];
			return *this;
		}

		Mat &operator*=(const T val)
		{
			values[0] *= val;
			values[1] *= val;
			values[2] *= val;
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
			values[2] *= inv;
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
				-values[0][0], -values[1][0], -values[2][0],
				-values[0][1], -values[1][1], -values[2][1],
				-values[0][2], -values[1][2], -values[2][2]
			);
		}

		// -- Binary operators --

		Mat operator+(const Mat &mat) const
		{
			return Mat(
				values[0][0] + mat[0][0], values[1][0] + mat[1][0], values[2][0] + mat[2][0],
				values[0][1] + mat[0][1], values[1][1] + mat[1][1], values[2][1] + mat[2][1],
				values[0][2] + mat[0][2], values[1][2] + mat[1][2], values[2][2] + mat[2][2]
			);
		}

		Mat operator-(const Mat &mat) const
		{
			return Mat(
				values[0][0] - mat[0][0], values[1][0] - mat[1][0], values[2][0] - mat[2][0],
				values[0][1] - mat[0][1], values[1][1] - mat[1][1], values[2][1] - mat[2][1],
				values[0][2] - mat[0][2], values[1][2] - mat[1][2], values[2][2] - mat[2][2]
			);
		}

		Mat operator*(const Mat &mat) const
		{
			return Mat(
				values[0][0] * mat[0][0] + values[1][0] * mat[0][1] + values[2][0] * mat[0][2],
				values[0][0] * mat[1][0] + values[1][0] * mat[1][1] + values[2][0] * mat[1][2],
				values[0][0] * mat[2][0] + values[1][0] * mat[2][1] + values[2][0] * mat[2][2],

				values[0][1] * mat[0][0] + values[1][1] * mat[0][1] + values[2][1] * mat[0][2],
				values[0][1] * mat[1][0] + values[1][1] * mat[1][1] + values[2][1] * mat[1][2],
				values[0][1] * mat[2][0] + values[1][1] * mat[2][1] + values[2][1] * mat[2][2],

				values[0][2] * mat[0][0] + values[1][2] * mat[0][1] + values[2][2] * mat[0][2],
				values[0][2] * mat[1][0] + values[1][2] * mat[1][1] + values[2][2] * mat[1][2],
				values[0][2] * mat[2][0] + values[1][2] * mat[2][1] + values[2][2] * mat[2][2]
			);
		}

		Col operator*(const Col &vec) const
		{
			return Col(
				values[0][0] * vec.x + values[1][0] * vec.y + values[2][0] * vec.z,
				values[0][1] * vec.x + values[1][1] * vec.y + values[2][1] * vec.z,
				values[0][2] * vec.x + values[1][2] * vec.y + values[2][2] * vec.z
			);
		}

		Mat operator*(const T val) const
		{
			return Mat(
				values[0][0] * val, values[1][0] * val, values[2][0] * val,
				values[0][1] * val, values[1][1] * val, values[2][1] * val,
				values[0][2] * val, values[1][2] * val, values[2][2] * val
			);
		}

		Matrix<1, 3, T, P> operator*(const Matrix<1, 3, T, P> &mat) const
		{
			Matrix<1, 3, T, P> res;
			res[0] = Col(
				values[0][0] * mat[0][0] + values[1][0] * mat[0][1] + values[2][0] * mat[0][2],
				values[0][1] * mat[0][0] + values[1][1] * mat[0][1] + values[2][1] * mat[0][2],
				values[0][2] * mat[0][0] + values[1][2] * mat[0][1] + values[2][2] * mat[0][2]
			);
			return res;
		}

		Matrix<2, 3, T, P> operator*(const Matrix<2, 3, T, P> &mat) const
		{
			Matrix<2, 3, T, P> res;
			res[0] = Col(
				values[0][0] * mat[0][0] + values[1][0] * mat[0][1] + values[2][0] * mat[0][2],
				values[0][1] * mat[0][0] + values[1][1] * mat[0][1] + values[2][1] * mat[0][2],
				values[0][2] * mat[0][0] + values[1][2] * mat[0][1] + values[2][2] * mat[0][2]
			);
			res[1] = Col(
				values[0][0] * mat[1][0] + values[1][0] * mat[1][1] + values[2][0] * mat[1][2],
				values[0][1] * mat[1][0] + values[1][1] * mat[1][1] + values[2][1] * mat[1][2],
				values[0][2] * mat[1][0] + values[1][2] * mat[1][1] + values[2][2] * mat[1][2]
			);
			return res;
		}

		Matrix<4, 3, T, P> operator*(const Matrix<4, 3, T, P> &mat) const
		{
			Matrix<4, 3, T, P> res;
			res[0] = Col(
				values[0][0] * mat[0][0] + values[1][0] * mat[0][1] + values[2][0] * mat[0][2],
				values[0][1] * mat[0][0] + values[1][1] * mat[0][1] + values[2][1] * mat[0][2],
				values[0][2] * mat[0][0] + values[1][2] * mat[0][1] + values[2][2] * mat[0][2]
			);
			res[1] = Col(
				values[0][0] * mat[1][0] + values[1][0] * mat[1][1] + values[2][0] * mat[1][2],
				values[0][1] * mat[1][0] + values[1][1] * mat[1][1] + values[2][1] * mat[1][2],
				values[0][2] * mat[1][0] + values[1][2] * mat[1][1] + values[2][2] * mat[1][2]
			);
			res[2] = Col(
				values[0][0] * mat[2][0] + values[1][0] * mat[2][1] + values[2][0] * mat[2][2],
				values[0][1] * mat[2][0] + values[1][1] * mat[2][1] + values[2][1] * mat[2][2],
				values[0][2] * mat[2][0] + values[1][2] * mat[2][1] + values[2][2] * mat[2][2]
			);
			res[3] = Col(
				values[0][0] * mat[3][0] + values[1][0] * mat[3][1] + values[2][0] * mat[3][2],
				values[0][1] * mat[3][0] + values[1][1] * mat[3][1] + values[2][1] * mat[3][2],
				values[0][2] * mat[3][0] + values[1][2] * mat[3][1] + values[2][2] * mat[3][2]
			);
			return res;
		}

		template<unsigned int S>
		Matrix<S, 3, T, P> operator*(const Matrix<S, 3, T, P> &mat) const
		{
			Matrix<S, 3, T, P> res;

			for (int i = 0; i < S; i++)
			{
				res[i] = Col(
					values[0][0] * mat[i][0] + values[1][0] * mat[i][1] + values[2][0] * mat[i][2],
					values[0][1] * mat[i][0] + values[1][1] * mat[i][1] + values[2][1] * mat[i][2],
					values[0][2] * mat[i][0] + values[1][2] * mat[i][1] + values[2][2] * mat[i][2]
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
				values[0][0] * inv, values[1][0] * inv, values[2][0] * inv,
				values[0][1] * inv, values[1][1] * inv, values[2][1] * inv,
				values[0][2] * inv, values[1][2] * inv, values[2][2] * inv
			);
		}

		// -- Boolean operators --

		bool operator==(const Mat &mat) const
		{
			return values[0] == mat[0] && values[1] == mat[1] && values[2] == mat[2];
		}

		bool operator!=(const Mat &mat) const
		{
			return values[0] != mat[0] || values[1] != mat[1] || values[2] != mat[2];
		}

		// -- Convertion operators --

		template<typename U, typename Q>
		operator Matrix<3, 3, U, Q>() const
		{
			return Matrix<3, 3, U, Q>(
				(U)values[0][0], (U)values[1][0], (U)values[2][0],
				(U)values[0][1], (U)values[1][1], (U)values[2][1],
				(U)values[0][2], (U)values[1][2], (U)values[2][2]
			);
		}

		// -- Stream operators --

		friend std::ostream &operator<<(std::ostream &ostream, const Mat &mat)
		{
			return ostream	<< mat[0][0] << ',' << mat[1][0] << ',' << mat[2][0] << '\n'
							<< mat[0][1] << ',' << mat[1][1] << ',' << mat[2][1] << '\n'
							<< mat[0][2] << ',' << mat[1][2] << ',' << mat[2][2] << '\n';
		}

		// -- Getters --

		Mat GetAdjugate() const
		{
			return Mat(
				values[1][1] * values[2][2] - values[2][1] * values[1][2],
				values[2][0] * values[1][2] - values[1][0] * values[2][2],
				values[1][0] * values[2][1] - values[2][0] * values[1][1],

				values[2][1] * values[0][2] - values[0][1] * values[2][2],
				values[0][0] * values[2][2] - values[2][0] * values[0][2],
				values[2][0] * values[0][1] - values[0][0] * values[2][1],

				values[0][1] * values[1][2] - values[1][1] * values[0][2],
				values[1][0] * values[0][2] - values[0][0] * values[1][2],
				values[0][0] * values[1][1] - values[1][0] * values[0][1]
			);
		}

		Mat GetCofactor() const
		{
			return Mat(
				values[1][1] * values[2][2] - values[2][1] * values[1][2],
				values[2][1] * values[0][2] - values[0][1] * values[2][2],
				values[0][1] * values[1][2] - values[1][1] * values[0][2],

				values[2][0] * values[1][2] - values[1][0] * values[2][2],
				values[0][0] * values[2][2] - values[2][0] * values[0][2],
				values[1][0] * values[0][2] - values[0][0] * values[1][2],

				values[1][0] * values[2][1] - values[2][0] * values[1][1],
				values[2][0] * values[0][1] - values[0][0] * values[2][1],
				values[0][0] * values[1][1] - values[1][0] * values[0][1]
			);
		}

		Mat GetDeterminant() const
		{
			return values[0][0] * (values[1][1] * values[2][2] - values[2][1] * values[1][2])
				-  values[1][0] * (values[0][1] * values[2][2] - values[2][1] * values[0][2])
				+  values[2][0] * (values[0][1] * values[1][2] - values[1][1] * values[0][2]);
		}

		Mat Inverted() const
		{
			return GetAdjugate() * (T(1) / GetDeterminant());
		}

		Mat Transposed() const
		{
			return Mat(
				values[0][0], values[0][1], values[0][2],
				values[1][0], values[1][1], values[1][2],
				values[2][0], values[2][1], values[2][2]
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

			values[0][2] += values[2][0];
			values[2][0] = values[0][2] - values[2][0];
			values[0][2] -= values[2][0];

			values[1][2] += values[2][1];
			values[2][1] = values[1][2] - values[2][1];
			values[1][2] -= values[2][1];

			return *this;
		}

		// -- Static getters --

		static inline Mat Identity()
		{
			return Mat(T(1));
		}
	};
}

