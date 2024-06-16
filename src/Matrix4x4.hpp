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
		Col values[4];

	public:
		// -- Constructors --

		Matrix()
			: values{ Col(), Col(), Col(), Col() }
		{

		}

		Matrix(const T val)
			: values{
				Col(val, T(0), T(0), T(0)),
				Col(T(0), val, T(0), T(0)),
				Col(T(0), T(0), val, T(0)),
				Col(T(0), T(0), T(0), val)
			}
		{

		}

		Matrix(
			const T x0, const T x1, const T x2, const T x3,
			const T x4, const T x5, const T x6, const T x7,
			const T x8, const T x9, const T x10, const T x11,
			const T x12, const T x13, const T x14, const T x15
		)
			: values{
				Col(x0, x4, x8, x12),
				Col(x1, x5, x9, x13),
				Col(x2, x6, x10, x14),
				Col(x3, x7, x11, x15)
			}
		{

		}

		Matrix(const T vals[16])
			: values{
				Col(vals[0], vals[4], vals[8], vals[12]),
				Col(vals[1], vals[5], vals[9], vals[13]),
				Col(vals[2], vals[6], vals[10], vals[14]),
				Col(vals[3], vals[7], vals[11], vals[15])
			}
		{

		}

		Matrix(const Col &x0, const Col &x1, const Col &x2, const Col &x3)
			: values{
				Col(x0[0], x1[0], x2[0], x3[0]),
				Col(x0[1], x1[1], x2[1], x3[1]),
				Col(x0[2], x1[2], x2[2], x3[2]),
				Col(x0[3], x1[3], x2[3], x3[3]),
			}
		{

		}

		Matrix(const Col vals[4])
			: values{
				Col(vals[0][0], vals[1][0], vals[2][0], vals[3][0]),
				Col(vals[0][1], vals[1][1], vals[2][1], vals[3][1]),
				Col(vals[0][2], vals[1][2], vals[2][2], vals[3][2]),
				Col(vals[0][3], vals[1][3], vals[2][3], vals[3][3])
			}
		{

		}

		Matrix (const Mat &mat)
			: values{ mat[0], mat[1], mat[2], mat[3] }
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
			values[3] += mat[3];
			return *this;
		}

		Mat &operator-=(const Mat &mat)
		{
			values[0] -= mat[0];
			values[1] -= mat[1];
			values[2] -= mat[2];
			values[3] -= mat[3];
			return *this;
		}

		Mat &operator*=(const Mat &mat)
		{
			Mat temp(*this);

			values[0][0] = temp[0][0] * mat[0][0] + temp[1][0] * mat[0][1] + temp[2][0] * mat[0][2] + temp[3][0] * mat[0][3];
			values[0][1] = temp[0][1] * mat[0][0] + temp[1][1] * mat[0][1] + temp[2][1] * mat[0][2] + temp[3][1] * mat[0][3];
			values[0][2] = temp[0][2] * mat[0][0] + temp[1][2] * mat[0][1] + temp[2][2] * mat[0][2] + temp[3][2] * mat[0][3];
			values[0][3] = temp[0][3] * mat[0][0] + temp[1][3] * mat[0][1] + temp[2][3] * mat[0][2] + temp[3][3] * mat[0][3];

			values[1][0] = temp[0][0] * mat[1][0] + temp[1][0] * mat[1][1] + temp[2][0] * mat[1][2] + temp[3][0] * mat[1][3];
			values[1][1] = temp[0][1] * mat[1][0] + temp[1][1] * mat[1][1] + temp[2][1] * mat[1][2] + temp[3][1] * mat[1][3];
			values[1][2] = temp[0][2] * mat[1][0] + temp[1][2] * mat[1][1] + temp[2][2] * mat[1][2] + temp[3][2] * mat[1][3];
			values[1][3] = temp[0][3] * mat[1][0] + temp[1][3] * mat[1][1] + temp[2][3] * mat[1][2] + temp[3][3] * mat[1][3];

			values[2][0] = temp[0][0] * mat[2][0] + temp[1][0] * mat[2][1] + temp[2][0] * mat[2][2] + temp[3][0] * mat[2][3];
			values[2][1] = temp[0][1] * mat[2][0] + temp[1][1] * mat[2][1] + temp[2][1] * mat[2][2] + temp[3][1] * mat[2][3];
			values[2][2] = temp[0][2] * mat[2][0] + temp[1][2] * mat[2][1] + temp[2][2] * mat[2][2] + temp[3][2] * mat[2][3];
			values[2][3] = temp[0][3] * mat[2][0] + temp[1][3] * mat[2][1] + temp[2][3] * mat[2][2] + temp[3][3] * mat[2][3];

			values[3][0] = temp[0][0] * mat[3][0] + temp[1][0] * mat[3][1] + temp[2][0] * mat[3][2] + temp[3][0] * mat[3][3];
			values[3][1] = temp[0][1] * mat[3][0] + temp[1][1] * mat[3][1] + temp[2][1] * mat[3][2] + temp[3][1] * mat[3][3];
			values[3][2] = temp[0][2] * mat[3][0] + temp[1][2] * mat[3][1] + temp[2][2] * mat[3][2] + temp[3][2] * mat[3][3];
			values[3][3] = temp[0][3] * mat[3][0] + temp[1][3] * mat[3][1] + temp[2][3] * mat[3][2] + temp[3][3] * mat[3][3];

			return *this;
		}

		Mat &operator*=(const T val)
		{
			values[0] *= val;
			values[1] *= val;
			values[2] *= val;
			values[3] *= val;
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
			values[3] *= inv;
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
				-values[0][0], -values[1][0], -values[2][0], -values[3][0],
				-values[0][1], -values[1][1], -values[2][1], -values[3][1],
				-values[0][2], -values[1][2], -values[2][2], -values[3][2],
				-values[0][3], -values[1][3], -values[2][3], -values[3][3]
			);
		}

		// -- Binary operators --

		Mat operator+(const Mat &mat) const
		{
			return Mat(
				values[0][0] + mat[0][0], values[1][0] + mat[1][0], values[2][0] + mat[2][0], values[3][0] + mat[3][0],
				values[0][1] + mat[0][1], values[1][1] + mat[1][1], values[2][1] + mat[2][1], values[3][1] + mat[3][1],
				values[0][2] + mat[0][2], values[1][2] + mat[1][2], values[2][2] + mat[2][2], values[3][2] + mat[3][2],
				values[0][3] + mat[0][3], values[1][3] + mat[1][3], values[2][3] + mat[2][3], values[3][3] + mat[3][3]
			);
		}

		Mat operator-(const Mat &mat) const
		{
			return Mat(
				values[0][0] - mat[0][0], values[1][0] - mat[1][0], values[2][0] - mat[2][0], values[3][0] - mat[3][0],
				values[0][1] - mat[0][1], values[1][1] - mat[1][1], values[2][1] - mat[2][1], values[3][1] - mat[3][1],
				values[0][2] - mat[0][2], values[1][2] - mat[1][2], values[2][2] - mat[2][2], values[3][2] - mat[3][2],
				values[0][3] - mat[0][3], values[1][3] - mat[1][3], values[2][3] - mat[2][3], values[3][3] - mat[3][3]
			);
		}

		Mat operator*(const Mat &mat) const
		{
			return Mat(
				// Row 0
				values[0][0] * mat[0][0] + values[1][0] * mat[0][1] + values[2][0] * mat[0][2] + values[3][0] * mat[0][3],
				values[0][0] * mat[1][0] + values[1][0] * mat[1][1] + values[2][0] * mat[1][2] + values[3][0] * mat[1][3],
				values[0][0] * mat[2][0] + values[1][0] * mat[2][1] + values[2][0] * mat[2][2] + values[3][0] * mat[2][3],
				values[0][0] * mat[3][0] + values[1][0] * mat[3][1] + values[2][0] * mat[3][2] + values[3][0] * mat[3][3],

				// Row 1
				values[0][1] * mat[0][0] + values[1][1] * mat[0][1] + values[2][1] * mat[0][2] + values[3][1] * mat[0][3],
				values[0][1] * mat[1][0] + values[1][1] * mat[1][1] + values[2][1] * mat[1][2] + values[3][1] * mat[1][3],
				values[0][1] * mat[2][0] + values[1][1] * mat[2][1] + values[2][1] * mat[2][2] + values[3][1] * mat[2][3],
				values[0][1] * mat[3][0] + values[1][1] * mat[3][1] + values[2][1] * mat[3][2] + values[3][1] * mat[3][3],

				// Row 2
				values[0][2] * mat[0][0] + values[1][2] * mat[0][1] + values[2][2] * mat[0][2] + values[3][2] * mat[0][3],
				values[0][2] * mat[1][0] + values[1][2] * mat[1][1] + values[2][2] * mat[1][2] + values[3][2] * mat[1][3],
				values[0][2] * mat[2][0] + values[1][2] * mat[2][1] + values[2][2] * mat[2][2] + values[3][2] * mat[2][3],
				values[0][2] * mat[3][0] + values[1][2] * mat[3][1] + values[2][2] * mat[3][2] + values[3][2] * mat[3][3],

				// Row 3
				values[0][3] * mat[0][0] + values[1][3] * mat[0][1] + values[2][3] * mat[0][2] + values[3][3] * mat[0][3],
				values[0][3] * mat[1][0] + values[1][3] * mat[1][1] + values[2][3] * mat[1][2] + values[3][3] * mat[1][3],
				values[0][3] * mat[2][0] + values[1][3] * mat[2][1] + values[2][3] * mat[2][2] + values[3][3] * mat[2][3],
				values[0][3] * mat[3][0] + values[1][3] * mat[3][1] + values[2][3] * mat[3][2] + values[3][3] * mat[3][3]
			);
		}

		Col operator*(const Col &vec) const
		{
			return Col(
				values[0][0] * vec.x + values[1][0] * vec.y + values[2][0] * vec.z + values[3][0] * vec.w,
				values[0][1] * vec.x + values[1][1] * vec.y + values[2][1] * vec.z + values[3][1] * vec.w,
				values[0][2] * vec.x + values[1][2] * vec.y + values[2][2] * vec.z + values[3][2] * vec.w,
				values[0][3] * vec.x + values[1][3] * vec.y + values[2][3] * vec.z + values[3][3] * vec.w
			);
		}

		Mat operator*(const T val) const
		{
			return Mat(
				values[0][0] * val, values[1][0] * val, values[2][0] * val, values[3][0] * val,
				values[0][1] * val, values[1][1] * val, values[2][1] * val, values[3][1] * val,
				values[0][2] * val, values[1][2] * val, values[2][2] * val, values[3][2] * val,
				values[0][3] * val, values[1][3] * val, values[2][3] * val, values[3][3] * val
			);
		}

		Matrix<1, 4, T, P> operator*(const Matrix<1, 4, T, P> &mat) const
		{
			Matrix<1, 4, T, P> res;
			res[0] = Col(
				values[0][0] * mat[0][0] + values[1][0] * mat[0][1] + values[2][0] * mat[0][2] + values[3][0] * mat[0][3],
				values[0][1] * mat[0][0] + values[1][1] * mat[0][1] + values[2][1] * mat[0][2] + values[3][1] * mat[0][3],
				values[0][2] * mat[0][0] + values[1][2] * mat[0][1] + values[2][2] * mat[0][2] + values[3][2] * mat[0][3],
				values[0][3] * mat[0][0] + values[1][3] * mat[0][1] + values[2][3] * mat[0][2] + values[3][3] * mat[0][3]
			);
			return res;
		}

		Matrix<2, 4, T, P> operator*(const Matrix<2, 4, T, P> &mat) const
		{
			Matrix<2, 4, T, P> res;
			res[0] = Col(
				values[0][0] * mat[0][0] + values[1][0] * mat[0][1] + values[2][0] * mat[0][2] + values[3][0] * mat[0][3],
				values[0][1] * mat[0][0] + values[1][1] * mat[0][1] + values[2][1] * mat[0][2] + values[3][1] * mat[0][3],
				values[0][2] * mat[0][0] + values[1][2] * mat[0][1] + values[2][2] * mat[0][2] + values[3][2] * mat[0][3],
				values[0][3] * mat[0][0] + values[1][3] * mat[0][1] + values[2][3] * mat[0][2] + values[3][3] * mat[0][3]
			);
			res[1] = Col(
				values[0][0] * mat[1][0] + values[1][0] * mat[1][1] + values[2][0] * mat[1][2] + values[3][0] * mat[1][3],
				values[0][1] * mat[1][0] + values[1][1] * mat[1][1] + values[2][1] * mat[1][2] + values[3][1] * mat[1][3],
				values[0][2] * mat[1][0] + values[1][2] * mat[1][1] + values[2][2] * mat[1][2] + values[3][2] * mat[1][3],
				values[0][3] * mat[1][0] + values[1][3] * mat[1][1] + values[2][3] * mat[1][2] + values[3][3] * mat[1][3]
			);
			return res;
		}

		Matrix<3, 4, T, P> operator*(const Matrix<3, 4, T, P> &mat) const
		{
			Matrix<3, 4, T, P> res;
			res[0] = Col(
				values[0][0] * mat[0][0] + values[1][0] * mat[0][1] + values[2][0] * mat[0][2] + values[3][0] * mat[0][3],
				values[0][1] * mat[0][0] + values[1][1] * mat[0][1] + values[2][1] * mat[0][2] + values[3][1] * mat[0][3],
				values[0][2] * mat[0][0] + values[1][2] * mat[0][1] + values[2][2] * mat[0][2] + values[3][2] * mat[0][3],
				values[0][3] * mat[0][0] + values[1][3] * mat[0][1] + values[2][3] * mat[0][2] + values[3][3] * mat[0][3]
			);
			res[1] = Col(
				values[0][0] * mat[1][0] + values[1][0] * mat[1][1] + values[2][0] * mat[1][2] + values[3][0] * mat[1][3],
				values[0][1] * mat[1][0] + values[1][1] * mat[1][1] + values[2][1] * mat[1][2] + values[3][1] * mat[1][3],
				values[0][2] * mat[1][0] + values[1][2] * mat[1][1] + values[2][2] * mat[1][2] + values[3][2] * mat[1][3],
				values[0][3] * mat[1][0] + values[1][3] * mat[1][1] + values[2][3] * mat[1][2] + values[3][3] * mat[1][3]
			);
			res[2] = Col(
				values[0][0] * mat[2][0] + values[1][0] * mat[2][1] + values[2][0] * mat[2][2] + values[3][0] * mat[2][3],
				values[0][1] * mat[2][0] + values[1][1] * mat[2][1] + values[2][1] * mat[2][2] + values[3][1] * mat[2][3],
				values[0][2] * mat[2][0] + values[1][2] * mat[2][1] + values[2][2] * mat[2][2] + values[3][2] * mat[2][3],
				values[0][3] * mat[2][0] + values[1][3] * mat[2][1] + values[2][3] * mat[2][2] + values[3][3] * mat[2][3]
			);
			return res;
		}

		template<unsigned int S>
		Matrix<S, 4, T, P> operator*(const Matrix<S, 4, T, P> &mat) const
		{
			Matrix<S, 4, T, P> res;

			for (int i = 0; i < S; i++)
			{
				res[i] = Col(
					values[0][0] * mat[i][0] + values[1][0] * mat[i][1] + values[2][0] * mat[i][2] + values[3][0] * mat[i][3],
					values[0][1] * mat[i][0] + values[1][1] * mat[i][1] + values[2][1] * mat[i][2] + values[3][1] * mat[i][3],
					values[0][2] * mat[i][0] + values[1][2] * mat[i][1] + values[2][2] * mat[i][2] + values[3][2] * mat[i][3],
					values[0][3] * mat[i][0] + values[1][3] * mat[i][1] + values[2][3] * mat[i][2] + values[3][3] * mat[i][3]
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
				values[0][0] * inv, values[1][0] * inv, values[2][0] * inv, values[3][0] * inv,
				values[0][1] * inv, values[1][1] * inv, values[2][1] * inv, values[3][1] * inv,
				values[0][2] * inv, values[1][2] * inv, values[2][2] * inv, values[3][2] * inv,
				values[0][3] * inv, values[1][3] * inv, values[2][3] * inv, values[3][3] * inv
			);
		}

		// -- Boolean operators --

		bool operator==(const Mat &mat) const
		{
			return values[0] == mat[0] && values[1] == mat[1] && values[2] == mat[2] && values[3] == mat[3];
		}

		bool operator!=(const Mat &mat) const
		{
			return values[0] != mat[0] || values[1] != mat[1] || values[2] != mat[2] || values[3] != mat[3];
		}

		// -- Convertion operators --

		template<typename U, typename Q>
		operator Matrix<4, 4, U, Q>() const
		{
			return Matrix<4, 4, U, Q>(
				(U)values[0][0], (U)values[1][0], (U)values[2][0], (U)values[3][0],
				(U)values[0][1], (U)values[1][1], (U)values[2][1], (U)values[3][1],
				(U)values[0][2], (U)values[1][2], (U)values[2][2], (U)values[3][2],
				(U)values[0][3], (U)values[1][3], (U)values[2][3], (U)values[3][3]
			);
		}

		// -- Stream operators --

		friend std::ostream &operator<<(std::ostream &ostream, const Mat &mat)
		{
			return ostream	<< mat[0][0] << ',' << mat[1][0] << ',' << mat[2][0] << ',' << mat[3][0] << '\n'
							<< mat[0][1] << ',' << mat[1][1] << ',' << mat[2][1] << ',' << mat[3][1] << '\n'
							<< mat[0][2] << ',' << mat[1][2] << ',' << mat[2][2] << ',' << mat[3][2] << '\n'
							<< mat[0][3] << ',' << mat[1][3] << ',' << mat[2][3] << ',' << mat[3][3] << '\n';
		}

		// -- Getters --

		Mat GetAdjugate() const
		{
			Mat res;

			// Col 1
			res[0][0] = values[0][0] * (
				  values[1][1] * (values[2][2] * values[3][3] - values[3][2] * values[2][3])
				- values[2][1] * (values[1][2] * values[3][3] - values[3][2] * values[1][3])
				+ values[3][1] * (values[1][2] * values[2][3] - values[2][2] * values[1][3])
			);
			res[0][1] = -values[1][0] * (
				  values[0][1] * (values[2][2] * values[3][3] - values[3][2] * values[2][3])
				- values[2][1] * (values[0][2] * values[3][3] - values[3][2] * values[0][3])
				+ values[3][1] * (values[0][2] * values[2][3] - values[2][2] * values[0][3])
			);
			res[0][2] = values[2][0] * (
				  values[0][1] * (values[1][2] * values[3][3] - values[3][2] * values[1][3])
				- values[1][1] * (values[0][2] * values[3][3] - values[3][2] * values[0][3])
				+ values[3][1] * (values[0][2] * values[1][3] - values[1][2] * values[0][3])
			);
			res[0][3] = -values[3][0] * (
				  values[0][1] * (values[1][2] * values[2][3] - values[2][2] * values[1][3])
				- values[1][1] * (values[0][2] * values[2][3] - values[2][2] * values[0][3])
				+ values[2][1] * (values[0][2] * values[1][3] - values[1][2] * values[0][3])
			);

			// Col 2
			res[1][0] = -values[0][1] * (
				  values[1][0] * (values[2][2] * values[3][3] - values[3][2] * values[2][3])
				- values[2][0] * (values[1][2] * values[3][3] - values[3][2] * values[1][3])
				+ values[3][0] * (values[1][2] * values[2][3] - values[2][2] * values[1][3])
			);
			res[1][1] = values[1][1] * (
				  values[0][0] * (values[2][2] * values[3][3] - values[3][2] * values[2][3])
				- values[2][0] * (values[0][2] * values[3][3] - values[3][2] * values[0][3])
				+ values[3][0] * (values[0][2] * values[2][3] - values[2][2] * values[0][3])
			);
			res[1][2] = -values[2][1] * (
				  values[0][0] * (values[1][2] * values[3][3] - values[3][2] * values[1][3])
				- values[1][0] * (values[0][2] * values[3][3] - values[3][2] * values[0][3])
				+ values[3][0] * (values[0][2] * values[1][3] - values[1][2] * values[0][3])
			);
			res[1][3] = values[3][1] * (
				  values[0][0] * (values[1][2] * values[2][3] - values[2][2] * values[1][3])
				- values[1][0] * (values[0][2] * values[2][3] - values[2][2] * values[0][3])
				+ values[2][0] * (values[0][2] * values[1][3] - values[1][2] * values[0][3])
			);

			// Col 3
			res[2][0] = values[0][2] * (
				  values[1][0] * (values[2][1] * values[3][3] - values[3][1] * values[2][3])
				- values[2][0] * (values[1][1] * values[3][3] - values[3][1] * values[1][3])
				+ values[3][0] * (values[1][1] * values[2][3] - values[2][1] * values[1][3])
			);
			res[2][1] = -values[1][2] * (
				  values[0][0] * (values[2][1] * values[3][3] - values[3][1] * values[2][3])
				- values[2][0] * (values[0][1] * values[3][3] - values[3][1] * values[0][3])
				+ values[3][0] * (values[0][1] * values[2][3] - values[2][1] * values[0][3])
			);
			res[2][2] = values[2][2] * (
				  values[0][0] * (values[1][1] * values[3][3] - values[3][1] * values[1][3])
				- values[1][0] * (values[0][1] * values[3][3] - values[3][1] * values[0][3])
				+ values[3][0] * (values[0][1] * values[1][3] - values[1][1] * values[0][3])
			);
			res[2][3] = -values[3][2] * (
				  values[0][0] * (values[2][1] * values[3][3] - values[3][1] * values[2][3])
				- values[1][0] * (values[1][1] * values[3][3] - values[3][1] * values[1][3])
				+ values[2][0] * (values[1][1] * values[2][3] - values[2][1] * values[1][3])
			);

			// Col 4
			res[3][0] = -values[0][3] * (
				  values[1][0] * (values[2][1] * values[3][2] - values[3][1] * values[2][2])
				- values[2][0] * (values[1][1] * values[3][2] - values[3][1] * values[1][2])
				+ values[3][0] * (values[1][1] * values[2][2] - values[2][1] * values[1][2])
			);
			res[3][1] = values[1][3] * (
				  values[0][0] * (values[2][1] * values[3][2] - values[3][1] * values[2][2])
				- values[2][0] * (values[0][1] * values[3][2] - values[3][1] * values[0][2])
				+ values[3][0] * (values[0][1] * values[2][2] - values[2][1] * values[0][2])
			);
			res[3][2] = -values[2][3] * (
				  values[0][0] * (values[1][1] * values[3][2] - values[3][1] * values[1][2])
				- values[1][0] * (values[0][1] * values[3][2] - values[3][1] * values[0][2])
				+ values[3][0] * (values[0][1] * values[1][2] - values[1][1] * values[0][2])
			);
			res[3][3] = values[3][3] * (
				  values[0][0] * (values[1][1] * values[2][2] - values[2][1] * values[1][2])
				- values[1][0] * (values[0][1] * values[2][2] - values[2][1] * values[0][2])
				+ values[2][0] * (values[0][1] * values[1][2] - values[1][1] * values[0][2])
			);

			return res;
		}

		Mat GetCofactor() const
		{
			Mat res;

			// Col 1
			res[0][0] = values[0][0] * (
				  values[1][1] * (values[2][2] * values[3][3] - values[3][2] * values[2][3])
				- values[2][1] * (values[1][2] * values[3][3] - values[3][2] * values[1][3])
				+ values[3][1] * (values[1][2] * values[2][3] - values[2][2] * values[1][3])
			);
			res[0][1] = -values[0][1] * (
				  values[1][0] * (values[2][2] * values[3][3] - values[3][2] * values[2][3])
				- values[2][0] * (values[1][2] * values[3][3] - values[3][2] * values[1][3])
				+ values[3][0] * (values[1][2] * values[2][3] - values[2][2] * values[1][3])
			);
			res[0][2] = values[0][2] * (
				  values[1][0] * (values[2][1] * values[3][3] - values[3][1] * values[2][3])
				- values[2][0] * (values[1][1] * values[3][3] - values[3][1] * values[1][3])
				+ values[3][0] * (values[1][1] * values[2][3] - values[2][1] * values[1][3])
			);
			res[0][3] = -values[0][3] * (
				  values[1][0] * (values[2][1] * values[3][2] - values[3][1] * values[2][2])
				- values[2][0] * (values[1][1] * values[3][2] - values[3][1] * values[1][2])
				+ values[3][0] * (values[1][1] * values[2][2] - values[2][1] * values[1][2])
			);

			// Col 2
			res[1][0] = -values[1][0] * (
				  values[0][1] * (values[2][2] * values[3][3] - values[3][2] * values[2][3])
				- values[2][1] * (values[0][2] * values[3][3] - values[3][2] * values[0][3])
				+ values[3][1] * (values[0][2] * values[2][3] - values[2][2] * values[0][3])
			);
			res[1][1] = values[1][1] * (
				  values[0][0] * (values[2][2] * values[3][3] - values[3][2] * values[2][3])
				- values[2][0] * (values[0][2] * values[3][3] - values[3][2] * values[0][3])
				+ values[3][0] * (values[0][2] * values[2][3] - values[2][2] * values[0][3])
			);
			res[1][2] = -values[1][2] * (
				  values[0][0] * (values[2][1] * values[3][3] - values[3][1] * values[2][3])
				- values[2][0] * (values[0][1] * values[3][3] - values[3][1] * values[0][3])
				+ values[3][0] * (values[0][1] * values[2][3] - values[2][1] * values[0][3])
			);
			res[1][3] = values[1][3] * (
				  values[0][0] * (values[2][1] * values[3][2] - values[3][1] * values[2][2])
				- values[2][0] * (values[0][1] * values[3][2] - values[3][1] * values[0][2])
				+ values[3][0] * (values[0][1] * values[2][2] - values[2][1] * values[0][2])
			);

			// Col 3
			res[2][0] = values[2][0] * (
				  values[0][1] * (values[1][2] * values[3][3] - values[3][2] * values[1][3])
				- values[1][1] * (values[0][2] * values[3][3] - values[3][2] * values[0][3])
				+ values[3][1] * (values[0][2] * values[1][3] - values[1][2] * values[0][3])
			);
			res[2][1] = -values[2][1] * (
				  values[0][0] * (values[1][2] * values[3][3] - values[3][2] * values[1][3])
				- values[1][0] * (values[0][2] * values[3][3] - values[3][2] * values[0][3])
				+ values[3][0] * (values[0][2] * values[1][3] - values[1][2] * values[0][3])
			);
			res[2][2] = values[2][2] * (
				  values[0][0] * (values[1][1] * values[3][3] - values[3][1] * values[1][3])
				- values[1][0] * (values[0][1] * values[3][3] - values[3][1] * values[0][3])
				+ values[3][0] * (values[0][1] * values[1][3] - values[1][1] * values[0][3])
			);
			res[2][3] = -values[2][3] * (
				  values[0][0] * (values[1][1] * values[3][2] - values[3][1] * values[1][2])
				- values[1][0] * (values[0][1] * values[3][2] - values[3][1] * values[0][2])
				+ values[3][0] * (values[0][1] * values[1][2] - values[1][1] * values[0][2])
			);

			// Col 4
			res[3][0] = -values[3][0] * (
				  values[0][1] * (values[1][2] * values[2][3] - values[2][2] * values[1][3])
				- values[1][1] * (values[0][2] * values[2][3] - values[2][2] * values[0][3])
				+ values[2][1] * (values[0][2] * values[1][3] - values[1][2] * values[0][3])
			);
			res[3][1] = values[3][1] * (
				  values[0][0] * (values[1][2] * values[2][3] - values[2][2] * values[1][3])
				- values[1][0] * (values[0][2] * values[2][3] - values[2][2] * values[0][3])
				+ values[2][0] * (values[0][2] * values[1][3] - values[1][2] * values[0][3])
			);
			res[3][2] = -values[3][2] * (
				  values[0][0] * (values[2][1] * values[3][3] - values[3][1] * values[2][3])
				- values[1][0] * (values[1][1] * values[3][3] - values[3][1] * values[1][3])
				+ values[2][0] * (values[1][1] * values[2][3] - values[2][1] * values[1][3])
			);
			res[3][3] = values[3][3] * (
				  values[0][0] * (values[1][1] * values[2][2] - values[2][1] * values[1][2])
				- values[1][0] * (values[0][1] * values[2][2] - values[2][1] * values[0][2])
				+ values[2][0] * (values[0][1] * values[1][2] - values[1][1] * values[0][2])
			);

			return res;
		}

		T GetDeterminant() const
		{
			return values[0][0] * (
				  values[1][1] * (values[2][2] * values[3][3] - values[3][2] * values[2][3])
				- values[2][1] * (values[1][2] * values[3][3] - values[3][2] * values[1][3])
				+ values[3][1] * (values[1][2] * values[2][3] - values[2][2] * values[1][3])
			)
			- values[1][0] * (
				  values[0][1] * (values[2][2] * values[3][3] - values[3][2] * values[2][3])
				- values[2][1] * (values[0][2] * values[3][3] - values[3][2] * values[0][3])
				+ values[3][1] * (values[0][2] * values[2][3] - values[2][2] * values[0][3])
			)
			+ values[2][0] * (
				  values[0][1] * (values[1][2] * values[3][3] - values[3][2] * values[1][3])
				- values[1][1] * (values[0][2] * values[3][3] - values[3][2] * values[0][3])
				+ values[3][1] * (values[0][2] * values[1][3] - values[1][2] * values[0][3])
			)
			- values[3][0] * (
				  values[0][1] * (values[1][2] * values[2][3] - values[2][2] * values[1][3])
				- values[1][1] * (values[0][2] * values[2][3] - values[2][2] * values[0][3])
				+ values[2][1] * (values[0][2] * values[1][3] - values[1][2] * values[0][3])
			);
		}

		Mat Inverted() const
		{
			return GetAdjugate() * (T(1) / GetDeterminant());
		}

		Mat Transposed() const
		{
			return Mat(
				values[0][0], values[0][1], values[0][2], values[0][3],
				values[1][0], values[1][1], values[1][2], values[1][3],
				values[2][0], values[2][1], values[2][2], values[2][3],
				values[3][0], values[3][1], values[3][2], values[3][3]
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

			values[0][3] += values[3][0];
			values[3][0] = values[0][3] - values[3][0];
			values[0][3] -= values[3][0];

			values[2][1] += values[1][2];
			values[1][2] = values[2][1] - values[1][2];
			values[2][1] -= values[1][2];

			values[3][1] += values[1][3];
			values[1][3] = values[3][1] - values[1][3];
			values[3][1] -= values[1][3];

			values[3][2] += values[2][3];
			values[2][3] = values[3][2] - values[2][3];
			values[3][2] -= values[2][3];
			return *this;
		}

		// -- Static getters --

		static inline Mat Identity()
		{
			return Mat(T(1));
		}
	};
}

