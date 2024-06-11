#pragma once

#include "Matrix.hpp"
#include "Vector.hpp"
#include <cmath>

namespace Math
{
	template<typename T>
	Matrix<4, 4, T> Translate(const Matrix<4, 4, T> &mat, const Vector<3, T, float> &vec)
	{
		Matrix<4, 4, T> res(mat);
		res[3] = mat[0] * vec[0] + mat[1] * vec[1] + mat[2] * vec[2] + mat[3];
		return res;
	}

	template<typename T>
	Matrix<4, 4, T> Rotate(const Matrix<4, 4, T> &mat, const T angle, const Vector<3, T, float> &vec)
	{
		const T cos = std::cos(angle);
		const T sin = std::sin(angle);

		Vector<3, T, float> axis = vec.Normalized();
		Vector<3, T, float> temp(axis * (T(1) - cos));

		Matrix<4, 4, T> rotate;
		rotate[0][0] = cos + temp[0] * axis[0];
		rotate[0][1] = temp[0] * axis[1] + sin * axis[2];
		rotate[0][2] = temp[0] * axis[2] - sin * axis[1];

		rotate[1][0] = temp[1] * axis[0] - sin * axis[2];
		rotate[1][1] = cos + temp[1] * axis[1];
		rotate[1][2] = temp[1] * axis[2] + sin * axis[0];

		rotate[2][0] = temp[2] * axis[0] + sin * axis[1];
		rotate[2][1] = temp[2] * axis[1] - sin * axis[0];
		rotate[2][2] = cos + temp[2] * axis[2];

		Matrix<4, 4, T> res;
		res[0] = mat[0] * rotate[0][0] + mat[1] * rotate[0][1] + mat[2] * rotate[0][2];
		res[1] = mat[0] * rotate[1][0] + mat[1] * rotate[1][1] + mat[2] * rotate[1][2];
		res[2] = mat[0] * rotate[2][0] + mat[1] * rotate[2][1] + mat[2] * rotate[2][2];
		res[3] = mat[3];
		return res;
	}

	template<typename T>
	Matrix<4, 4, T> Scale(const Matrix<4, 4, T> &mat, const Vector<3, T, float> &vec)
	{
		Matrix<4, 4, T> res;
		res[0] = mat[0] * vec[0];
		res[1] = mat[1] * vec[1];
		res[2] = mat[2] * vec[2];
		res[3] = mat[3];
		return res;
	}
}

