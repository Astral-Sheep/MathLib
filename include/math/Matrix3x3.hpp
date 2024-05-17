#pragma once

#include "Matrix.hpp"

namespace Math
{
	template<>
	struct Matrix<3, 3, float>
	{
		MATRIX_BODY(3, 3, float, float)
		SQUARE_MAT_FUNCTIONS(3, float)

	public:
		// -- Constructors --
		Matrix(
			const float x0, const float x1, const float x2,
			const float x3, const float x4, const float x5,
			const float x6, const float x7, const float x8
		);
		Matrix(const Vector3 &x0, const Vector3 &x1, const Vector3 &x2);

		// -- Binary operators --
		Vector3 operator*(const Vector3 &matrix) const;
		Matrix1x3 operator*(const Matrix1x3 &matrix) const;
		Matrix2x3 operator*(const Matrix2x3 &matrix) const;
		Matrix4x3 operator*(const Matrix4x3 &matrix) const;

		// -- Static getters --
		static inline Matrix3x3 Identity()
		{
			return Matrix3x3(
				1, 0, 0,
				0, 1, 0,
				0, 0, 1
			);
		}
	};
}

