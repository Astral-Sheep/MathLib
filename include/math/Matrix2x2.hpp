#pragma once

#include "Matrix.hpp"

namespace Math
{
	template<>
	struct Matrix<2, 2, float>
	{
		MATRIX_BODY(2, 2, float, float)
		SQUARE_MAT_FUNCTIONS(2, float);

	public:
		// -- Constructors --
		Matrix(const float x0, const float x1, const float x2, const float x3);
		Matrix(const Vector2 &x0, const Vector2 &x1);

		// -- Binary operators --
		Vector2 operator*(const Vector2 &vector) const;
		Matrix1x2 operator*(const Matrix1x2 &matrix) const;
		Matrix3x2 operator*(const Matrix3x2 &matrix) const;
		Matrix4x2 operator*(const Matrix4x2 &matrix) const;

		// -- Static getters --
		static inline Matrix2x2 Identity()
		{
			return Matrix2x2(
				1, 0,
				0, 1
			);
		}
	};
}

