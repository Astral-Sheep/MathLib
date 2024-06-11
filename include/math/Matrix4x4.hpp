#pragma once

#include "Core.h"
#include "Matrix.hpp"

namespace Math
{
	template<>
	struct MATHLIB Matrix<4, 4, float>
	{
		MATRIX_BODY(4, 4, float, float)
		SQUARE_MAT_FUNCTIONS(4, float)

	public:
		// -- Constructors --
		Matrix(
			const float x0, const float x1, const float x2, const float x3,
			const float x4, const float x5, const float x6, const float x7,
			const float x8, const float x9, const float x10, const float x11,
			const float x12, const float x13, const float x14, const float x15
		);
		Matrix(const Vector4 &x0, const Vector4 &x1, const Vector4 &x2, const Vector4 &x3);

		// -- Binary operators --
		Vector4 operator*(const Vector4 &vector) const;
		Matrix1x4 operator*(const Matrix1x4 &matrix) const;
		Matrix2x4 operator*(const Matrix2x4 &matrix) const;
		Matrix3x4 operator*(const Matrix3x4 &matrix) const;

		// -- Static getters --
		static inline Matrix4x4 Identity()
		{
			return Matrix4x4(
				1.0f, 0.0f, 0.0f, 0.0f,
				0.0f, 1.0f, 0.0f, 0.0f,
				0.0f, 0.0f, 1.0f, 0.0f,
				0.0f, 0.0f, 0.0f, 1.0f
			);
		}
	};
}

