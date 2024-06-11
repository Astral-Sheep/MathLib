#include "Vector4.hpp"
#include "Matrix4x4.hpp"

namespace Math
{
	// -- Constructors --

	Matrix4x4::Matrix()
		: values{ Vector4(), Vector4(), Vector4(), Vector4() }
	{

	}

	Matrix4x4::Matrix(const float scalar)
		: values{
			Vector4(scalar, 0.0f, 0.0f, 0.0f),
			Vector4(0.0f, scalar, 0.0f, 0.0f),
			Vector4(0.0f, 0.0f, scalar, 0.0f),
			Vector4(0.0f, 0.0f, 0.0f, scalar)
		}
	{

	}

	Matrix4x4::Matrix(
		const float x0, const float x1, const float x2, const float x3,
		const float x4, const float x5, const float x6, const float x7,
		const float x8, const float x9, const float x10, const float x11,
		const float x12, const float x13, const float x14, const float x15
	)
		: values{
			Vector4(x0, x1, x2, x3),
			Vector4(x4, x5, x6, x7),
			Vector4(x8, x9, x10, x11),
			Vector4(x12, x13, x14, x15)
		}
	{

	}

	Matrix4x4::Matrix(const Vector4 columns[4])
		: values{ columns[0], columns[1], columns[2], columns[3] }
	{

	}

	Matrix4x4::Matrix(const float values[16])
		: values{
			Vector4(values[0], values[1], values[2], values[3]),
			Vector4(values[4], values[5], values[6], values[7]),
			Vector4(values[8], values[9], values[10], values[11]),
			Vector4(values[12], values[13], values[14], values[15])
		}
	{

	}

	Matrix4x4::Matrix(const Vector4 &x0, const Vector4 &x1, const Vector4 &x2, const Vector4 &x3)
		: values{ x0, x1, x2, x3 }
	{

	}

	Matrix4x4::Matrix(const Matrix4x4 &matrix)
		: values{
			matrix[0],
			matrix[1],
			matrix[2],
			matrix[3]
		}
	{

	}

	// -- Accesses --

	const Vector4 &Matrix4x4::operator[](const int index) const noexcept
	{
		return values[index];
	}

	Vector4 &Matrix4x4::operator[](const int index) noexcept
	{
		return values[index];
	}

	// -- Unary arithmetic operators --

	Matrix4x4 &Matrix4x4::operator+=(const Matrix4x4 &matrix)
	{
		values[0] += matrix[0];
		values[1] += matrix[1];
		values[2] += matrix[2];
		values[3] += matrix[3];
		return *this;
	}

	Matrix4x4 &Matrix4x4::operator-=(const Matrix4x4 &matrix)
	{
		values[0] -= matrix[0];
		values[1] -= matrix[1];
		values[2] -= matrix[2];
		values[3] -= matrix[3];
		return *this;
	}

	Matrix4x4 &Matrix4x4::operator*=(const Matrix4x4 &matrix)
	{
		Matrix4x4 temp(*this);
		values[0][0] = temp[0][0] * matrix[0][0] + temp[1][0] * matrix[0][1] + temp[2][0] * matrix[0][2] + temp[3][0] * matrix[0][3];
		values[0][1] = temp[0][1] * matrix[0][0] + temp[1][1] * matrix[0][1] + temp[2][1] * matrix[0][2] + temp[3][1] * matrix[0][3];
		values[0][2] = temp[0][2] * matrix[0][0] + temp[1][2] * matrix[0][1] + temp[2][2] * matrix[0][2] + temp[3][2] * matrix[0][3];
		values[0][3] = temp[0][3] * matrix[0][0] + temp[1][3] * matrix[0][1] + temp[2][3] * matrix[0][2] + temp[3][3] * matrix[0][3];

		values[1][0] = temp[0][0] * matrix[1][0] + temp[1][0] * matrix[1][1] + temp[2][0] * matrix[1][2] + temp[3][0] * matrix[1][3];
		values[1][1] = temp[0][1] * matrix[1][0] + temp[1][1] * matrix[1][1] + temp[2][1] * matrix[1][2] + temp[3][1] * matrix[1][3];
		values[1][2] = temp[0][2] * matrix[1][0] + temp[1][2] * matrix[1][1] + temp[2][2] * matrix[1][2] + temp[3][2] * matrix[1][3];
		values[1][3] = temp[0][3] * matrix[1][0] + temp[1][3] * matrix[1][1] + temp[2][3] * matrix[1][2] + temp[3][3] * matrix[1][3];

		values[2][0] = temp[0][0] * matrix[2][0] + temp[1][0] * matrix[2][1] + temp[2][0] * matrix[2][2] + temp[3][0] * matrix[2][3];
		values[2][1] = temp[0][1] * matrix[2][0] + temp[1][1] * matrix[2][1] + temp[2][1] * matrix[2][2] + temp[3][1] * matrix[2][3];
		values[2][2] = temp[0][2] * matrix[2][0] + temp[1][2] * matrix[2][1] + temp[2][2] * matrix[2][2] + temp[3][2] * matrix[2][3];
		values[2][3] = temp[0][3] * matrix[2][0] + temp[1][3] * matrix[2][1] + temp[2][3] * matrix[2][2] + temp[3][3] * matrix[2][3];

		values[3][0] = temp[0][0] * matrix[3][0] + temp[1][0] * matrix[3][1] + temp[2][0] * matrix[3][2] + temp[3][0] * matrix[3][3];
		values[3][1] = temp[0][1] * matrix[3][0] + temp[1][1] * matrix[3][1] + temp[2][1] * matrix[3][2] + temp[3][1] * matrix[3][3];
		values[3][2] = temp[0][2] * matrix[3][0] + temp[1][2] * matrix[3][1] + temp[2][2] * matrix[3][2] + temp[3][2] * matrix[3][3];
		values[3][3] = temp[0][3] * matrix[3][0] + temp[1][3] * matrix[3][1] + temp[2][3] * matrix[3][2] + temp[3][3] * matrix[3][3];
		return *this;
	}

	Matrix4x4 &Matrix4x4::operator*=(const float scalar)
	{
		values[0] *= scalar;
		values[1] *= scalar;
		values[2] *= scalar;
		values[3] *= scalar;
		return *this;
	}

	Matrix4x4 &Matrix4x4::operator/=(const Matrix4x4 &matrix)
	{
		return *this *= matrix.Inverted();
	}

	Matrix4x4 &Matrix4x4::operator/=(const float scalar)
	{
		values[0] /= scalar;
		values[1] /= scalar;
		values[2] /= scalar;
		values[3] /= scalar;
		return *this;
	}

	// -- Unary operators --

	Matrix4x4 Matrix4x4::operator-() const
	{
		return Matrix4x4(
			-values[0],
			-values[1],
			-values[2],
			-values[3]
		);
	}

	// -- Binary operators --

	Matrix4x4 Matrix4x4::operator+(const Matrix4x4 &matrix) const
	{
		return Matrix4x4(
			values[0] + matrix[0],
			values[1] + matrix[1],
			values[2] + matrix[2],
			values[3] + matrix[3]
		);
	}

	Matrix4x4 Matrix4x4::operator-(const Matrix4x4 &matrix) const
	{
		return Matrix4x4(
			values[0] - matrix[0],
			values[1] - matrix[1],
			values[2] - matrix[2],
			values[3] - matrix[3]
		);
	}

	Matrix4x4 Matrix4x4::operator*(const Matrix4x4 &matrix) const
	{
		return Matrix4x4(
			Vector4(
				values[0][0] * matrix[0][0] + values[1][0] * matrix[0][1] + values[2][0] * matrix[0][2] + values[3][0] * matrix[0][3],
				values[0][1] * matrix[0][0] + values[1][1] * matrix[0][1] + values[2][1] * matrix[0][2] + values[3][1] * matrix[0][3],
				values[0][2] * matrix[0][0] + values[1][2] * matrix[0][1] + values[2][2] * matrix[0][2] + values[3][2] * matrix[0][3],
				values[0][3] * matrix[0][0] + values[1][3] * matrix[0][1] + values[2][3] * matrix[0][2] + values[3][3] * matrix[0][3]
			),
			Vector4(
				values[0][0] * matrix[1][0] + values[1][0] * matrix[1][1] + values[2][0] * matrix[1][2] + values[3][0] * matrix[1][3],
				values[0][1] * matrix[1][0] + values[1][1] * matrix[1][1] + values[2][1] * matrix[1][2] + values[3][1] * matrix[1][3],
				values[0][2] * matrix[1][0] + values[1][2] * matrix[1][1] + values[2][2] * matrix[1][2] + values[3][2] * matrix[1][3],
				values[0][3] * matrix[1][0] + values[1][3] * matrix[1][1] + values[2][3] * matrix[1][2] + values[3][3] * matrix[1][3]
			),
			Vector4(
				values[0][0] * matrix[2][0] + values[1][0] * matrix[2][1] + values[2][0] * matrix[2][2] + values[3][0] * matrix[2][3],
				values[0][1] * matrix[2][0] + values[1][1] * matrix[2][1] + values[2][1] * matrix[2][2] + values[3][1] * matrix[2][3],
				values[0][2] * matrix[2][0] + values[1][2] * matrix[2][1] + values[2][2] * matrix[2][2] + values[3][2] * matrix[2][3],
				values[0][3] * matrix[2][0] + values[1][3] * matrix[2][1] + values[2][3] * matrix[2][2] + values[3][3] * matrix[2][3]
			),
			Vector4(
				values[0][0] * matrix[3][0] + values[1][0] * matrix[3][1] + values[2][0] * matrix[3][2] + values[3][0] * matrix[3][3],
				values[0][1] * matrix[3][0] + values[1][1] * matrix[3][1] + values[2][1] * matrix[3][2] + values[3][1] * matrix[3][3],
				values[0][2] * matrix[3][0] + values[1][2] * matrix[3][1] + values[2][2] * matrix[3][2] + values[3][2] * matrix[3][3],
				values[0][3] * matrix[3][0] + values[1][3] * matrix[3][1] + values[2][3] * matrix[3][2] + values[3][3] * matrix[3][3]
			)
		);
	}

	Vector4 Matrix4x4::operator*(const Vector4 &vector) const
	{
		return Vector4(
			values[0][0] * vector.x + values[1][0] * vector.y + values[2][0] * vector.z + values[3][0] * vector.w,
			values[0][1] * vector.x + values[1][1] * vector.y + values[2][1] * vector.z + values[3][1] * vector.w,
			values[0][2] * vector.x + values[1][2] * vector.y + values[2][2] * vector.z + values[3][2] * vector.w,
			values[0][3] * vector.x + values[1][3] * vector.y + values[2][3] * vector.z + values[3][3] * vector.w
		);
	}

	Matrix1x4 Matrix4x4::operator*(const Matrix1x4 &matrix) const
	{
		Matrix1x4 res;
		res[0] = Vector4(
			values[0][0] * matrix[0][0] + values[1][0] * matrix[0][1] + values[2][0] * matrix[0][2] + values[3][0] * matrix[0][3],
			values[0][1] * matrix[0][0] + values[1][1] * matrix[0][1] + values[2][1] * matrix[0][2] + values[3][1] * matrix[0][3],
			values[0][2] * matrix[0][0] + values[1][2] * matrix[0][1] + values[2][2] * matrix[0][2] + values[3][2] * matrix[0][3],
			values[0][3] * matrix[0][0] + values[1][3] * matrix[0][1] + values[2][3] * matrix[0][2] + values[3][3] * matrix[0][3]
		);
		return res;
	}

	Matrix2x4 Matrix4x4::operator*(const Matrix2x4 &matrix) const
	{
		Matrix2x4 res;
		res[0] = Vector4(
			values[0][0] * matrix[0][0] + values[1][0] * matrix[0][1] + values[2][0] * matrix[0][2] + values[3][0] * matrix[0][3],
			values[0][1] * matrix[0][0] + values[1][1] * matrix[0][1] + values[2][1] * matrix[0][2] + values[3][1] * matrix[0][3],
			values[0][2] * matrix[0][0] + values[1][2] * matrix[0][1] + values[2][2] * matrix[0][2] + values[3][2] * matrix[0][3],
			values[0][3] * matrix[0][0] + values[1][3] * matrix[0][1] + values[2][3] * matrix[0][2] + values[3][3] * matrix[0][3]
		);
		res[1] = Vector4(
			values[0][0] * matrix[1][0] + values[1][0] * matrix[1][1] + values[2][0] * matrix[1][2] + values[3][0] * matrix[1][3],
			values[0][1] * matrix[1][0] + values[1][1] * matrix[1][1] + values[2][1] * matrix[1][2] + values[3][1] * matrix[1][3],
			values[0][2] * matrix[1][0] + values[1][2] * matrix[1][1] + values[2][2] * matrix[1][2] + values[3][2] * matrix[1][3],
			values[0][3] * matrix[1][0] + values[1][3] * matrix[1][1] + values[2][3] * matrix[1][2] + values[3][3] * matrix[1][3]
		);
		return res;
	}

	Matrix3x4 Matrix4x4::operator*(const Matrix3x4 &matrix) const
	{
		Matrix3x4 res;
		res[0] = Vector4(
			values[0][0] * matrix[0][0] + values[1][0] * matrix[0][1] + values[2][0] * matrix[0][2] + values[3][0] * matrix[0][3],
			values[0][1] * matrix[0][0] + values[1][1] * matrix[0][1] + values[2][1] * matrix[0][2] + values[3][1] * matrix[0][3],
			values[0][2] * matrix[0][0] + values[1][2] * matrix[0][1] + values[2][2] * matrix[0][2] + values[3][2] * matrix[0][3],
			values[0][3] * matrix[0][0] + values[1][3] * matrix[0][1] + values[2][3] * matrix[0][2] + values[3][3] * matrix[0][3]
		);
		res[1] = Vector4(
			values[0][0] * matrix[1][0] + values[1][0] * matrix[1][1] + values[2][0] * matrix[1][2] + values[3][0] * matrix[1][3],
			values[0][1] * matrix[1][0] + values[1][1] * matrix[1][1] + values[2][1] * matrix[1][2] + values[3][1] * matrix[1][3],
			values[0][2] * matrix[1][0] + values[1][2] * matrix[1][1] + values[2][2] * matrix[1][2] + values[3][2] * matrix[1][3],
			values[0][3] * matrix[1][0] + values[1][3] * matrix[1][1] + values[2][3] * matrix[1][2] + values[3][3] * matrix[1][3]
		);
		res[2] = Vector4(
			values[0][0] * matrix[2][0] + values[1][0] * matrix[2][1] + values[2][0] * matrix[2][2] + values[3][0] * matrix[2][3],
			values[0][1] * matrix[2][0] + values[1][1] * matrix[2][1] + values[2][1] * matrix[2][2] + values[3][1] * matrix[2][3],
			values[0][2] * matrix[2][0] + values[1][2] * matrix[2][1] + values[2][2] * matrix[2][2] + values[3][2] * matrix[2][3],
			values[0][3] * matrix[2][0] + values[1][3] * matrix[2][1] + values[2][3] * matrix[2][2] + values[3][3] * matrix[2][3]
		);
		return res;
	}

	template<unsigned int S>
	Matrix<S, 4, float> Matrix4x4::operator*(const Matrix<S, 4, float> &matrix) const
	{
		Matrix<S, 4, float> res;

		for (int i = 0; i < S; i++)
		{
			res[i] = Vector4(
				values[0][0] * matrix[i][0] + values[1][0] * matrix[i][1] + values[2][0] * matrix[i][2] + values[3][0] * matrix[i][3],
				values[0][1] * matrix[i][0] + values[1][1] * matrix[i][1] + values[2][1] * matrix[i][2] + values[3][1] * matrix[i][3],
				values[0][2] * matrix[i][0] + values[1][2] * matrix[i][1] + values[2][2] * matrix[i][2] + values[3][2] * matrix[i][3],
				values[0][3] * matrix[i][0] + values[1][3] * matrix[i][1] + values[2][3] * matrix[i][2] + values[3][3] * matrix[i][3]
			);
		}

		return res;
	}

	Matrix4x4 Matrix4x4::operator*(const float scalar) const
	{
		return Matrix4x4(
			values[0] * scalar,
			values[1] * scalar,
			values[2] * scalar,
			values[3] * scalar
		);
	}

	Matrix4x4 Matrix4x4::operator/(const Matrix4x4 &matrix) const
	{
		return *this * matrix.Inverted();
	}

	Matrix4x4 Matrix4x4::operator/(const float scalar) const
	{
		return Matrix4x4(
			values[0] / scalar,
			values[1] / scalar,
			values[2] / scalar,
			values[3] / scalar
		);
	}

	// -- Boolean operators --

	bool Matrix4x4::operator==(const Matrix4x4 &matrix) const
	{
		return !(*this != matrix);
	}

	bool Matrix4x4::operator!=(const Matrix4x4 &matrix) const
	{
		return values[0] != matrix[0] || values[1] != matrix[1] || values[2] != matrix[2] || values[3] != matrix[3];
	}

	// -- Stream operators --

	std::ostream &operator<<(std::ostream &ostream, const Matrix4x4 &matrix)
	{
		return ostream	<< matrix[0][0] << ',' << matrix[1][0] << ',' << matrix[2][0] << ',' << matrix[3][0] << '\n'
						<< matrix[0][1] << ',' << matrix[1][1] << ',' << matrix[2][1] << ',' << matrix[3][1] << '\n'
						<< matrix[0][2] << ',' << matrix[1][2] << ',' << matrix[2][2] << ',' << matrix[3][2] << '\n'
						<< matrix[0][3] << ',' << matrix[1][3] << ',' << matrix[2][3] << ',' << matrix[3][3] << '\n';
	}

	// -- Getters --

	Matrix4x4 Matrix4x4::GetAdjugate() const
	{
		Matrix4x4 res;

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

	Matrix4x4 Matrix4x4::GetCofactor() const
	{
		Matrix4x4 res;

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

	float Matrix4x4::GetDeterminant() const
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

	Matrix4x4 Matrix4x4::Inverted() const
	{
		return GetAdjugate() / GetDeterminant();
	}

	Matrix4x4 Matrix4x4::Transposed() const
	{
		return Matrix4x4(
			values[0][0], values[1][0], values[2][0], values[3][0],
			values[0][1], values[1][1], values[2][1], values[3][1],
			values[0][2], values[1][2], values[2][2], values[3][2],
			values[0][3], values[1][3], values[2][3], values[3][3]
		);
	}

	// -- Transformations --

	Matrix4x4 &Matrix4x4::Invert()
	{
		*this = GetAdjugate() / GetDeterminant();
		return *this;
	}

	Matrix4x4 &Matrix4x4::Transpose()
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
}

