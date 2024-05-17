#include "Vector2.hpp"
#include "Matrix2x2.hpp"

namespace Math
{
	// -- Constructors --

	Matrix2x2::Matrix()
		: values{ Vector2(), Vector2() }
	{

	}

	Matrix2x2::Matrix(const float value)
		: values{ Vector2(value, 0.0f), Vector2(0.0f, value) }
	{

	}

	Matrix2x2::Matrix(const float x0, const float x1, const float x2, const float x3)
		: values{ Vector2(x0, x1), Vector2(x2, x3) }
	{

	}

	Matrix2x2::Matrix(const Vector2 &x0, const Vector2 &x1)
		: values{ x0, x1 }
	{

	}

	Matrix2x2::Matrix(const Vector2 values[2])
		: values{ values[0], values[1] }
	{

	}

	Matrix2x2::Matrix(const float values[4])
		: values{ Vector2(values[0], values[1]), Vector2(values[2], values[3]) }
	{

	}

	Matrix2x2::Matrix(const Matrix2x2 &matrix)
		: values{ matrix[0], matrix[1] }
	{

	}

	// -- Accesses --

	const Vector2 &Matrix2x2::operator[](const int index) const noexcept
	{
		return values[index];
	}

	Vector2 &Matrix2x2::operator[](const int index) noexcept
	{
		return values[index];
	}

	// -- Unary arithmetic operators --

	Matrix2x2 &Matrix2x2::operator+=(const Matrix2x2 &matrix)
	{
		values[0][0] += matrix[0][0];
		values[0][1] += matrix[0][1];
		values[1][0] += matrix[1][0];
		values[1][1] += matrix[1][1];
		return *this;
	}

	Matrix2x2 &Matrix2x2::operator-=(const Matrix2x2 &matrix)
	{
		values[0][0] -= matrix[0][0];
		values[0][1] -= matrix[0][1];
		values[1][0] -= matrix[1][0];
		values[1][1] -= matrix[1][1];
		return *this;
	}

	Matrix2x2 &Matrix2x2::operator*=(const Matrix2x2 &matrix)
	{
		Matrix2x2 temp(*this);
		values[0][0] = temp[0][0] * matrix[0][0] + temp[1][0] * matrix[0][1];
		values[0][1] = temp[0][1] * matrix[0][0] + temp[1][1] * matrix[0][1];
		values[1][0] = temp[0][0] * matrix[1][0] + temp[1][0] * matrix[1][1];
		values[1][1] = temp[0][1] * matrix[1][0] + temp[1][1] * matrix[1][1];
		return *this;
	}

	Matrix2x2 &Matrix2x2::operator*=(const float scalar)
	{
		values[0][0] *= scalar;
		values[0][1] *= scalar;
		values[1][0] *= scalar;
		values[1][1] *= scalar;
		return *this;
	}

	Matrix2x2 &Matrix2x2::operator/=(const Matrix2x2 &matrix)
	{
		return *this *= matrix.Inverted();
	}

	Matrix2x2 &Matrix2x2::operator/=(const float scalar)
	{
		values[0][0] /= scalar;
		values[0][1] /= scalar;
		values[1][0] /= scalar;
		values[1][1] /= scalar;
		return *this;
	}

	// -- Unary operators --

	Matrix2x2 Matrix2x2::operator-() const
	{
		return Matrix2x2(-values[0], -values[1]);
	}

	// -- Binary operators --

	Matrix2x2 Matrix2x2::operator+(const Matrix2x2 &matrix) const
	{
		return Matrix2x2(values[0] + matrix[0], values[1] + matrix[1]);
	}

	Matrix2x2 Matrix2x2::operator-(const Matrix2x2 &matrix) const
	{
		return Matrix2x2(values[0] - matrix[0], values[1] - matrix[1]);
	}

	Matrix2x2 Matrix2x2::operator*(const Matrix2x2 &matrix) const
	{
		return Matrix2x2(
			values[0][0] * matrix[0][0] + values[1][0] * matrix[0][1],
			values[0][1] * matrix[0][0] + values[1][1] * matrix[0][1],
			values[0][0] * matrix[1][0] + values[1][0] * matrix[1][1],
			values[0][1] * matrix[1][0] + values[1][1] * matrix[1][1]
		);
	}

	Vector2 Matrix2x2::operator*(const Vector2 &vector) const
	{
		return Vector2(
			values[0][0] * vector[0] + values[1][0] * vector[1],
			values[0][1] * vector[0] + values[1][1] * vector[1]
		);
	}

	Matrix1x2 Matrix2x2::operator*(const Matrix1x2 &matrix) const
	{
		Vector2 columns[1] = {
			Vector2(
				values[0][0] * matrix[0][0] + values[1][0] * matrix[0][1],
				values[0][1] * matrix[0][0] + values[1][1] * matrix[0][1]
			)
		};
		return Matrix1x2(columns);
	}

	Matrix3x2 Matrix2x2::operator*(const Matrix3x2 &matrix) const
	{
		Vector2 columns[3] = {
			Vector2(
				values[0][0] * matrix[0][0] + values[1][0] * matrix[0][1],
				values[0][1] * matrix[0][0] + values[1][1] * matrix[0][1]
			),
			Vector2(
				values[0][0] * matrix[1][0] + values[1][0] * matrix[1][1],
				values[0][1] * matrix[1][0] + values[1][1] * matrix[1][1]
			),
			Vector2(
				values[0][0] * matrix[2][0] + values[1][0] * matrix[2][1],
				values[0][1] * matrix[2][0] + values[1][1] * matrix[2][1]
			)
		};
		return Matrix3x2(columns);
	}

	Matrix4x2 Matrix2x2::operator*(const Matrix4x2 &matrix) const
	{
		Vector2 columns[4] = {
			Vector2(
				values[0][0] * matrix[0][0] + values[1][0] * matrix[0][1],
				values[0][1] * matrix[0][0] + values[1][1] * matrix[0][1]
			),
			Vector2(
				values[0][0] * matrix[1][0] + values[1][0] * matrix[1][1],
				values[0][1] * matrix[1][0] + values[1][1] * matrix[1][1]
			),
			Vector2(
				values[0][0] * matrix[2][0] + values[1][0] * matrix[2][1],
				values[0][1] * matrix[2][0] + values[1][1] * matrix[2][1]
			),
			Vector2(
				values[0][0] * matrix[3][0] + values[1][0] * matrix[3][1],
				values[0][1] * matrix[3][0] + values[1][1] * matrix[3][1]
			)
		};

		return Matrix4x2(columns);
	}

	template<unsigned int S>
	Matrix<S, 2, float> Matrix2x2::operator*(const Matrix<S, 2, float> &matrix) const
	{
		Vector2 columns[S];

		for (int i = 0; i < S; i++)
		{
			columns[i] = Vector2(
				values[0][0] * matrix[i][0] + values[1][0] * matrix[i][1],
				values[0][1] * matrix[i][0] + values[1][1] * matrix[i][1]
			);
		}

		return Matrix<S, 2, float>(columns);
	}

	Matrix2x2 Matrix2x2::operator*(const float scalar) const
	{
		return Matrix2x2(
			values[0] * scalar,
			values[1] * scalar
		);
	}

	Matrix2x2 Matrix2x2::operator/(const Matrix2x2 &matrix) const
	{
		return *this * matrix.Inverted();
	}

	Matrix2x2 Matrix2x2::operator/(const float scalar) const
	{
		return Matrix2x2(
			values[0] / scalar,
			values[1] / scalar
		);
	}

	// -- Boolean operators --

	bool Matrix2x2::operator==(const Matrix2x2 &matrix) const
	{
		return !(*this != matrix);
	}

	bool Matrix2x2::operator!=(const Matrix2x2 &matrix) const
	{
		return values[0][0] != matrix[0][0] || values[0][1] != matrix[0][1]
			|| values[1][0] != matrix[1][0] || values[1][1] != matrix[1][1];
	}

	// -- Stream operators --

	std::ostream &operator<<(std::ostream &ostream, const Matrix2x2 &matrix)
	{
		return ostream << matrix[0][0] << ',' << matrix[1][0] << '\n'
					   << matrix[0][1] << ',' << matrix[1][1] << '\n';
	}

	// -- Getters --

	Matrix2x2 Matrix2x2::GetAdjugate() const
	{
		return Matrix2x2(values[1][1], -values[0][1], -values[1][0], values[0][0]);
	}

	Matrix2x2 Matrix2x2::GetCofactor() const
	{
		return Matrix2x2(values[1][1], -values[1][0], -values[0][1], values[0][0]);
	}

	float Matrix2x2::GetDeterminant() const
	{
		return values[0][0] * values[1][1] - values[1][0] * values[0][1];
	}

	Matrix2x2 Matrix2x2::Inverted() const
	{
		return GetAdjugate() * (1.0f / GetDeterminant());
	}

	Matrix2x2 Matrix2x2::Transposed() const
	{
		return Matrix2x2(values[0][0], values[1][0], values[0][1], values[1][1]);
	}

	// -- Transformations --

	Matrix2x2 &Matrix2x2::Invert()
	{
		*this = GetAdjugate() / GetDeterminant();
		return *this;
	}

	Matrix2x2 &Matrix2x2::Transpose()
	{
		values[0][1] += values[1][0];
		values[1][0] = values[0][1] - values[1][0];
		values[0][1] -= values[1][0];
		return *this;
	}
}

