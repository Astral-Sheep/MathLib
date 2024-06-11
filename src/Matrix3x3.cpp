#include "Vector3.hpp"
#include "Matrix3x3.hpp"

namespace Math
{
	// -- Constructors --

	Matrix3x3::Matrix()
		: values{ Vector3(), Vector3(), Vector3() }
	{

	}

	Matrix3x3::Matrix(const float scalar)
		: values{
			Vector3(scalar, 0.0f, 0.0f),
			Vector3(0.0f, scalar, 0.0f),
			Vector3(0.0f, 0.0f, scalar)
		}
	{

	}

	Matrix3x3::Matrix(
		const float x0, const float x1, const float x2,
		const float x3, const float x4, const float x5,
		const float x6, const float x7, const float x8
	) : values{
			Vector3(x0, x1, x2),
			Vector3(x3, x4, x5),
			Vector3(x6, x7, x8)
		}
	{

	}

	Matrix3x3::Matrix(const Vector3 columns[3])
		: values{ columns[0], columns[1], columns[2] }
	{

	}

	Matrix3x3::Matrix(const float columns[9])
		: values{
			Vector3(columns[0], columns[1], columns[2]),
			Vector3(columns[3], columns[4], columns[5]),
			Vector3(columns[6], columns[7], columns[8]),
		}
	{

	}

	Matrix3x3::Matrix(const Vector3 &x0, const Vector3 &x1, const Vector3 &x2)
		: values{ x0, x1, x2 }
	{

	}

	Matrix3x3::Matrix(const Matrix3x3 &matrix)
		: values{ matrix[0], matrix[1], matrix[2] }
	{

	}

	// -- Accessors --

	const Vector3 &Matrix3x3::operator[](const int index) const noexcept
	{
		return values[index];
	}

	Vector3 &Matrix3x3::operator[](const int index) noexcept
	{
		return values[index];
	}

	// -- Unary arithmetic operators --

	Matrix3x3 &Matrix3x3::operator+=(const Matrix3x3 &matrix)
	{
		values[0] += matrix[0];
		values[1] += matrix[1];
		values[2] += matrix[2];
		return *this;
	}

	Matrix3x3 &Matrix3x3::operator-=(const Matrix3x3 &matrix)
	{
		values[0] -= matrix[0];
		values[1] -= matrix[1];
		values[2] -= matrix[2];
		return *this;
	}

	Matrix3x3 &Matrix3x3::operator*=(const Matrix3x3 &matrix)
	{
		Matrix3x3 temp(*this);

		values[0][0] = temp[0][0] * matrix[0][0] + temp[1][0] * matrix[0][1] + temp[2][0] * matrix[0][2];
		values[0][1] = temp[0][1] * matrix[0][0] + temp[1][1] * matrix[0][1] + temp[2][1] * matrix[0][2];
		values[0][2] = temp[0][2] * matrix[0][0] + temp[1][2] * matrix[0][1] + temp[2][2] * matrix[0][2];

		values[1][0] = temp[0][0] * matrix[1][0] + temp[1][0] * matrix[1][1] + temp[2][0] * matrix[1][2];
		values[1][1] = temp[0][1] * matrix[1][0] + temp[1][1] * matrix[1][1] + temp[2][1] * matrix[1][2];
		values[1][2] = temp[0][2] * matrix[1][0] + temp[1][2] * matrix[1][1] + temp[2][2] * matrix[1][2];

		values[2][0] = temp[0][0] * matrix[2][0] + temp[1][0] * matrix[2][1] + temp[2][0] * matrix[2][2];
		values[2][1] = temp[0][1] * matrix[2][0] + temp[1][1] * matrix[2][1] + temp[2][1] * matrix[2][2];
		values[2][2] = temp[0][2] * matrix[2][0] + temp[1][2] * matrix[2][1] + temp[2][2] * matrix[2][2];
		return *this;
	}

	Matrix3x3 &Matrix3x3::operator*=(const float scalar)
	{
		values[0] *= scalar;
		values[1] *= scalar;
		values[2] *= scalar;
		return *this;
	}

	Matrix3x3 &Matrix3x3::operator/=(const Matrix3x3 &matrix)
	{
		return *this *= matrix.Inverted();
	}

	Matrix3x3 &Matrix3x3::operator/=(const float scalar)
	{
		values[0] /= scalar;
		values[1] /= scalar;
		values[2] /= scalar;
		return *this;
	}

	// -- Unary operators --

	Matrix3x3 Matrix3x3::operator-() const
	{
		return Matrix3x3(-values[0], -values[1], -values[2]);
	}

	// -- Binary operators --

	Matrix3x3 Matrix3x3::operator+(const Matrix3x3 &matrix) const
	{
		return Matrix3x3(
			values[0] + matrix[0],
			values[1] + matrix[1],
			values[2] + matrix[2]
		);
	}

	Matrix3x3 Matrix3x3::operator-(const Matrix3x3 &matrix) const
	{
		return Matrix3x3(
			values[0] - matrix[0],
			values[1] - matrix[1],
			values[2] - matrix[2]
		);
	}

	Matrix3x3 Matrix3x3::operator*(const Matrix3x3 &matrix) const
	{
		return Matrix3x3(
			values[0][0] * matrix[0][0] + values[1][0] * matrix[0][1] + values[2][0] * matrix[0][2],
			values[0][1] * matrix[0][0] + values[1][1] * matrix[0][1] + values[2][1] * matrix[0][2],
			values[0][2] * matrix[0][0] + values[1][2] * matrix[0][1] + values[2][2] * matrix[0][2],

			values[0][0] * matrix[1][0] + values[1][0] * matrix[1][1] + values[2][0] * matrix[1][2],
			values[0][1] * matrix[1][0] + values[1][1] * matrix[1][1] + values[2][1] * matrix[1][2],
			values[0][2] * matrix[1][0] + values[1][2] * matrix[1][1] + values[2][2] * matrix[1][2],

			values[0][0] * matrix[2][0] + values[1][0] * matrix[2][1] + values[2][0] * matrix[2][2],
			values[0][1] * matrix[2][0] + values[1][1] * matrix[2][1] + values[2][1] * matrix[2][2],
			values[0][2] * matrix[2][0] + values[1][2] * matrix[2][1] + values[2][2] * matrix[2][2]
		);
	}

	Vector3 Matrix3x3::operator*(const Vector3 &vector) const
	{
		return Vector3(
			values[0][0] * vector.x + values[1][0] * vector.y + values[2][0] * vector.z,
			values[0][1] * vector.x + values[1][1] * vector.y + values[2][1] * vector.z,
			values[0][2] * vector.x + values[1][2] * vector.y + values[2][2] * vector.z
		);
	}

	Matrix1x3 Matrix3x3::operator*(const Matrix1x3 &matrix) const
	{
		Matrix1x3 res;
		res[0][1] = values[0][0] * matrix[0][0] + values[1][0] * matrix[0][1] + values[2][0] * matrix[0][2];
		res[0][1] = values[0][1] * matrix[0][0] + values[1][1] * matrix[0][1] + values[2][1] * matrix[0][2];
		res[0][2] = values[0][2] * matrix[0][0] + values[1][2] * matrix[0][1] + values[2][2] * matrix[0][2];
		return res;
	}

	Matrix2x3 Matrix3x3::operator*(const Matrix2x3 &matrix) const
	{
		Matrix2x3 res;
		res[0] = Vector3(
			values[0][0] * matrix[0][0] + values[1][0] * matrix[0][1] + values[2][0] * matrix[0][2],
			values[0][1] * matrix[0][0] + values[1][1] * matrix[0][1] + values[2][1] * matrix[0][2],
			values[0][2] * matrix[0][0] + values[1][2] * matrix[0][1] + values[2][2] * matrix[0][2]
		);
		res[1] = Vector3(
			values[0][0] * matrix[1][0] + values[1][0] * matrix[1][1] + values[2][0] * matrix[1][2],
			values[0][1] * matrix[1][0] + values[1][1] * matrix[1][1] + values[2][1] * matrix[1][2],
			values[0][2] * matrix[1][0] + values[1][2] * matrix[1][1] + values[2][2] * matrix[1][2]
		);
		return res;
	}

	Matrix4x3 Matrix3x3::operator*(const Matrix4x3 &matrix) const
	{
		Matrix4x3 res;
		res[0] = Vector3(
			values[0][0] * matrix[0][0] + values[1][0] * matrix[0][1] + values[2][0] * matrix[0][2],
			values[0][1] * matrix[0][0] + values[1][1] * matrix[0][1] + values[2][1] * matrix[0][2],
			values[0][2] * matrix[0][0] + values[1][2] * matrix[0][1] + values[2][2] * matrix[0][2]
		);
		res[1] = Vector3(
			values[0][0] * matrix[1][0] + values[1][0] * matrix[1][1] + values[2][0] * matrix[1][2],
			values[0][1] * matrix[1][0] + values[1][1] * matrix[1][1] + values[2][1] * matrix[1][2],
			values[0][2] * matrix[1][0] + values[1][2] * matrix[1][1] + values[2][2] * matrix[1][2]
		);
		res[2] = Vector3(
			values[0][0] * matrix[2][0] + values[1][0] * matrix[2][1] + values[2][0] * matrix[2][2],
			values[0][1] * matrix[2][0] + values[1][1] * matrix[2][1] + values[2][1] * matrix[2][2],
			values[0][2] * matrix[2][0] + values[1][2] * matrix[2][1] + values[2][2] * matrix[2][2]
		);
		res[3] = Vector3(
			values[0][0] * matrix[3][0] + values[1][0] * matrix[3][1] + values[2][0] * matrix[2][2],
			values[0][1] * matrix[3][0] + values[1][1] * matrix[3][1] + values[2][1] * matrix[2][2],
			values[0][2] * matrix[3][0] + values[1][2] * matrix[3][1] + values[2][2] * matrix[2][2]
		);
		return res;
	}

	template<unsigned int S>
	Matrix<S, 3, float> Matrix3x3::operator*(const Matrix<S, 3, float> &matrix) const
	{
		Matrix<S, 3, float> res;

		for (int i = 0; i < S; i++)
		{
			res[i] = Vector3(
				values[0][0] * matrix[i][0] + values[1][0] * matrix[i][1] + values[2][0] * matrix[i][2],
				values[0][1] * matrix[i][0] + values[1][1] * matrix[i][1] + values[2][1] * matrix[i][2],
				values[0][2] * matrix[i][0] + values[1][2] * matrix[i][1] + values[2][2] * matrix[i][2]
			);
		}

		return res;
	}

	Matrix3x3 Matrix3x3::operator*(const float scalar) const
	{
		return Matrix3x3(
			values[0] * scalar,
			values[1] * scalar,
			values[2] * scalar
		);
	}

	Matrix3x3 Matrix3x3::operator/(const Matrix3x3 &matrix) const
	{
		return *this * matrix.Inverted();
	}

	Matrix3x3 Matrix3x3::operator/(const float scalar) const
	{
		return Matrix3x3(
			values[0] / scalar,
			values[1] / scalar,
			values[2] / scalar
		);
	}

	// -- Boolean operators --

	bool Matrix3x3::operator==(const Matrix3x3 &matrix) const
	{
		return !(*this != matrix);
	}

	bool Matrix3x3::operator!=(const Matrix3x3 &matrix) const
	{
		return values[0][0] != matrix[0][0] || values[0][1] != matrix[0][1] || values[0][2] != matrix[0][2]
			|| values[1][0] != matrix[1][0] || values[1][1] != matrix[1][1] || values[1][2] != matrix[1][2]
			|| values[2][0] != matrix[2][0] || values[2][1] != matrix[2][1] || values[2][2] != matrix[2][2];
	}

	// -- Stream operators --

	std::ostream &operator<<(std::ostream &ostream, const Matrix3x3 &matrix)
	{
		return ostream	<< matrix[0][0] << ',' << matrix[1][0] << ',' << matrix[2][0] << '\n'
						<< matrix[0][1] << ',' << matrix[1][1] << ',' << matrix[2][1] << '\n'
						<< matrix[0][2] << ',' << matrix[1][2] << ',' << matrix[2][2] << '\n';
	}

	// -- Getters --

	Matrix3x3 Matrix3x3::GetAdjugate() const
	{
		return Matrix3x3(
			// Col 1
			values[1][1] * values[2][2] - values[2][1] * values[1][2],
			values[2][1] * values[0][2] - values[0][1] * values[2][2],
			values[0][1] * values[1][2] - values[1][1] * values[0][2],
			// Col 2
			values[2][0] * values[1][2] - values[1][0] * values[2][2],
			values[0][0] * values[2][2] - values[2][0] * values[0][2],
			values[1][0] * values[0][2] - values[0][0] * values[1][2],
			// Col 3
			values[1][0] * values[2][1] - values[2][0] * values[1][1],
			values[2][0] * values[0][1] - values[0][0] * values[2][1],
			values[0][0] * values[1][1] - values[1][0] * values[0][1]
		);
	}

	Matrix3x3 Matrix3x3::GetCofactor() const
	{
		return Matrix3x3(
			// Col 0
			values[1][1] * values[2][2] - values[2][1] * values[1][2],
			values[2][0] * values[1][2] - values[1][0] * values[2][2],
			values[1][0] * values[2][1] - values[2][0] * values[1][1],
			// Col 1
			values[2][1] * values[0][2] - values[0][1] * values[2][2],
			values[0][0] * values[2][2] - values[2][0] * values[0][2],
			values[2][0] * values[0][1] - values[0][0] * values[2][1],
			// Col 2
			values[0][1] * values[1][2] - values[1][1] * values[0][2],
			values[1][0] * values[0][2] - values[0][0] * values[1][2],
			values[0][0] * values[1][1] - values[1][0] * values[0][1]
		);
	}

	float Matrix3x3::GetDeterminant() const
	{
		return values[0][0] * (values[1][1] * values[2][2] - values[2][1] * values[1][2])
			-  values[1][0] * (values[0][1] * values[2][2] - values[2][1] * values[0][2])
			+  values[2][0] * (values[0][1] * values[1][2] - values[1][1] * values[0][2]);
	}

	Matrix3x3 Matrix3x3::Inverted() const
	{
		return GetAdjugate() / GetDeterminant();
	}

	Matrix3x3 Matrix3x3::Transposed() const
	{
		return Matrix3x3(
			values[0][0], values[1][0], values[2][0],
			values[0][1], values[1][1], values[2][1],
			values[0][2], values[1][2], values[2][2]
		);
	}

	// -- Transformations --

	Matrix3x3 &Matrix3x3::Invert()
	{
		*this = GetAdjugate() / GetDeterminant();
		return *this;
	}

	Matrix3x3 &Matrix3x3::Transpose()
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
}

