#include "QuaternionD.hpp"
#include "Quaternion.hpp"
#include <cmath>

namespace Math
{
	// -- Constructors --

	QuaternionD::QuaternionD()
		: x(0.0), y(0.0), z(0.0), w(0.0)
	{

	}

	QuaternionD::QuaternionD(const double x, const double y, const double z, const double w)
		: x(x), y(y), z(z), w(w)
	{

	}

	QuaternionD::QuaternionD(const Vector3D &vector, const double scalar)
		: x(vector.x), y(vector.y), z(vector.z), w(scalar)
	{

	}

	QuaternionD::QuaternionD(const QuaternionD &quaternion)
		: x(quaternion.x), y(quaternion.y), z(quaternion.z), w(quaternion.w)
	{

	}

	// -- Unary arithmetic operators --

	QuaternionD &QuaternionD::operator+=(const QuaternionD &quaternion)
	{
		x += quaternion.x;
		y += quaternion.y;
		z += quaternion.z;
		w += quaternion.w;
		return *this;
	}

	QuaternionD &QuaternionD::operator-=(const QuaternionD &quaternion)
	{
		x -= quaternion.x;
		y -= quaternion.y;
		z -= quaternion.z;
		w -= quaternion.w;
		return *this;
	}

	QuaternionD &QuaternionD::operator*=(const QuaternionD &quaternion)
	{
		const double temp_x = x;
		const double temp_y = y;
		const double temp_z = z;
		const double temp_w = w;

		x = temp_w * quaternion.x + temp_x * quaternion.w + temp_y * quaternion.z - temp_z * quaternion.y;
		y = temp_w * quaternion.y - temp_x * quaternion.z + temp_y * quaternion.w + temp_z * quaternion.x;
		z = temp_w * quaternion.z + temp_x * quaternion.y - temp_y * quaternion.x + temp_z * quaternion.w;
		w = temp_w * quaternion.w - temp_x * quaternion.x - temp_y * quaternion.y - temp_z * quaternion.z;
		return *this;
	}

	QuaternionD &QuaternionD::operator*=(const double scalar)
	{
		x *= scalar;
		y *= scalar;
		z *= scalar;
		w *= scalar;
		return *this;
	}

	QuaternionD &QuaternionD::operator/=(const QuaternionD &quaternion)
	{
		*this *= quaternion.Inversed();
		return *this;
	}

	QuaternionD &QuaternionD::operator/=(const double scalar)
	{
		x /= scalar;
		y /= scalar;
		z /= scalar;
		w /= scalar;
		return *this;
	}

	// -- Unary operators --

	QuaternionD QuaternionD::operator-() const
	{
		return QuaternionD(-x, -y, -z, -w);
	}

	// -- Binary operators --

	QuaternionD QuaternionD::operator+(const QuaternionD &quaternion) const
	{
		return QuaternionD(x + quaternion.x, y + quaternion.y, z + quaternion.z, w + quaternion.w);
	}

	QuaternionD QuaternionD::operator-(const QuaternionD &quaternion) const
	{
		return QuaternionD(x - quaternion.x, y - quaternion.y, z - quaternion.z, w - quaternion.w);
	}

	QuaternionD QuaternionD::operator*(const QuaternionD &quaternion) const
	{
		return QuaternionD(
			w * quaternion.x + x * quaternion.w + y * quaternion.z - z * quaternion.y,
			w * quaternion.y - x * quaternion.z + y * quaternion.w + z * quaternion.x,
			w * quaternion.z + x * quaternion.y - y * quaternion.x + z * quaternion.w,
			w * quaternion.w - x * quaternion.x - y * quaternion.y - z * quaternion.z
		);
	}

	Vector3D QuaternionD::operator*(const Vector3D &vector) const
	{
		return Vector3D(
			(1.0 - 2.0 * y * y - 2.0 * z * z) * vector.x + 2.0 * (x * y + w * z) * vector.y + 2.0 * (x * z - w * y) * vector.z,
			2.0 * (x * y - w * z) * vector.x + (1.0 - 2.0 * x * x - 2.0 * z * z) * vector.y + 2.0 * (x * z + w * y) * vector.z,
			2.0 * (x * z + w * y) * vector.x + 2.0 * (y * z - w * x) * vector.y + (1.0 - 2.0 * x * x - 2.0 * y * y) * vector.z
		);
	}

	QuaternionD QuaternionD::operator*(const double scalar) const
	{
		return QuaternionD(x * scalar, y * scalar, z * scalar, w * scalar);
	}

	QuaternionD QuaternionD::operator/(const QuaternionD &quaternion) const
	{
		return *this * quaternion.Inversed();
	}

	QuaternionD QuaternionD::operator/(const double scalar) const
	{
		return QuaternionD(x / scalar, y / scalar, z / scalar, w / scalar);
	}

	// -- Boolean operators --

	bool QuaternionD::operator==(const QuaternionD &quaternion) const
	{
		return x == quaternion.x && y == quaternion.y && z == quaternion.z && w == quaternion.w;
	}

	bool QuaternionD::operator!=(const QuaternionD &quaternion) const
	{
		return !(*this == quaternion);
	}

	// -- Convertion operators

	QuaternionD::operator Quaternion() const
	{
		return Quaternion(x, y, z, w);
	}

	// -- Stream operators --

	std::ostream &operator<<(std::ostream &ostream, const QuaternionD &quaternion)
	{
		return ostream << "(" << quaternion.x << "," << quaternion.y << "," << quaternion.z << "," << quaternion.w << ")";
	}

	// -- Static methods --

	QuaternionD QuaternionD::AngleAxis(Vector3D axis, const double angle)
	{
		if (!axis.IsNormalized())
		{
			axis.Normalize();
		}

		const double sin = std::sin(angle * 0.5);
		return QuaternionD(axis.x * sin, axis.y * sin, axis.z * sin, std::cos(angle * 0.5));
	}

	// -- Getters --

	QuaternionD QuaternionD::Inversed() const
	{
		const double sqrlength = MagnitudeSquared();
		return QuaternionD(
			-x / sqrlength,
			-y / sqrlength,
			-z / sqrlength,
			w / sqrlength
		);
	}

	double QuaternionD::Magnitude() const
	{
		return std::sqrt(x * x + y * y + z * z + w * w);
	}

	double QuaternionD::MagnitudeSquared() const
	{
		return x * x + y * y + z * z + w * w;
	}

	// -- Transformations --

	QuaternionD &QuaternionD::Conjugate()
	{
		x = -x;
		y = -y;
		z = -z;
		return *this;
	}

	QuaternionD &QuaternionD::Inverse()
	{
		const double sqrlength = MagnitudeSquared();
		x = -x / sqrlength;
		y = -y / sqrlength;
		z = -z / sqrlength;
		w /= sqrlength;
		return *this;
	}
}

