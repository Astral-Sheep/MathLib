#include "Quaternion.hpp"
#include "QuaternionD.hpp"
#include <cmath>

namespace Math
{
	// -- Constructors --

	Quaternion::Quaternion()
		: x(0.0f), y(0.0f), z(0.0f), w(0.0f)
	{

	}

	Quaternion::Quaternion(const float x, const float y, const float z, const float w)
		: x(x), y(y), z(z), w(w)
	{

	}

	Quaternion::Quaternion(const Vector3 &vector, const float scalar)
		: x(vector.x), y(vector.y), z(vector.z), w(scalar)
	{

	}

	Quaternion::Quaternion(const Quaternion &quaternion)
		: x(quaternion.x), y(quaternion.y), z(quaternion.z), w(quaternion.w)
	{

	}

	// -- Unary arithmetic operators --

	Quaternion &Quaternion::operator+=(const Quaternion &quaternion)
	{
		x += quaternion.x;
		y += quaternion.y;
		z += quaternion.z;
		w += quaternion.w;
		return *this;
	}

	Quaternion &Quaternion::operator-=(const Quaternion &quaternion)
	{
		x -= quaternion.x;
		y -= quaternion.y;
		z -= quaternion.z;
		w -= quaternion.w;
		return *this;
	}

	Quaternion &Quaternion::operator*=(const Quaternion &quaternion)
	{
		const float temp_x = x;
		const float temp_y = y;
		const float temp_z = z;
		const float temp_w = w;

		x = temp_w * quaternion.x + temp_x * quaternion.w + temp_y * quaternion.z - temp_z * quaternion.y;
		y = temp_w * quaternion.y - temp_x * quaternion.z + temp_y * quaternion.w + temp_z * quaternion.x;
		z = temp_w * quaternion.z + temp_x * quaternion.y - temp_y * quaternion.x + temp_z * quaternion.w;
		w = temp_w * quaternion.w - temp_x * quaternion.x - temp_y * quaternion.y - temp_z * quaternion.z;
		return *this;
	}

	Quaternion &Quaternion::operator*=(const float scalar)
	{
		x *= scalar;
		y *= scalar;
		z *= scalar;
		w *= scalar;
		return *this;
	}

	Quaternion &Quaternion::operator/=(const Quaternion &quaternion)
	{
		*this *= quaternion.Inversed();
		return *this;
	}

	Quaternion &Quaternion::operator/=(const float scalar)
	{
		x /= scalar;
		y /= scalar;
		z /= scalar;
		w /= scalar;
		return *this;
	}

	// -- Unary operators --

	Quaternion Quaternion::operator-() const
	{
		return Quaternion(-x, -y, -z, -w);
	}

	// -- Binary operators --

	Quaternion Quaternion::operator+(const Quaternion &quaternion) const
	{
		return Quaternion(x + quaternion.x, y + quaternion.y, z + quaternion.z, w + quaternion.w);
	}

	Quaternion Quaternion::operator-(const Quaternion &quaternion) const
	{
		return Quaternion(x - quaternion.x, y - quaternion.y, z - quaternion.z, w - quaternion.w);
	}

	Quaternion Quaternion::operator*(const Quaternion &quaternion) const
	{
		return Quaternion(
			w * quaternion.x + x * quaternion.w + y * quaternion.z - z * quaternion.y,
			w * quaternion.y - x * quaternion.z + y * quaternion.w + z * quaternion.x,
			w * quaternion.z + x * quaternion.y - y * quaternion.x + z * quaternion.w,
			w * quaternion.w - x * quaternion.x - y * quaternion.y - z * quaternion.z
		);
	}

	Vector3 Quaternion::operator*(const Vector3 &vector) const
	{
		return Vector3(
			(1.0f - 2.0f * y * y - 2.0f * z * z) * vector.x + 2.0f * (x * y + w * z) * vector.y + 2.0f * (x * z - w * y) * vector.z,
			2.0f * (x * y - w * z) * vector.x + (1.0f - 2.0f * x * x - 2.0f * z * z) * vector.y + 2.0f * (x * z + w * y) * vector.z,
			2.0f * (x * z + w * y) * vector.x + 2.0f * (y * z - w * x) * vector.y + (1.0f - 2.0f * x * x - 2.0f * y * y) * vector.z
		);
	}

	Quaternion Quaternion::operator*(const float scalar) const
	{
		return Quaternion(x * scalar, y * scalar, z * scalar, w * scalar);
	}

	Quaternion Quaternion::operator/(const Quaternion &quaternion) const
	{
		return *this * quaternion.Inversed();
	}

	Quaternion Quaternion::operator/(const float scalar) const
	{
		return Quaternion(x / scalar, y / scalar, z / scalar, w / scalar);
	}

	// -- Boolean operators --

	bool Quaternion::operator==(const Quaternion &quaternion) const
	{
		return x == quaternion.x && y == quaternion.y && z == quaternion.z && w == quaternion.w;
	}

	bool Quaternion::operator!=(const Quaternion &quaternion) const
	{
		return !(*this == quaternion);
	}

	// -- Convertion operators --

	Quaternion::operator QuaternionD() const
	{
		return QuaternionD(x, y, z, w);
	}

	// -- Stream operators --

	std::ostream &operator<<(std::ostream &ostream, const Quaternion &quaternion)
	{
		return ostream << "(" << quaternion.x << "," << quaternion.y << "," << quaternion.z << "," << quaternion.w << ")";
	}

	// -- Static methods --

	Quaternion Quaternion::AngleAxis(Vector3 axis, const float angle)
	{
		if (!axis.IsNormalized())
		{
			axis.Normalize();
		}

		const float sin = std::sin(angle * 0.5f);
		return Quaternion(axis.x * sin, axis.y * sin, axis.z * sin, std::cos(angle * 0.5f));
	}

	// -- Getters --

	Quaternion Quaternion::Inversed() const
	{
		const float sqrlength = MagnitudeSquared();
		return Quaternion(-x / sqrlength, -y / sqrlength, -z / sqrlength, w / sqrlength);
	}

	float Quaternion::Magnitude() const
	{
		return std::sqrt(x * x + y * y + z * z + w * w);
	}

	float Quaternion::MagnitudeSquared() const
	{
		return x * x + y * y + z * z + w * w;
	}

	// -- Transformations --

	Quaternion &Quaternion::Conjugate()
	{
		x = -x;
		y = -y;
		z = -z;
		return *this;
	}

	Quaternion &Quaternion::Inverse()
	{
		const float sqrlength = MagnitudeSquared();
		x = -x / sqrlength;
		y = -y / sqrlength;
		z = -z / sqrlength;
		w /= sqrlength;
		return *this;
	}
}

