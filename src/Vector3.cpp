#include "Vector3.hpp"
#include "Math.hpp"

namespace Math
{
#define IMPLEMENT_VEC3(Name, Type, Precision)\
	/* -- Constructors -- */\
	\
	Name::Vector()\
		: x(0), y(0), z(0)\
	{\
	\
	}\
	\
	Name::Vector(const Type x, const Type y, const Type z)\
		: x(x), y(y), z(z)\
	{\
	\
	}\
	\
	Name::Vector(const Type fields[3])\
		: x(fields[0]), y(fields[1]), z(fields[2])\
	{\
	\
	}\
	\
	Name::Vector(const Name &vector)\
		: x(vector.x), y(vector.y), z(vector.z)\
	{\
	\
	}\
	\
	/* -- Accesses -- */\
	\
	const Type &Name::operator[](const int index) const noexcept\
	{\
		return *(&x + index);\
	}\
	\
	Type &Name::operator[](const int index) noexcept\
	{\
		return *(&x + index);\
	}\
	\
	/* -- Unary arithmetic operators -- */\
	\
	Name &Name::operator+=(const Name &vector)\
	{\
		x += vector.x;\
		y += vector.y;\
		z += vector.z;\
		return *this;\
	}\
	\
	Name &Name::operator-=(const Name &vector)\
	{\
		x -= vector.x;\
		y -= vector.y;\
		z -= vector.z;\
		return *this;\
	}\
	\
	Name &Name::operator*=(const Name &vector)\
	{\
		x *= vector.x;\
		y *= vector.y;\
		z *= vector.z;\
		return *this;\
	}\
	\
	Name &Name::operator*=(const Type scalar)\
	{\
		x *= scalar;\
		y *= scalar;\
		z *= scalar;\
		return *this;\
	}\
	\
	Name &Name::operator/=(const Name &vector)\
	{\
		x /= vector.x;\
		y /= vector.y;\
		z /= vector.z;\
		return *this;\
	}\
	\
	Name &Name::operator/=(const Type scalar)\
	{\
		x /= scalar;\
		y /= scalar;\
		z /= scalar;\
		return *this;\
	}\
	\
	/* -- Unary operators -- */\
	\
	Name Name::operator-() const\
	{\
		return Name(-x, -y, -z);\
	}\
	\
	/* -- Binary operators -- */\
	\
	Name Name::operator+(const Name &vector) const\
	{\
		return Name(x + vector.x, y + vector.y, z + vector.z);\
	}\
	\
	Name Name::operator-(const Name &vector) const\
	{\
		return Name(x - vector.x, y - vector.y, z - vector.z);\
	}\
	\
	Name Name::operator*(const Name &vector) const\
	{\
		return Name(x * vector.x, y * vector.y, z * vector.z);\
	}\
	\
	Name Name::operator*(const Type scalar) const\
	{\
		return Name(x * scalar, y * scalar, z * scalar);\
	}\
	\
	Name Name::operator/(const Name &vector) const\
	{\
		return Name(x / vector.x, y / vector.y, z / vector.z);\
	}\
	\
	Name Name::operator/(const Type scalar) const\
	{\
		return Name(x / scalar, y / scalar, z / scalar);\
	}\
	\
	/* -- Boolean operators -- */\
	\
	bool Name::operator==(const Name &vector) const\
	{\
		return x == vector.x && y == vector.y && z == vector.z;\
	}\
	\
	bool Name::operator!=(const Name &vector) const\
	{\
		return !(*this == vector);\
	}\
	\
	/* -- Stream operators -- */\
	\
	std::ostream &operator<<(std::ostream &ostream, const Name &vector)\
	{\
		return ostream << "(" << vector.x << "," << vector.y << "," << vector.z << ")";\
	}\
	\
	/* -- Static methods -- */\
	\
	Name Name::Lerp(const Name &from, const Name &to, const Precision t)\
	{\
		return from + (to - from) * t;\
	}\
	\
	Name Name::LerpClamped(const Name &from, const Name &to, const Precision t)\
	{\
		return from + (to - from) * Clamp01(t);\
	}\
	\
	/* -- Getters -- */\
	\
	Name Name::Abs() const\
	{\
		return Name(std::abs(x), std::abs(y), std::abs(z));\
	}\
	\
	Name Name::Cross(const Name &vector) const\
	{\
		return Name(\
			y * vector.z - z * vector.y,\
			z * vector.x - x * vector.z,\
			x * vector.y - y * vector.x\
		);\
	}\
	\
	Precision Name::Distance(const Name &vector) const\
	{\
		const Precision to_x = vector.x - x;\
		const Precision to_y = vector.y - y;\
		const Precision to_z = vector.z - z;\
		return std::sqrt(to_x * to_x + to_y * to_y + to_z * to_z);\
	}\
	\
	Precision Name::DistanceSquared(const Name &vector) const\
	{\
		const Precision to_x = vector.x - x;\
		const Precision to_y = vector.y - y;\
		const Precision to_z = vector.z - z;\
		return to_x * to_x + to_y * to_y + to_z * to_z;\
	}\
	\
	Precision Name::Dot(const Name &vector) const\
	{\
		return x * vector.x + y * vector.y + z * vector.z;\
	}\
	\
	bool Name::IsNormalized() const\
	{\
		return LengthSquared() == 1;\
	}\
	\
	Precision Name::Length() const\
	{\
		return std::sqrt(x * x + y * y + z * z);\
	}\
	\
	Precision Name::LengthSquared() const\
	{\
		return x * x + y * y + z * z;\
	}\
	\
	Name Name::Normalized() const\
	{\
		Precision length = LengthSquared();\
		\
		if (length == 0 || length == 1)\
		{\
			return *this;\
		}\
		\
		length = std::sqrt(length);\
		return Name(x / length, y / length, z / length);\
	}\
	\
	Name Name::Rotated(const Vector<3, Precision, Precision> &rotation, const EulerAngleOrder order) const\
	{\
		Name copy = Name(*this);\
		copy.Rotate(rotation, order);\
		return copy;\
	}\
	\
	Name Name::Sign() const\
	{\
		return Name(Math::Sign(x), Math::Sign(y), Math::Sign(z));\
	}\
	\
	size_t Name::SizeOfField() const\
	{\
		return sizeof(Type);\
	}\
	\
	/* -- Transformations -- */\
	\
	Name &Name::Normalize()\
	{\
		Precision length = LengthSquared();\
		\
		if (length == 0 || length == 1)\
		{\
			return *this;\
		}\
		\
		length = std::sqrt(length);\
		x /= length;\
		y /= length;\
		z /= length;\
		return *this;\
	}\
	\
	Name &Name::RotateAroundX(const Precision cosAngle, const Precision sinAngle)\
	{\
		Type temp_y = y;\
		Type temp_z = z;\
		y = cosAngle * temp_y - sinAngle * temp_z;\
		z = sinAngle * temp_y + cosAngle * temp_z;\
		return *this;\
	}\
	\
	Name &Name::RotateAroundX(const Precision angle)\
	{\
		return RotateAroundX(std::cos(angle), std::sin(angle));\
	}\
	\
	Name &Name::RotateAroundY(const Precision cosAngle, const Precision sinAngle)\
	{\
		Type temp_x = x;\
		Type temp_z = z;\
		x = cosAngle * temp_x + sinAngle * temp_z;\
		z = -sinAngle * temp_x + cosAngle * temp_z;\
		return *this;\
	}\
	\
	Name &Name::RotateAroundY(const Precision angle)\
	{\
		return RotateAroundY(std::cos(angle), std::sin(angle));\
	}\
	\
	Name &Name::RotateAroundZ(const Precision cosAngle, const Precision sinAngle)\
	{\
		Type temp_x = x;\
		Type temp_y = y;\
		x = cosAngle * temp_x - sinAngle * temp_y;\
		y = sinAngle * temp_x + cosAngle * temp_y;\
		return *this;\
	}\
	\
	Name &Name::RotateAroundZ(const Precision angle)\
	{\
		return RotateAroundZ(std::cos(angle), std::sin(angle));\
	}\
	\
	Name &Name::Rotate(const Vector<3, Precision, Precision> &rotation, const EulerAngleOrder order)\
	{\
		/*To refactor*/\
		switch (order)\
		{\
			case EulerAngleOrder::XYZ:\
				RotateAroundX(rotation.x);\
				RotateAroundY(rotation.y);\
				RotateAroundZ(rotation.z);\
				break;\
			case EulerAngleOrder::XZY:\
				RotateAroundX(rotation.x);\
				RotateAroundZ(rotation.z);\
				RotateAroundY(rotation.y);\
				break;\
			case EulerAngleOrder::YXZ:\
				RotateAroundY(rotation.y);\
				RotateAroundX(rotation.x);\
				RotateAroundZ(rotation.z);\
				break;\
			case EulerAngleOrder::YZX:\
				RotateAroundY(rotation.y);\
				RotateAroundZ(rotation.z);\
				RotateAroundX(rotation.x);\
				break;\
			case EulerAngleOrder::ZXY:\
				RotateAroundZ(rotation.z);\
				RotateAroundX(rotation.x);\
				RotateAroundY(rotation.y);\
				break;\
			case EulerAngleOrder::ZYX:\
				RotateAroundZ(rotation.z);\
				RotateAroundY(rotation.y);\
				RotateAroundX(rotation.x);\
				break;\
			default:\
				break;\
		}\
		\
		return *this;\
	}

	IMPLEMENT_VEC3(Vector3I, int, float)
	IMPLEMENT_VEC3(Vector3, float, float)
	IMPLEMENT_VEC3(Vector3D, double, double)
}

