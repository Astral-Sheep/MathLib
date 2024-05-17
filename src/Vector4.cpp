#include "Vector4.hpp"
#include "Math.hpp"

namespace Math
{
#define IMPLEMENT_VEC4(Name, Type, Precision)\
	/* -- Constructors -- */\
	\
	Name::Vector()\
		: x(0), y(0), z(0), w(0)\
	{\
	\
	}\
	\
	Name::Vector(const Type x, const Type y, const Type z, const Type w)\
		: x(x), y(y), z(z), w(w)\
	{\
	\
	}\
	\
	Name::Vector(const Type fields[4])\
		: x(fields[0]), y(fields[1]), z(fields[2]), w(fields[3])\
	{\
	\
	}\
	\
	Name::Vector(const Name &vector)\
		: x(vector.x), y(vector.y), z(vector.z), w(vector.w)\
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
		w += vector.w;\
		return *this;\
	}\
	\
	Name &Name::operator-=(const Name &vector)\
	{\
		x -= vector.x;\
		y -= vector.y;\
		z -= vector.z;\
		w -= vector.w;\
		return *this;\
	}\
	\
	Name &Name::operator*=(const Name &vector)\
	{\
		x *= vector.x;\
		y *= vector.y;\
		z *= vector.z;\
		w *= vector.w;\
		return *this;\
	}\
	\
	Name &Name::operator*=(const Type scalar)\
	{\
		x *= scalar;\
		y *= scalar;\
		z *= scalar;\
		w *= scalar;\
		return *this;\
	}\
	\
	Name &Name::operator/=(const Name &vector)\
	{\
		x /= vector.x;\
		y /= vector.y;\
		z /= vector.z;\
		w /= vector.w;\
		return *this;\
	}\
	\
	Name &Name::operator/=(const Type scalar)\
	{\
		x /= scalar;\
		y /= scalar;\
		z /= scalar;\
		w /= scalar;\
		return *this;\
	}\
	\
	/* -- Unary operators -- */\
	\
	Name Name::operator-() const\
	{\
		return Name(-x, -y, -z, -w);\
	}\
	\
	/* -- Binary operators -- */\
	\
	Name Name::operator+(const Name &vector) const\
	{\
		return Name(x + vector.x, y + vector.y, z + vector.z, w + vector.w);\
	}\
	\
	Name Name::operator-(const Name &vector) const\
	{\
		return Name(x - vector.x, y - vector.y, z - vector.z, w - vector.w);\
	}\
	\
	Name Name::operator*(const Name &vector) const\
	{\
		return Name(x * vector.x, y * vector.y, z * vector.z, w * vector.w);\
	}\
	\
	Name Name::operator*(const Type scalar) const\
	{\
		return Name(x * scalar, y * scalar, z * scalar, w * scalar);\
	}\
	\
	Name Name::operator/(const Name &vector) const\
	{\
		return Name(x / vector.x, y / vector.y, z / vector.z, w / vector.w);\
	}\
	\
	Name Name::operator/(const Type scalar) const\
	{\
		return Name(x / scalar, y / scalar, z / scalar, w / scalar);\
	}\
	\
	/* -- Boolean operators -- */\
	\
	bool Name::operator==(const Name &vector) const\
	{\
		return x == vector.x && y == vector.y && z == vector.z && w == vector.w;\
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
		return ostream << "(" << vector.x << "," << vector.y << "," << vector.z << "," << vector.w << ")";\
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
		return Name(std::abs(x), std::abs(y), std::abs(z), std::abs(w));\
	}\
	\
	Precision Name::Distance(const Name &vector) const\
	{\
		const Precision to_x = vector.x - x;\
		const Precision to_y = vector.y - y;\
		const Precision to_z = vector.z - z;\
		const Precision to_w = vector.w - w;\
		return std::sqrt(to_x * to_x + to_y * to_y + to_z * to_z + to_w * to_w);\
	}\
	\
	Precision Name::DistanceSquared(const Name &vector) const\
	{\
		const Precision to_x = vector.x - x;\
		const Precision to_y = vector.y - y;\
		const Precision to_z = vector.z - z;\
		const Precision to_w = vector.w - w;\
		return to_x * to_x + to_y * to_y + to_z * to_z + to_w * to_w;\
	}\
	\
	Precision Name::Dot(const Name &vector) const\
	{\
		return x * vector.x + y * vector.y + z * vector.z + w * vector.w;\
	}\
	\
	bool Name::IsNormalized() const\
	{\
		return LengthSquared() == 1;\
	}\
	\
	Precision Name::Length() const\
	{\
		return std::sqrt(x * x + y * y + z * z + w * w);\
	}\
	\
	Precision Name::LengthSquared() const\
	{\
		return x * x + y * y + z * z + w * w;\
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
		return Name(x / length, y / length, z / length, w / length);\
	}\
	\
	Name Name::Sign() const\
	{\
		return Name(Math::Sign(x), Math::Sign(y), Math::Sign(z), Math::Sign(w));\
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
		w /= length;\
		return *this;\
	}\

	IMPLEMENT_VEC4(Vector4I, int, float)
	IMPLEMENT_VEC4(Vector4, float, float)
	IMPLEMENT_VEC4(Vector4D, double, double)
}

