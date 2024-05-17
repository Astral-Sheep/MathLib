#include "Vector2.hpp"
#include "Math.hpp"

namespace Math
{
#define IMPLEMENT_VEC2(Name, Type, Precision)\
	/* -- Constructors -- */\
	\
	Name::Vector()\
		: x(0), y(0)\
	{\
	\
	}\
	\
	Name::Vector(const Type x, const Type y)\
		: x(x), y(y)\
	{\
	\
	}\
	\
	Name::Vector(const Type fields[2])\
		: x(fields[0]), y(fields[1])\
	{\
	\
	}\
	\
	Name::Vector(const Name &vector)\
		: x(vector.x), y(vector.y)\
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
		return *this;\
	}\
	\
	Name &Name::operator-=(const Name &vector)\
	{\
		x -= vector.x;\
		y -= vector.y;\
		return *this;\
	}\
	\
	Name &Name::operator*=(const Name &vector)\
	{\
		x *= vector.x;\
		y *= vector.y;\
		return *this;\
	}\
	\
	Name &Name::operator*=(const Type scalar)\
	{\
		x *= scalar;\
		y *= scalar;\
		return *this;\
	}\
	\
	Name &Name::operator/=(const Name &vector)\
	{\
		x /= vector.x;\
		y /= vector.y;\
		return *this;\
	}\
	\
	Name &Name::operator/=(const Type scalar)\
	{\
		x /= scalar;\
		y /= scalar;\
		return *this;\
	}\
	\
	/* -- Unary operators -- */\
	\
	Name Name::operator-() const\
	{\
		return Name(-x, -y);\
	}\
	\
	/* -- Binary operators -- */\
	\
	Name Name::operator+(const Name &vector) const\
	{\
		return Name(x + vector.x, y + vector.y);\
	}\
	\
	Name Name::operator-(const Name &vector) const\
	{\
		return Name(x - vector.x, y - vector.y);\
	}\
	\
	Name Name::operator*(const Name &vector) const\
	{\
		return Name(x * vector.x, y * vector.y);\
	}\
	\
	Name Name::operator*(const Type scalar) const\
	{\
		return Name(x * scalar, y * scalar);\
	}\
	\
	Name Name::operator/(const Name &vector) const\
	{\
		return Name(x / vector.x, y / vector.y);\
	}\
	\
	Name Name::operator/(const Type scalar) const\
	{\
		return Name(x / scalar, y / scalar);\
	}\
	\
	/* -- Boolean operators -- */\
	\
	bool Name::operator==(const Name &vector) const\
	{\
		return x == vector.x && y == vector.y;\
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
		return ostream << "(" << vector.x << "," << vector.y << ")";\
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
		return Name(std::abs(x), std::abs(y));\
	}\
	\
	Precision Name::Cross(const Name &vector) const\
	{\
		return x * vector.y - y * vector.x;\
	}\
	\
	Precision Name::Distance(const Name &vector) const\
	{\
		const Precision to_x = vector.x - x;\
		const Precision to_y = vector.y - y;\
		return std::sqrt(to_x * to_x + to_y * to_y);\
	}\
	\
	Precision Name::DistanceSquared(const Name &vector) const\
	{\
		const Precision to_x = vector.x - x;\
		const Precision to_y = vector.y - y;\
		return to_x * to_x + to_y * to_y;\
	}\
	\
	Precision Name::Dot(const Name &vector) const\
	{\
		return x * vector.x + y * vector.y;\
	}\
	\
	bool Name::IsNormalized() const\
	{\
		return LengthSquared() == 1;\
	}\
	\
	Precision Name::Length() const\
	{\
		return std::sqrt(x * x + y * y);\
	}\
	\
	Precision Name::LengthSquared() const\
	{\
		return x * x + y * y;\
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
		return Name(x / length, y / length);\
	}\
	\
	Name Name::Rotated(const Precision angle) const\
	{\
		Precision cos = std::cos(angle);\
		Precision sin = std::sin(angle);\
		return Name(cos * x - sin * y, sin * x + cos * y);\
	}\
	\
	Name Name::Sign() const\
	{\
		return Name(Math::Sign(x), Math::Sign(y));\
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
		return *this;\
	}\
	\
	Name &Name::Rotate(const Precision angle)\
	{\
		Precision cos = std::cos(angle);\
		Precision sin = std::sin(angle);\
		Precision temp_x = x;\
		Precision temp_y = y;\
		x = cos * temp_x - sin * temp_y;\
		y = sin * temp_x + cos * temp_y;\
		return *this;\
	}

	IMPLEMENT_VEC2(Vector2I, int, float)
	IMPLEMENT_VEC2(Vector2, float, float)
	IMPLEMENT_VEC2(Vector2D, double, double)
}
