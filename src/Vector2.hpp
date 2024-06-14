#pragma once

#include "Math.hpp"
#include "Vector.hpp"

namespace Math
{
	template<typename T, typename P>
	struct Vector<2, T, P>
	{
	private:
		typedef Vector<2, T, P> Vec;

	public:
		union { T x, r, s; };
		union { T y, g, t; };

		// -- Constructors --

		Vector<2, T, P>() : x(T(0)), y(T(0)) {}

		Vector<2, T, P>(const T val)
			: x(val), y(val) {}

		Vector<2, T, P>(const T x, const T y)
			: x(x), y(y) {}

		Vector<2, T, P>(const T vals[2])
			: x(vals[0]), y(vals[1]) {}

		Vector<2, T, P>(const Vector<2, T, P> &vec)
			: x(vec.x), y(vec.y) {}

		// -- Accesses --

		inline const T &operator[](const int index) const noexcept
		{
			return *(&x + index * sizeof(T));
		}

		inline T &operator[](const int index) noexcept
		{
			return *(&x + index * sizeof(T));
		}

		// -- Unary arithmetic operators --

		inline Vec &operator+=(const Vec &vec)
		{
			x += vec.x;
			y += vec.y;
			return *this;
		}

		inline Vec &operator-=(const Vec &vec)
		{
			x -= vec.x;
			y -= vec.y;
			return *this;
		}

		inline Vec &operator*=(const Vec &vec)
		{
			x *= vec.x;
			y *= vec.y;
			return *this;
		}

		inline Vec &operator*=(const T val)
		{
			x *= val;
			y *= val;
			return *this;
		}

		inline Vec &operator/=(const Vec &vec)
		{
			x /= vec.x;
			y /= vec.y;
			return *this;
		}

		inline Vec &operator/=(const T val)
		{
			x /= val;
			y /= val;
			return *this;
		}

		// -- Unary operators --

		inline Vec operator+() const
		{
			return *this;
		}

		inline Vec operator-() const
		{
			return Vector<2, T, P>(-x, -y);
		}

		// -- Binary operators --

		inline Vec operator+(const Vec &vec) const
		{
			return Vector<2, T, P>(x + vec.x, y + vec.y);
		}

		inline Vec operator-(const Vec &vec) const
		{
			return Vector<2, T, P>(x - vec.x, y - vec.y);
		}

		inline Vec operator*(const Vec &vec) const
		{
			return Vector<2, T, P>(x * vec.x, y * vec.y);
		}

		inline Vec operator*(const T val) const
		{
			return Vector<2, T, P>(x * val, y * val);
		}

		inline Vec operator/(const Vec &vec) const
		{
			return Vector<2, T, P>(x / vec.x, y / vec.y);
		}

		inline Vec operator/(const T val) const
		{
			return Vector<2, T, P>(x / val, y / val);
		}

		// -- Boolean operators --

		inline bool operator==(const Vec &vec) const
		{
			return x == vec.x && y == vec.y;
		}

		inline bool operator!=(const Vec &vec) const
		{
			return x != vec.x || y != vec.y;
		}

		// -- Convertion operators --

		template<typename U, typename Q>
		inline explicit operator Vector<2, U, Q>() const
		{
			return Vector<2, U, Q>((U)x, (U)y);
		}

		// -- Stream operators --

		friend std::ostream &operator<<(std::ostream &ostream, const Vec &vec)
		{
			return ostream << "(" << vec.x << ", " << vec.y << ")";
		}

		// -- Getters --

		Vec Abs() const
		{
			return Vector<2, T, P>(Math::Abs(x), Math::Abs(y));
		}

		T Cross(const Vec &vec) const
		{
			return x * vec.y - y * vec.x;
		}

		P Distance(const Vec &vec) const
		{
			const P to_x = vec.x - x;
			const P to_y = vec.y - y;
			return std::sqrt(to_x * to_x - to_y * to_y);
		}

		P DistanceSquared(const Vec &vec) const
		{
			const P to_x = vec.x - x;
			const P to_y = vec.y - y;
			return to_x * to_x + to_y * to_y;
		}

		T Dot(const Vec &vec) const
		{
			x * vec.x + y * vec.y;
		}

		bool IsNormalized() const
		{
			return x * x + y * y == T(1);
		}

		P Length() const
		{
			return std::sqrt(x * x + y * y);
		}

		T LengthSquared() const
		{
			return x * x + y * y;
		}

		Vec Normalized() const
		{
			P length = LengthSquared();

			if (length == 0 || length == 1)
			{
				return *this;
			}

			length = std::sqrt(length);
			return Vector<2, T, P>(x / length, y / length);
		}

		Vec Rotated(const P angle) const
		{
			const P cos = std::cos(angle);
			const P sin = std::sin(angle);
			return Vector<2, T, P>(cos * x - sin * y, sin * x + cos * y);
		}

		Vec Sign() const
		{
			return Vector<2, T, P>(Math::Sign(x), Math::Sign(y));
		}

		inline size_t SizeofField() const
		{
			return sizeof(T);
		}

		// -- Transformations --

		Vec &Normalize()
		{
			P length = LengthSquared();

			if (length == 0 || length == 1)
			{
				return *this;
			}

			length = std::sqrt(length);
			x /= length;
			y /= length;
			return *this;
		}

		Vec &Rotate(const P angle)
		{
			const P cos = std::cos(angle);
			const P sin = std::sin(angle);
			const T temp_x = x;
			const T temp_y = y;

			x = cos * temp_x - sin * temp_y;
			y = sin * temp_x + cos * temp_y;
			return *this;
		}

		// -- Static getters --

		static Vec Zero()
		{
			return Vector<2, T, P>(T(0), T(0));
		}

		static Vec One()
		{
			return Vector<2, T, P>(T(1), T(1));
		}

		static Vec NegOne()
		{
			return Vector<2, T, P>(T(-1), T(-1));
		}

		static Vec Left()
		{
			return Vector<2, T, P>(T(-1), T(0));
		}

		static Vec Right()
		{
			return Vector<2, T, P>(T(1), T(0));
		}

		static Vec Up()
		{
			return Vector<2, T, P>(T(0), T(1));
		}

		static Vec Down()
		{
			return Vector<2, T, P>(T(0), T(-1));
		}

		// -- Static methods --

		static Vec Lerp(const Vec &lhs, const Vec &rhs, const P t)
		{
			return Vector<2, T, P>(
				Math::Lerp(lhs.x, rhs.x, t),
				Math::Lerp(lhs.y, rhs.y, t)
			);
		}

		static Vec LerpClamped(const Vec &lhs, const Vec &rhs, const P t)
		{
			return Vector<2, T, P>(
				Math::LerpClamped(lhs.x, rhs.x, t),
				Math::LerpClamped(lhs.y, rhs.y, t)
			);
		}
	};

	typedef Vector<2, float, float> Vector2;
	typedef Vector<2, int, float> Vector2I;
	typedef Vector<2, double, double> Vector2D;
}

