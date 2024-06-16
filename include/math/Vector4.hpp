#pragma once

#include "Math.hpp"
#include "Vector.hpp"

namespace Math
{
	template<typename T, typename P>
	struct Vector<4, T, P>
	{
	private:
		typedef Vector<4, T, P> Vec;

	public:
		union { T x, r, s; };
		union { T y, g, t; };
		union { T z, b, p; };
		union { T w, a, q; };

		// -- Constructors --

		Vector<4, T, P>()
			: x(T(0)), y(T(0)), z(T(0)), w(T(0)) {}

		explicit Vector<4, T, P>(const T val)
			: x(val), y(val), z(val), w(val) {}

		Vector <4, T, P>(const T x, const T y, const T z, const T w)
			: x(x), y(y), z(z), w(w) {}

		Vector<4, T, P>(const T vals[4])
			: x(vals[0]), y(vals[1]), z(vals[2]), w(vals[3]) {}

		Vector<4, T, P>(const Vec &vec)
			: x(vec.x), y(vec.y), z(vec.z), w(vec.w) {}

		// -- Accesses --

		inline const T &operator[](const int index) const noexcept
		{
			return *(&x + index);
		}

		inline T &operator[](const int index) noexcept
		{
			return *(&x + index);
		}

		// -- Unary arithmetic operators --

		inline Vec &operator+=(const Vec &vec)
		{
			x += vec.x;
			y += vec.y;
			z += vec.z;
			w += vec.w;
			return *this;
		}

		inline Vec &operator-=(const Vec &vec)
		{
			x -= vec.x;
			y -= vec.y;
			z -= vec.z;
			w -= vec.w;
			return *this;
		}

		inline Vec &operator*=(const Vec &vec)
		{
			x *= vec.x;
			y *= vec.y;
			z *= vec.z;
			w *= vec.w;
			return *this;
		}

		inline Vec &operator*=(const T val)
		{
			x *= val;
			y *= val;
			z *= val;
			w *= val;
			return *this;
		}

		inline Vec &operator/=(const Vec &vec)
		{
			x /= vec.x;
			y /= vec.y;
			z /= vec.z;
			w /= vec.w;
			return *this;
		}

		inline Vec &operator/=(const T val)
		{
			x /= val;
			y /= val;
			z /= val;
			w /= val;
			return *this;
		}

		// -- Unary operators --

		inline Vec operator+() const
		{
			return *this;
		}

		inline Vec operator-() const
		{
			return Vector<4, T, P>(-x, -y, -z, -w);
		}

		// -- Binary operators --

		inline Vec operator+(const Vec &vec) const
		{
			return Vector<4, T, P>(x + vec.x, y + vec.y, z + vec.z, w + vec.w);
		}

		inline Vec operator-(const Vec &vec) const
		{
			return Vector<4, T, P>(x - vec.x, y - vec.y, z - vec.z, w - vec.w);
		}

		inline Vec operator*(const Vec &vec) const
		{
			return Vector<4, T, P>(x * vec.x, y * vec.y, z * vec.z, w * vec.w);
		}

		inline Vec operator*(const T val) const
		{
			return Vector<4, T, P>(x * val, y * val, z * val, w * val);
		}

		inline Vec operator/(const Vec &vec) const
		{
			return Vector<4, T, P>(x / vec.x, y / vec.y, z / vec.z, w / vec.w);
		}

		inline Vec operator/(const T val) const
		{
			return Vector<4, T, P>(x / val, y / val, z / val, w / val);
		}

		// -- Boolean operators --

		inline bool operator==(const Vec &vec) const
		{
			return x == vec.x && y == vec.y && z == vec.z && w == vec.w;
		}

		inline bool operator!=(const Vec &vec) const
		{
			return x != vec.x || y != vec.y || z != vec.z || w != vec.w;
		}

		// -- Convertion operators --

		template<typename U, typename Q>
		inline explicit operator Vector<4, U, Q>() const
		{
			return Vector<4, U, Q>((U)x, (U)y, (U)z, (U)w);
		}

		// -- Stream operators --

		friend std::ostream &operator<<(std::ostream &ostream, const Vec &vec)
		{
			return ostream << "(" << vec.x << ", " << vec.y << ", " << vec.z << ", " << vec.w << ")";
		}

		// -- Getters --

		Vec Abs() const
		{
			return Vector<4, T, P>(Math::Abs(x), Math::Abs(y), Math::Abs(z), Math::Abs(w));
		}

		P Distance(const Vec &vec) const
		{
			const P to_x = vec.x - x;
			const P to_y = vec.y - y;
			const P to_z = vec.z - z;
			const P to_w = vec.w - w;
			return std::sqrt(to_x * to_x + to_y * to_y + to_z * to_z + to_w * to_w);
		}

		P DistanceSquared(const Vec &vec) const
		{
			const P to_x = vec.x - x;
			const P to_y = vec.y - y;
			const P to_z = vec.z - z;
			const P to_w = vec.w - w;
			return to_x * to_x + to_y * to_y + to_z * to_z + to_w * to_w;
		}

		T Dot(const Vec &vec) const
		{
			return x * vec.x + y * vec.y + z * vec.z + w * vec.w;
		}

		bool IsNormalized() const
		{
			return x * x + y * y + z * z + w * w == T(1);
		}

		P Length() const
		{
			return std::sqrt(x * x + y * y + z * z + w * w);
		}

		T LengthSquared() const
		{
			return x * x + y * y + z * z + w * w;
		}

		Vec Normalized() const
		{
			P length = LengthSquared();

			if (length == 0 || length == 1)
			{
				return *this;
			}

			length = std::sqrt(length);
			return Vector<4, T, P>(x / length, y / length, z / length, w / length);
		}

		Vec Sign() const
		{
			return Vector<4, T, P>(Math::Sign(x), Math::Sign(y), Math::Sign(z), Math::Sign(w));
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
			z /= length;
			w /= length;
			return *this;
		}

		// -- Static getters --

		static inline Vec Zero()
		{
			return Vector<4, T, P>(T(0), T(0), T(0), T(0));
		}

		static inline Vec One()
		{
			return Vector<4, T, P>(T(1), T(1), T(1), T(1));
		}

		static inline Vec NegOne()
		{
			return Vector<4, T, P>(T(-1), T(-1), T(-1), T(-1));
		}

		static inline Vec PosX()
		{
			return Vector<4, T, P>(T(1), T(0), T(0), T(0));
		}

		static inline Vec NegX()
		{
			return Vector<4, T, P>(T(-1), T(0), T(0), T(0));
		}

		static inline Vec PosY()
		{
			return Vector<4, T, P>(T(0), T(1), T(0), T(0));
		}

		static inline Vec NegY()
		{
			return Vector<4, T, P>(T(0), T(-1), T(0), T(0));
		}

		static inline Vec PosZ()
		{
			return Vector<4, T, P>(T(0), T(0), T(1), T(0));
		}

		static inline Vec NegZ()
		{
			return Vector<4, T, P>(T(0), T(0), T(-1), T(0));
		}

		static inline Vec PosW()
		{
			return Vector<4, T, P>(T(0), T(0), T(0), T(1));
		}

		static inline Vec NegW()
		{
			return Vector<4, T, P>(T(0), T(0), T(0), T(-1));
		}

		// -- Static methods --

		static Vec Lerp(const Vec &lhs, const Vec &rhs, const P t)
		{
			return Vector<4, T, P>(
				Math::Lerp(lhs.x, rhs.x, t),
				Math::Lerp(lhs.y, rhs.y, t),
				Math::Lerp(lhs.z, rhs.z, t),
				Math::Lerp(lhs.w, rhs.w, t)
			);
		}

		static Vec LerpClamped(const Vec &lhs, const Vec &rhs, const P t)
		{
			return Vector<4, T, P>(
				Math::LerpClamped(lhs.x, rhs.x, t),
				Math::LerpClamped(lhs.y, rhs.y, t),
				Math::LerpClamped(lhs.z, rhs.z, t),
				Math::LerpClamped(lhs.w, rhs.w, t)
			);
		}
	};

	typedef Vector<4, float, float> Vector4;
	typedef Vector<4, int, float> Vector4I;
	typedef Vector<4, double, double> Vector4D;
}

