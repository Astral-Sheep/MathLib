#pragma once

#include "Math.hpp"
#include "Vector.hpp"

namespace Math
{
	enum class EulerAngleOrder
	{
		XYZ,
		XZY,
		YXZ,
		YZX,
		ZXY,
		ZYX,
	};

	template<typename T, typename P>
	struct Vector<3, T, P>
	{
	private:
		typedef Vector<3, T, P> Vec;

		Vec &RotateAroundX(const P cos, const P sin)
		{
			T temp_y = y;
			T temp_z = z;
			y = cos * temp_y - sin * temp_z;
			z = sin * temp_y + cos * temp_z;
			return *this;
		}

		Vec &RotateAroundY(const P cos, const P sin)
		{
			T temp_x = x;
			T temp_z = z;
			x = cos * temp_x + sin * temp_z;
			z = -sin * temp_x + cos * temp_z;
			return *this;
		}

		Vec &RotateAroundZ(const P cos, const P sin)
		{
			T temp_x = x;
			T temp_y = y;
			x = cos * temp_x - sin * temp_y;
			y = sin * temp_x + cos * temp_y;
			return *this;
		}

	public:
		union { T x, r, s; };
		union { T y, g, t; };
		union { T z, b, p; };

		// -- Constructors --

		Vector<3, T, P>()
			: x(T(0)), y(T(0)), z(T(0)) {}

		explicit Vector<3, T, P>(const T val)
			: x(val), y(val), z(val) {}

		Vector<3, T, P>(const T x, const T y, const T z)
			: x(x), y(y), z(z) {}

		Vector<3, T, P>(const T vals[3])
			: x(vals[0]), y(vals[1]), z(vals[2]) {}

		Vector<3, T, P>(const Vec &vec)
			: x(vec.x), y(vec.y), z(vec.z) {}

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
			return *this;
		}

		inline Vec &operator-=(const Vec &vec)
		{
			x -= vec.x;
			y -= vec.y;
			z -= vec.z;
			return *this;
		}

		inline Vec &operator*=(const Vec &vec)
		{
			x *= vec.x;
			y *= vec.y;
			z *= vec.z;
			return *this;
		}

		inline Vec &operator*=(const T val)
		{
			x *= val;
			y *= val;
			z *= val;
			return *this;
		}

		inline Vec &operator/=(const Vec &vec)
		{
			x /= vec.x;
			y /= vec.y;
			z /= vec.z;
			return *this;
		}

		inline Vec &operator/=(const T val)
		{
			x /= val;
			y /= val;
			z /= val;
			return *this;
		}

		// -- Unary operators --

		inline Vec operator+() const
		{
			return *this;
		}

		inline Vec operator-() const
		{
			return Vector<3, T, P>(-x, -y, -z);
		}

		// -- Binary operators --

		inline Vec operator+(const Vec &vec) const
		{
			return Vector<3, T, P>(x + vec.x, y + vec.y, z + vec.z);
		}

		inline Vec operator-(const Vec &vec) const
		{
			return Vector<3, T, P>(x - vec.x, y - vec.y, z - vec.z);
		}

		inline Vec operator*(const Vec &vec) const
		{
			return Vector<3, T, P>(x * vec.x, y * vec.y, z * vec.z);
		}

		inline Vec operator*(const T val) const
		{
			return Vector<3, T, P>(x * val, y * val, z * val);
		}

		inline Vec operator/(const Vec &vec) const
		{
			return Vector<3, T, P>(x / vec.x, y / vec.y, z / vec.z);
		}

		inline Vec operator/(const T val) const
		{
			return Vector<3, T, P>(x / val, y / val, z / val);
		}

		// -- Boolean operators --

		inline bool operator==(const Vec &vec) const
		{
			return x == vec.x && y == vec.y && z == vec.z;
		}

		inline bool operator!=(const Vec &vec) const
		{
			return x != vec.x || y != vec.y || z != vec.z;
		}

		// -- Convertion operators --

		template<typename U, typename Q>
		inline explicit operator Vector<3, U, Q>() const
		{
			return Vector<3, U, Q>((U)x, (U)y, (U)z);
		}

		// -- Stream operators --

		friend std::ostream &operator<<(std::ostream &ostream, const Vec &vec)
		{
			return ostream << "(" << vec.x << ", " << vec.y << ", " << vec.z << ")";
		}

		// -- Getters --

		Vec Abs() const
		{
			return Vector<3, T, P>(Math::Abs(x), Math::Abs(y), Math::Abs(z));
		}

		Vec Cross(const Vec &vec) const
		{
			return Vector<3, T, P>(
				y * vec.z - z * vec.y,
				z * vec.x - x * vec.z,
				x * vec.y - y * vec.x
			);
		}

		P Distance(const Vec &vec) const
		{
			const P to_x = vec.x - x;
			const P to_y = vec.y - y;
			const P to_z = vec.z - z;
			return std::sqrt(to_x * to_x + to_y * to_y + to_z * to_z);
		}

		P DistanceSquared(const Vec &vec) const
		{
			const P to_x = vec.x - x;
			const P to_y = vec.y - y;
			const P to_z = vec.z - z;
			return to_x * to_x + to_y * to_y + to_z * to_z;
		}

		T Dot(const Vec &vec) const
		{
			return x * vec.x + y * vec.y + z * vec.z;
		}

		bool IsNormalized() const
		{
			return x * x + y * y + z * z == T(1);
		}

		P Length() const
		{
			return std::sqrt(x * x + y * y + z * z);
		}

		T LengthSquared() const
		{
			return x * x + y * y + z * z;
		}

		Vec Normalized() const
		{
			P length = LengthSquared();

			if (length == 0 || length == 1)
			{
				return *this;
			}

			length = std::sqrt(length);
			return Vector<3, T, P>(x / length, y / length, z / length);
		}

		Vec Sign() const
		{
			return Vector<3, T, P>(Math::Sign(x), Math::Sign(y), Math::Sign(z));
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
			return *this;
		}

		Vec &RotateAroundX(const P angle)
		{
			return RotateAroundX(std::cos(angle), std::sin(angle));
		}

		Vec &RotateAroundY(const P angle)
		{
			return RotateAroundY(std::cos(angle), std::sin(angle));
		}

		Vec &RotateAroundZ(const P angle)
		{
			return RotateAroundZ(std::cos(angle), std::sin(angle));
		}

		Vec &Rotate(const Vector<3, P, P> &rot, const EulerAngleOrder order = EulerAngleOrder::YXZ)
		{
			/*To refactor*/
			switch (order)
			{
				case EulerAngleOrder::XYZ:
					RotateAroundX(rot.x);
					RotateAroundY(rot.y);
					RotateAroundZ(rot.z);
					break;
				case EulerAngleOrder::XZY:
					RotateAroundX(rot.x);
					RotateAroundZ(rot.z);
					RotateAroundY(rot.y);
					break;
				case EulerAngleOrder::YXZ:
					RotateAroundY(rot.y);
					RotateAroundX(rot.x);
					RotateAroundZ(rot.z);
					break;
				case EulerAngleOrder::YZX:
					RotateAroundY(rot.y);
					RotateAroundZ(rot.z);
					RotateAroundX(rot.x);
					break;
				case EulerAngleOrder::ZXY:
					RotateAroundZ(rot.z);
					RotateAroundX(rot.x);
					RotateAroundY(rot.y);
					break;
				case EulerAngleOrder::ZYX:
					RotateAroundZ(rot.z);
					RotateAroundY(rot.y);
					RotateAroundX(rot.x);
					break;
				default:
					break;
			}

			return *this;
		}

		Vec Rotated(const Vector<3, P, P> &rot, const EulerAngleOrder order = EulerAngleOrder::YXZ) const
		{
			Vec res(*this);
			res.Rotate(rot, order);
			return res;
		}

		// -- Static getters --

		static Vec Zero()
		{
			return Vector<3, T, P>(T(0), T(0), T(0));
		}

		static Vec One()
		{
			return Vector<3, T, P>(T(1), T(1), T(1));
		}

		static Vec NegOne()
		{
			return Vector<3, T, P>(T(-1), T(-1), T(-1));
		}

		static Vec Left()
		{
			return Vector<3, T, P>(T(-1), T(0), T(0));
		}

		static Vec Right()
		{
			return Vector<3, T, P>(T(1), T(0), T(0));
		}

		static Vec Up()
		{
			return Vector<3, T, P>(T(0), T(1), T(0));
		}

		static Vec Down()
		{
			return Vector<3, T, P>(T(0), T(-1), T(0));
		}

		static Vec Forward()
		{
			return Vector<3, T, P>(T(0), T(0), T(1));
		}

		static Vec Backward()
		{
			return Vector<3, T, P>(T(0), T(0), T(-1));
		}

		// -- Static methods --

		static Vec Lerp(const Vec &lhs, const Vec &rhs, const P t)
		{
			return Vector<3, T, P>(
				Math::Lerp(lhs.x, rhs.x, t),
				Math::Lerp(lhs.y, rhs.y, t),
				Math::Lerp(lhs.z, rhs.z, t)
			);
		}

		static Vec LerpClamped(const Vec &lhs, const Vec &rhs, const P t)
		{
			return Vector<3, T, P>(
				Math::LerpClamped(lhs.x, rhs.x, t),
				Math::LerpClamped(lhs.y, rhs.y, t),
				Math::LerpClamped(lhs.z, rhs.z, t)
			);
		}
	};

	typedef Vector<3, float, float> Vector3;
	typedef Vector<3, int, float> Vector3I;
	typedef Vector<3, double, double> Vector3D;
}

