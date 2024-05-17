#pragma once

#include <ostream>

namespace Math
{
	template<unsigned int S, typename T, typename P>
	struct Vector
	{
	private:
		typedef Vector<S, T, P> Vec;
		T fields[S];

	public:
		// -- Constructors --

		Vector()
		{
			for (int i = 0; i < S; i++)
			{
				fields[i] = T(0);
			}
		}

		explicit Vector(const T value)
		{
			for (int i = 0; i < S; i++)
			{
				fields[i] = value;
			}
		}

		Vector(const T fields[S])
		{
			for (int i = 0; i < S; i++)
			{
				this->fields[i] = fields[i];
			}
		}

		Vector(const Vector &vector)
		{
			for (int i = 0; i < S; i++)
			{
				fields[i] = vector[i];
			}
		}

		// -- Accesses --

		const T &operator[](const int index) const noexcept
		{
			return fields[index];
		}

		T &operator[](const int index) noexcept
		{
			return fields[index];
		}

		// -- Unary arithmetic operators --

		Vec &operator+=(const Vec &vector)
		{
			for (int i = 0; i < S; i++)
			{
				fields[i] += vector[i];
			}
		}

		Vec &operator-=(const Vec &vector)
		{
			for (int i = 0; i < S; i++)
			{
				fields[i] -= vector[i];
			}
		}

		Vec &operator*=(const Vec &vector)
		{
			for (int i = 0; i < S; i++)
			{
				fields[i] *= vector[i];
			}
		}

		Vec &operator*=(const T scalar)
		{
			for (int i = 0; i < S; i++)
			{
				fields[i] *= scalar;
			}
		}

		Vec &operator/=(const Vec &vector)
		{
			for (int i = 0; i < S; i++)
			{
				fields[i] /= vector[i];
			}
		}

		Vec &operator/=(const T scalar)
		{
			for (int i = 0; i < S; i++)
			{
				fields[i] /= scalar;
			}
		}

		// -- Unary operators --

		Vec operator-() const
		{
			Vec res;

			for (int i = 0; i < S; i++)
			{
				res[i] = fields[i];
			}

			return res;
		}

		// -- Binary operators --

		Vec operator+(const Vec &vector) const
		{
			Vec res;

			for (int i = 0; i < S; i++)
			{
				res[i] = fields[i] + vector[i];
			}

			return res;
		}

		Vec operator-(const Vec &vector) const
		{
			Vec res;

			for (int i = 0; i < S; i++)
			{
				res[i] = fields[i] - vector[i];
			}

			return res;
		}

		Vec operator*(const Vec &vector) const
		{
			Vec res;

			for (int i = 0; i < S; i++)
			{
				res[i] = fields[i] * vector[i];
			}

			return res;
		}

		Vec operator*(const T scalar) const
		{
			Vec res;

			for (int i = 0; i < S; i++)
			{
				res[i] = fields[i] * scalar;
			}

			return res;
		}

		Vec operator/(const Vec &vector) const
		{
			Vec res;

			for (int i = 0; i < S; i++)
			{
				res[i] = fields[i] / vector[i];
			}

			return res;
		}

		Vec operator/(const T scalar) const
		{
			Vec res;

			for (int i = 0; i < S; i++)
			{
				res[i] = fields[i] / scalar;
			}

			return res;
		}

		// -- Boolean operators --

		bool operator==(const Vec &vector) const
		{
			for (int i = 0; i < S; i++)
			{
				if (fields[i] != vector[i])
				{
					return false;
				}
			}

			return true;
		}

		bool operator!=(const Vec &vector) const
		{
			for (int i = 0; i < S; i++)
			{
				if (fields[i] == vector[i])
				{
					return false;
				}
			}

			return true;
		}

		// -- Stream operators --

		friend std::ostream &operator<<(std::ostream &ostream, const Vec &vector)
		{
			for (int i = 0; i < S; i++)
			{
				ostream << vector[i];

				if (i < S - 1)
				{
					ostream << ',';
				}
			}

			return ostream;
		}

		// -- Static methods --

		static Vec Lerp(const Vec &from, const Vec &to, const float t)
		{
			return from + (to - from) * t;
		}

		static Vec LerpClamped(const Vec &from, const Vec &to, const float t)
		{
			return from + (to - from) * Clamp01(t);
		}

		// -- Getters --

		Vec Abs() const
		{
			Vec res;

			for (int i = 0; i < S; i++)
			{
				res[i] = Abs(fields[i]);
			}

			return res;
		}

		P Distance(const Vec &vector) const
		{
			P sqrdist = P(0);

			for (int i = 0; i < S; i++)
			{
				sqrdist += (vector[i] - fields[i]) * (vector[i] - fields[i]);
			}

			return std::sqrt(sqrdist);
		}

		P DistanceSquared(const Vec &vector) const
		{
			P sqrdist = P(0);

			for (int i = 0; i < S; i++)
			{
				sqrdist += (vector[i] - fields[i]) * (vector[i] - fields[i]);
			}

			return sqrdist;
		}

		P Dot(const Vec &vector) const
		{
			P dot = P(0);

			for (int i = 0; i < S; i++)
			{
				dot = fields[i] * vector[i];
			}

			return dot;
		}

		inline bool IsNormalized() const
		{
			return LengthSquared() == P(1);
		}

		P Length() const
		{
			P sqrLength = P(0);

			for (int i = 0; i < S; i++)
			{
				sqrLength += fields[i] * fields[i];
			}

			return std::sqrt(sqrLength);
		}

		P LengthSquared() const
		{
			P sqrLength = P(0);

			for (int i = 0; i < S; i++)
			{
				sqrLength += fields[i] * fields[i];
			}

			return sqrLength;
		}

		Vec Normalized() const
		{
			P length = LengthSquared();

			if (length == 0 || length == 1)
			{
				return *this;
			}

			Vec res;
			length = std::sqrt(length);

			for (int i = 0; i < S; i++)
			{
				res[i] = fields[i] / length;
			}

			return res;
		}

		Vec Sign() const
		{
			Vec res;

			for (int i = 0; i < S; i++)
			{
				res[i] = Sign(fields[i]);
			}

			return res;
		}

		inline size_t SizeOfField() const
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

			for (int i = 0; i < S; i++)
			{
				fields[i] /= length;
			}

			return *this;
		}
	};

	typedef Vector<2, float, float> Vector2;
	typedef Vector<2, int, float> Vector2I;
	typedef Vector<2, double, double> Vector2D;

	typedef Vector<3, float, float> Vector3;
	typedef Vector<3, int, float> Vector3I;
	typedef Vector<3, double, double> Vector3D;

	typedef Vector<4, float, float> Vector4;
	typedef Vector<4, int, float> Vector4I;
	typedef Vector<4, double, double> Vector4D;

#define VECTOR_BODY(S, T, P)\
	private:\
		typedef Vector<S, T, P> Vec##S;\
	\
	public:\
		/* -- Constructors -- */\
		Vector<S, T, P>();\
		Vector<S, T, P>(const T value);\
		Vector<S, T, P>(const T fields[S]);\
		Vector<S, T, P>(const Vec##S &vector);\
		\
		/* -- Accesses -- */\
		const T &operator[](const int index) const noexcept;\
		T &operator[](const int index) noexcept;\
		\
		/* -- Unary arithmetic operators -- */\
		Vec##S &operator+=(const Vec##S &vector);\
		Vec##S &operator-=(const Vec##S &vector);\
		Vec##S &operator*=(const Vec##S &vector);\
		Vec##S &operator*=(const T scalar);\
		Vec##S &operator/=(const Vec##S &vector);\
		Vec##S &operator/=(const T scalar);\
		\
		/* -- Unary operators -- */\
		Vec##S operator-() const;\
		\
		/* -- Binary operators -- */\
		Vec##S operator+(const Vec##S &vector) const;\
		Vec##S operator-(const Vec##S &vector) const;\
		Vec##S operator*(const Vec##S &vector) const;\
		Vec##S operator*(const T scalar) const;\
		Vec##S operator/(const Vec##S &vector) const;\
		Vec##S operator/(const T scalar) const;\
		\
		/* -- Boolean operators -- */\
		bool operator==(const Vec##S &vector) const;\
		bool operator!=(const Vec##S &vector) const;\
		\
		/* -- Stream operators -- */\
		friend std::ostream &operator<<(std::ostream &ostream, const Vec##S &vector);\
		\
		/* -- Static methods -- */\
		static Vec##S Lerp(const Vec##S &from, const Vec##S &to, const P t);\
		static Vec##S LerpClamped(const Vec##S &from, const Vec##S &to, const P t);\
		\
		/* -- Getters -- */\
		Vec##S Abs() const;\
		P Distance(const Vec##S &vector) const;\
		P DistanceSquared(const Vec##S &vector) const;\
		P Dot(const Vec##S &vector) const;\
		bool IsNormalized() const;\
		P Length() const;\
		P LengthSquared() const;\
		Vec##S Normalized() const;\
		Vec##S Sign() const;\
		size_t SizeOfField() const;\
		\
		/* -- Transformations -- */\
		Vec##S &Normalize();
}

