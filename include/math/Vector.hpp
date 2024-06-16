#pragma once

#include "Math.hpp"
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

		inline const T &operator[](const int index) const noexcept
		{
			return fields[index];
		}

		inline T &operator[](const int index) noexcept
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

			return *this;
		}

		Vec &operator-=(const Vec &vector)
		{
			for (int i = 0; i < S; i++)
			{
				fields[i] -= vector[i];
			}

			return *this;
		}

		Vec &operator*=(const Vec &vector)
		{
			for (int i = 0; i < S; i++)
			{
				fields[i] *= vector[i];
			}

			return *this;
		}

		Vec &operator*=(const T scalar)
		{
			for (int i = 0; i < S; i++)
			{
				fields[i] *= scalar;
			}

			return *this;
		}

		Vec &operator/=(const Vec &vector)
		{
			for (int i = 0; i < S; i++)
			{
				fields[i] /= vector[i];
			}

			return *this;
		}

		Vec &operator/=(const T scalar)
		{
			for (int i = 0; i < S; i++)
			{
				fields[i] /= scalar;
			}

			return *this;
		}

		// -- Unary operators --

		Vec operator+() const
		{
			return *this;
		}

		Vec operator-() const
		{
			Vec res;

			for (int i = 0; i < S; i++)
			{
				res[i] = -fields[i];
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

		// -- Convertion operators --

		template<typename U, typename Q>
		operator Vector<S, U, Q>() const
		{
			Vector<S, U, Q> res;

			for (int i = 0; i < S; i++)
			{
				res[i] = (U)fields[i];
			}

			return res;
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
			return from + (to - from) * Math::Clamp01(t);
		}

		// -- Getters --

		Vec Abs() const
		{
			Vec res;

			for (int i = 0; i < S; i++)
			{
				res[i] = Math::Abs(fields[i]);
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
				res[i] = Math::Sign(fields[i]);
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

	template<unsigned int S>
	using VectorI = Vector<S, int, float>;

	template<unsigned int S>
	using VectorF = Vector<S, float, float>;

	template<unsigned int S>
	using VectorD = Vector<S, double, double>;
}

