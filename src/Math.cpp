#include "Math.hpp"

#include <algorithm>
#include <limits.h>
#include <cmath>

namespace Math
{
	int Clamp(const int a, const int lhs, const int rhs)
	{
		return std::clamp(a, lhs, rhs);
	}

	float Clamp(const float a, const float lhs, const float rhs)
	{
		return std::clamp(a, lhs, rhs);
	}

	double Clamp(const double a, const double lhs, const double rhs)
	{
		return std::clamp(a, lhs, rhs);
	}

	int Clamp01(const int a)
	{
		return std::clamp(a, 0, 1);
	}

	float Clamp01(const float a)
	{
		return std::clamp(a, 0.0f, 1.0f);
	}

	double Clamp01(const double a)
	{
		return std::clamp(a, 0.0, 1.0);
	}

#define DEFAULT_EUCLIDIAN_REMAINDER(a, b) (a) - floor((a) / (b)) * (b)

	int EuclidianRemainder(const int a, const int b)
	{
#if (-1 >> 1) == -1
		#define OFFSET ((a % b >> (sizeof(int) * CHAR_BIT - 1)) & b)
#else
		#define OFFSET ((a % b < 0) * b)
#endif //ARITHMETIC SHIFT CHECK

		return a % b + OFFSET;
	}

	float EuclidianRemainder(const float a, const float b)
	{
		return DEFAULT_EUCLIDIAN_REMAINDER(a, b);
	}

	double EuclidianRemainder(const double a, const double b)
	{
		return DEFAULT_EUCLIDIAN_REMAINDER(a, b);
	}

#define EQUAL_APPROX(a, b, threshold) std::abs((a) - (b)) <= (threshold)

	bool EqualApprox(const float a, const float b, const float threshold) noexcept
	{
		return EQUAL_APPROX(a, b, threshold);
	}

	bool EqualApprox(const double a, const double b, const double threshold) noexcept
	{
		return EQUAL_APPROX(a, b, threshold);
	}

#define LERP(from, to, t) (from) + ((to) - (from)) * (t)

	int Lerp(const int from, const int to, const float t)
	{
		return LERP(from, to, t);
	}

	float Lerp(const float from, const float to, const float t)
	{
		return LERP(from, to, t);
	}

	double Lerp(const double from, const double to, const float t)
	{
		return LERP(from, to, t);
	}

	int LerpClamped(const int from, const int to, const float t)
	{
		return LERP(from, to, Clamp(t, 0.0f, 1.0f));
	}

	float LerpClamped(const float from, const float to, const float t)
	{
		return LERP(from, to, Clamp(t, 0.0f, 1.0f));
	}

	double LerpClamped(const double from, const double to, const float t)
	{
		return LERP(from, to, Clamp(t, 0.0f, 1.0f));
	}
}

