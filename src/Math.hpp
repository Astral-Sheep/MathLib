#pragma once

#include <cmath>

namespace Math
{
	constexpr float PI = 3.1415926535f;
	constexpr float TAU = 6.2831853071f;
	constexpr float RAD2DEG = 180.0f / PI;
	constexpr float DEG2RAD = PI / 180.0f;

	/*
	 * Return the absolute value of v
	 */
	template<typename T>
	T Abs(const T &v)
	{
		return v >= 0 ? v : -v;
	}

	/*
	 *	Return a value between low and high. v if v >= low and v <= high, low if v < low, high if v > high.
	 */
	template<typename T>
	T Clamp(const T &v, T low, T high)
	{
		if (high < low)
		{
			low += high;
			high = low - high;
			low -= high;
		}

		return low <= v ? low : v <= high ? v : high;
	}

	/*
	 *	Return v clamped between 0 and 1.
	 */
	template<typename T>
	T Clamp01(const T &v)
	{
		return T(0) <= v ? T(0) : v <= T(1) ? v : T(1);
	}

	/*
	 * Return v rounded towards positive infinity
	 */
	template<typename T>
	T Ceil(const T &v)
	{
		return std::ceil(v);
	}

	/*
	 * Return v rounded towards negative infinity
	 */
	template<typename T>
	T Floor(const T &v)
	{
		return std::floor(v);
	}

	/*
	 * Return v without its decimal part
	 */
	template<typename T>
	T Truncate(const T &v)
	{
		return (T)(long long)v;
	}

	/*
	 * Return the closest integer to v
	 */
	template<typename T>
	T Round(const T &v)
	{
		return std::round(v);
	}

	/*
	 *	Return the euclidian remainder of a / b.
	 *	Unlike the % operator, this function only returns positive values.
	 */
	template<typename T>
	T EuclidianRemainder(const T &a, const T &b)
	{
		return a - std::floor(a / b) * b;
	}

	/*
	 *	Return if a is approximately equal to b by a delta of threshold.
	 */
	template<typename T>
	bool EqualApprox(const T &a, const T &b, const T &threshold) noexcept
	{
		return Abs(a - b) <= threshold;
	}

	/*
	 *	Linear interpolation between a and b by a ratio t.
	 */
	template<typename T>
	T Lerp(const T &from, const T &to, const T &t)
	{
		return from + (to - from) * t;
	}

	/*
	 *	Linear interpolation between a and b by a ratio t clamped between 0 and 1
	 */
	template<typename T>
	T LerpClamped(const T &from, const T &to, const T &t)
	{
		return from + (to - from) * Clamp01(t);
	}

	/*
	 *	Return the sign of the given value.
	 */
	template<typename T>
	T Sign(const T &value)
	{
		return (T(0) < value) - (value < T(0));
	}
}

