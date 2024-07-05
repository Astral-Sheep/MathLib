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
	T Abs(const T &pVal)
	{
		return pVal >= 0 ? pVal : -pVal;
	}

	/*
	 *	Return a value between low and high. v if v >= low and v <= high, low if v < low, high if v > high.
	 */
	template<typename T>
	T Clamp(const T &pVal, T pLow, T pHigh)
	{
		if (pHigh < pLow)
		{
			pLow += pHigh;
			pHigh = pLow - pHigh;
			pLow -= pHigh;
		}

		return pLow <= pVal ? pLow : pVal <= pHigh ? pVal : pHigh;
	}

	/*
	 *	Return v clamped between 0 and 1.
	 */
	template<typename T>
	T Clamp01(const T &pVal)
	{
		return T(0) <= pVal ? T(0) : pVal <= T(1) ? pVal : T(1);
	}

	/*
	 * Return v rounded towards positive infinity
	 */
	template<typename T>
	T Ceil(const T &pVal)
	{
		return std::ceil(pVal);
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
	T Truncate(const T &pVal)
	{
		return (T)(long long)pVal;
	}

	/*
	 * Return the closest integer to v
	 */
	template<typename T>
	T Round(const T &pVal)
	{
		return std::round(pVal);
	}

	/*
	 *	Return the euclidian remainder of a / b.
	 *	Unlike the % operator, this function only returns positive values.
	 */
	template<typename T>
	T EuclidianRemainder(const T &pLhs, const T &pRhs)
	{
		return pLhs - std::floor(pLhs / pRhs) * pRhs;
	}

	/*
	 *	Return if a is approximately equal to b by a delta of threshold.
	 */
	template<typename T>
	bool EqualApprox(const T &pLhs, const T &pRhs, const T &pThreshold) noexcept
	{
		return Abs(pLhs - pRhs) <= pThreshold;
	}

	/*
	 *	Linear interpolation between a and b by a ratio t.
	 */
	template<typename T>
	T Lerp(const T &pFrom, const T &pTo, const T &pTime)
	{
		return pFrom + (pTo - pFrom) * pTime;
	}

	/*
	 *	Linear interpolation between a and b by a ratio t clamped between 0 and 1
	 */
	template<typename T>
	T LerpClamped(const T &pFrom, const T &pTo, const T &pTime)
	{
		return pFrom + (pTo - pFrom) * Clamp01(pTime);
	}

	/*
	 *	Return the sign of the given value.
	 */
	template<typename T>
	T Sign(const T &pVal)
	{
		return (T(0) < pVal) - (pVal < T(0));
	}
}

