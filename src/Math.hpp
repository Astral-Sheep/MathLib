#pragma once

namespace Math
{
	constexpr float PI = 3.1415926535f;
	constexpr float RAD2DEG = 180.0f / PI;
	constexpr float DEG2RAD = PI / 180.0f;

	/*
	 *	Return a value between lhs and rhs. a if a >= lhs and a <= rhs, lhs if a < lhs, rhs if a > rhs.
	 */
	int Clamp(const int a, const int lhs, const int rhs);

	/*
	 *	Return a value between lhs and rhs. a if a >= lhs and a <= rhs, lhs if a < lhs, rhs if a > rhs.
	 */
	float Clamp(const float a, const float lhs, const float rhs);

	/*
	 *	Return a value between lhs and rhs. a if a >= lhs and a <= rhs, lhs if a < lhs, rhs if a > rhs.
	 */
	double Clamp(const double a, const double lhs, const double rhs);

	/*
	 *	Return a clamped between 0 and 1.
	 */
	int Clamp01(const int a);

	/*
	 *	Return a clamped between 0 and 1.
	 */
	float Clamp01(const float a);

	/*
	 *	Return a clamped between 0 and 1.
	 */
	double Clamp01(const double a);

	/*
	 *	Return the euclidian remainder of a / b.
	 *	Unlike the % operator, this function only returns positive values.
	 */
	int EuclidianRemainder(const int a, const int b);

	/*
	 *	Return the euclidian remainder of a / b.
	 *	Unlike the % operator, this function only returns positive values.
	 */
	float EuclidianRemainder(const float a, const float b);

	/*
	 *	Return the euclidian remainder of a / b.
	 *	Unlike the % operator, this function only returns positive values.
	 */
	double EuclidianRemainder(const double a, const double b);

	/*
	 *	Return if a is approximately equal to be by a delta of threshold.
	 */
	bool EqualApprox(const float a, const float b, const float threshold) noexcept;

	/*
	 *	Return if a is approximately equal to be by a delta of threshold.
	 */
	bool EqualApprox(const double a, const double b, const double threshold) noexcept;

	/*
	 *	Linear interpolation between a and b by a ratio t.
	 */
	int Lerp(const int from, const int to, const float t);

	/*
	 *	Linear interpolation between a and b by a ratio t.
	 */
	float Lerp(const float from, const float to, const float t);

	/*
	 *	Linear interpolation between a and b by a ratio t.
	 */
	double Lerp(const double from, const double to, const float t);

	/*
	 *	Linear interpolation between a and b by a ratio t clamped between 0 and 1
	 */
	int LerpClamped(const int from, const int to, const float t);

	/*
	 *	Linear interpolation between a and b by a ratio t clamped between 0 and 1
	 */
	float LerpClamped(const float from, const float to, const float t);

	/*
	 *	Linear interpolation between a and b by a ratio t clamped between 0 and 1
	 */
	double LerpClamped(const double from, const double to, const float t);

	/*
	 *	Return the sign of the given value.
	 */
	template<typename T>
	T Sign(const T &value)
	{
		return (T(0) < value) - (value < T(0));
	}
}

