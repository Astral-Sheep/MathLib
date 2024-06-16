#pragma once

#include "Matrix.hpp"
#include <cassert>
#include <cmath>

namespace Math
{
	template<typename T>
	Matrix<4, 4, T, T> Ortho(const T left, const T right, const T bottom, const T top, const T near, const T far)
	{
		Matrix<4, 4, T, T> res(T(1));
		res[0][0] = T(2) / (right - left);
		res[1][1] = T(2) / (top - bottom);
		res[2][2] = T(1) / (far - near);
		res[3][0] = -(right + left) / (right - left);
		res[3][1] = -(top + bottom) / (top - bottom);
		res[3][2] = -near / (far - near);
		return res;
	}

	template<typename T>
	Matrix<4, 4, T, T> Frustum(const T left, const T right, const T bottom, const T top, const T near, const T far)
	{
		Matrix<4, 4, T, T> res(T(0));
		res[0][0] = (T(2) * near) / (right - left);
		res[1][1] = (T(2) * near) / (top - bottom);
		res[2][0] = -(right + left) / (right - left);
		res[2][1] = -(top + bottom) / (top - bottom);
		res[2][2] = far / (far - near);
		res[2][3] = T(1);
		res[3][2] = -(far * near) / (far - near);
		return res;
	}

	template<typename T>
	Matrix<4, 4, T, T> Perspective(const T fovY, const T aspect, const T near, const T far)
	{
		assert(abs(aspect - std::numeric_limits<T>::epsilon()) > T(0));

		const T tanHalfFovY = tan(fovY / T(2));

		Matrix<4, 4, T, T> res(T(0));
		res[0][0] = T(1) / (aspect * tanHalfFovY);
		res[1][1] = T(1) / (tanHalfFovY);
		res[2][2] = far / (far - near);
		res[2][3] = T(1);
		res[3][2] = -(far * near) / (far - near);
		return res;
	}

	template<typename T>
	Matrix<4, 4, T, T> PerspectiveFov(const T fov, const T width, const T height, const T near, const T far)
	{
		assert(width > T(0));
		assert(height > T(0));
		assert(fov > T(0));

		const T h = std::cos(T(0.5) * fov) / std::sin(T(0.5) * fov);
		const T w = h * height / width;

		Matrix<4, 4, T, T> res(T(0));
		res[0][0] = w;
		res[1][1] = h;
		res[2][2] = far / (far - near);
		res[2][3] = T(1);
		res[3][2] = -(far * near) / (far - near);
		return res;
	}

	template<typename T>
	Matrix<4, 4, T, T> InfinitePerspective(const T fovY, const T aspect, const T near)
	{
		const T range = std::tan(fovY / T(2)) * near;
		const T left = -range * aspect;
		const T right = range * aspect;
		const T bottom = -range;
		const T top = range;

		Matrix<4, 4, T, T> res(T(0));
		res[0][0] = (T(2) * near) / (right - left);
		res[1][1] = (T(2) * near) / (top - bottom);
		res[2][2] = T(1);
		res[2][3] = T(1);
		res[3][2] = -near;
		return res;
	}
}

