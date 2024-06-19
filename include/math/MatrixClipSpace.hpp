#pragma once

#include "Matrix.hpp"
#include "Math.hpp"
#include <cassert>
#include <cmath>

namespace Math
{
	template<typename T>
	Matrix<4, 4, T, T> Ortho(const T pLeft, const T pRight, const T pBottom, const T pTop, const T pNear, const T pFar)
	{
		Matrix<4, 4, T, T> lRes(T(1));
		lRes[0][0] = T(2) / (pRight - pLeft);
		lRes[1][1] = T(2) / (pTop - pBottom);
		lRes[2][2] = T(1) / (pFar - pNear);
		lRes[3][0] = -(pRight + pLeft) / (pRight - pLeft);
		lRes[3][1] = -(pTop + pBottom) / (pTop - pBottom);
		lRes[3][2] = -pNear / (pFar - pNear);
		return lRes;
	}

	template<typename T>
	Matrix<4, 4, T, T> Frustum(const T pLeft, const T pRight, const T pBottom, const T pTop, const T pNear, const T pFar)
	{
		Matrix<4, 4, T, T> lRes(T(0));
		lRes[0][0] = (T(2) * pNear) / (pRight - pLeft);
		lRes[1][1] = (T(2) * pNear) / (pTop - pBottom);
		lRes[2][0] = -(pRight + pLeft) / (pRight - pLeft);
		lRes[2][1] = -(pTop + pBottom) / (pTop - pBottom);
		lRes[2][2] = pFar / (pFar - pNear);
		lRes[2][3] = T(1);
		lRes[3][2] = -(pFar * pNear) / (pFar - pNear);
		return lRes;
	}

	template<typename T>
	Matrix<4, 4, T, T> Perspective(const T pFovY, const T pAspect, const T pNear, const T pFar)
	{
		assert(Abs(pAspect - std::numeric_limits<T>::epsilon()) > T(0));

		const T lTanHalfFovY = tan(pFovY / T(2));

		Matrix<4, 4, T, T> lRes(T(0));
		lRes[0][0] = T(1) / (pAspect * lTanHalfFovY);
		lRes[1][1] = T(1) / (lTanHalfFovY);
		lRes[2][2] = pFar / (pFar - pNear);
		lRes[2][3] = T(1);
		lRes[3][2] = -(pFar * pNear) / (pFar - pNear);
		return lRes;
	}

	template<typename T>
	Matrix<4, 4, T, T> PerspectiveFov(const T pFov, const T pWidth, const T pHeight, const T pNear, const T pFar)
	{
		assert(pWidth > T(0));
		assert(pHeight > T(0));
		assert(pFov > T(0));

		const T lH = std::cos(T(0.5) * pFov) / std::sin(T(0.5) * pFov);
		const T lW = lH * pHeight / pWidth;

		Matrix<4, 4, T, T> lRes(T(0));
		lRes[0][0] = lW;
		lRes[1][1] = lH;
		lRes[2][2] = pFar / (pFar - pNear);
		lRes[2][3] = T(1);
		lRes[3][2] = -(pFar * pNear) / (pFar - pNear);
		return lRes;
	}

	template<typename T>
	Matrix<4, 4, T, T> InfinitePerspective(const T pFovY, const T pAspect, const T pNear)
	{
		const T lRange = std::tan(pFovY / T(2)) * pNear;
		const T lLeft = -lRange * pAspect;
		const T lRight = lRange * pAspect;
		const T lBottom = -lRange;
		const T lTop = lRange;

		Matrix<4, 4, T, T> lRes(T(0));
		lRes[0][0] = (T(2) * pNear) / (lRight - lLeft);
		lRes[1][1] = (T(2) * pNear) / (lTop - lBottom);
		lRes[2][2] = T(1);
		lRes[2][3] = T(1);
		lRes[3][2] = -pNear;
		return lRes;
	}
}

