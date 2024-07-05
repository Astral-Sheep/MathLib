#pragma once

#include "Matrix.hpp"
#include "Vector.hpp"
#include "Vector3.hpp"
#include "Quaternion.hpp"
#include "Matrix4.hpp"
#include <cmath>

namespace Math
{
	template<typename T>
	Matrix<4, 4, T, float> MakeTransform(const Vector<3, T, float> &pTrans, const Vector<3, T, float> &pRot, const Vector<3, T, float> &pScale)
	{
		Matrix<4, 4, T, float> lTransMat(
			1.0f, 0.0f, 0.0f, pTrans.x,
			0.0f, 1.0f, 0.0f, pTrans.y,
			0.0f, 0.0f, 1.0f, pTrans.z,
			0.0f, 0.0f, 0.0f, 1.0f
		);

		Matrix<4, 4, T, float> lRotMat = Rotate(
			Rotate(
				Rotate(Matrix<4, 4, T, float>(T(1)), pRot.z, Vector<3, T, float>::Up()),
				pRot.y,
				Vector<3, T, float>::Right()
			),
			pRot.x,
			Vector<3, T, float>::Forward()
		);

		Matrix<4, 4, T, float> lScaleMat(
			pScale.x, 0.0f, 0.0f, 0.0f,
			0.0f, pScale.y, 0.0f, 0.0f,
			0.0f, 0.0f, pScale.z, 0.0f,
			0.0f, 0.0f, 0.0f, 1.0f
		);

		return lTransMat * lRotMat * lScaleMat;
	}

	template<typename T>
	Matrix<4, 4, T, float> Translate(const Matrix<4, 4, T, float> &pMat, const Vector<3, T, float> &pVec)
	{
		Matrix<4, 4, T, float> lRes(pMat);
		lRes[3] = pMat[0] * pVec[0] + pMat[1] * pVec[1] + pMat[2] * pVec[2] + pMat[3];
		return lRes;
	}

	template<typename T>
	Matrix<4, 4, T, float> Rotate(const Matrix<4, 4, T, float> &pMat, const T pAngle, const Vector<3, T, float> &pVec)
	{
		const T lCos = std::cos(pAngle);
		const T lSin = std::sin(pAngle);

		const Vector<3, T, float> lAxis = pVec.Normalized();
		const Vector<3, T, float> lTemp(lAxis * (T(1) - lCos));

		Matrix<4, 4, T, float> rotate;
		rotate[0][0] = lCos + lTemp[0] * lAxis[0];
		rotate[0][1] = lTemp[0] * lAxis[1] + lSin * lAxis[2];
		rotate[0][2] = lTemp[0] * lAxis[2] - lSin * lAxis[1];

		rotate[1][0] = lTemp[1] * lAxis[0] - lSin * lAxis[2];
		rotate[1][1] = lCos + lTemp[1] * lAxis[1];
		rotate[1][2] = lTemp[1] * lAxis[2] + lSin * lAxis[0];

		rotate[2][0] = lTemp[2] * lAxis[0] + lSin * lAxis[1];
		rotate[2][1] = lTemp[2] * lAxis[1] - lSin * lAxis[0];
		rotate[2][2] = lCos + lTemp[2] * lAxis[2];

		Matrix<4, 4, T, float> lRes;
		lRes[0] = pMat[0] * rotate[0][0] + pMat[1] * rotate[0][1] + pMat[2] * rotate[0][2];
		lRes[1] = pMat[0] * rotate[1][0] + pMat[1] * rotate[1][1] + pMat[2] * rotate[1][2];
		lRes[2] = pMat[0] * rotate[2][0] + pMat[1] * rotate[2][1] + pMat[2] * rotate[2][2];
		lRes[3] = pMat[3];
		return lRes;
	}

	template<typename T>
	Matrix<4, 4, T, float> Scale(const Matrix<4, 4, T, float> &pMat, const Vector<3, T, float> &pVec)
	{
		Matrix<4, 4, T, float> lRes;
		lRes[0] = pMat[0] * pVec[0];
		lRes[1] = pMat[1] * pVec[1];
		lRes[2] = pMat[2] * pVec[2];
		lRes[3] = pMat[3];
		return lRes;
	}

	template<typename T>
	Matrix<4, 4, T, float> LookAt(const Vector<3, T, float> &pEye, const Vector<3, T, float> &pCenter, const Vector<3, T, float> &pUp)
	{
		const Vector<3, T, float> f = (pCenter - pEye).Normalize();
		const Vector<3, T, float> s = f.Cross(pUp).Normalize();
		const Vector<3, T, float> u = s.Cross(f);

		Matrix<4, 4, T, float> res(T(1));
		res[0][0] = s.x;
		res[1][0] = s.y;
		res[2][0] = s.z;
		res[0][1] = u.x;
		res[1][1] = u.y;
		res[2][1] = u.z;
		res[0][2] =-f.x;
		res[1][2] =-f.y;
		res[2][2] =-f.z;
		res[3][0] =-s.Dot(pEye);
		res[3][1] =-u.Dot(pEye);
		res[3][2] = f.Dot(pEye);
		return res;
	}

	template<typename T>
	Vector<3, T, float> GetTranslation(const Matrix<4, 4, T, float> &pMat)
	{
		Vector<3, T, float> lRes;
		lRes[0] = pMat[3][0];
		lRes[1] = pMat[3][1];
		lRes[2] = pMat[3][2];
		return lRes;
	}

	template<typename T>
	Vector<3, T, float> GetEulerRotation(const Matrix<4, 4, T, float> &pMat)
	{
		const Vector<3, T, float> lMat0 = pMat[0].Normalized();
		const Vector<3, T, float> lMat1 = pMat[1].Normalized();
		const Vector<3, T, float> lMat2 = pMat[2].Normalized();

		Vector<3, T, float> lRes;
		lRes[0] = std::atan2(lMat1[2], lMat2[2]);
		lRes[1] = std::atan2(-lMat0[2], std::sqrt(lMat1[2] * lMat1[2] + lMat2[2] * lMat2[2]));
		lRes[2] = std::atan2(lMat0[1], lMat0[0]);
		return lRes;
	}

	template<typename T>
	Quaternion GetQuatRotation(const Matrix<4, 4, T, float> &pMat)
	{
		const Vector<3, T, float> lMat0 = pMat[0].Normalized();
		const Vector<3, T, float> lMat1 = pMat[1].Normalized();
		const Vector<3, T, float> lMat2 = pMat[2].Normalized();

		return Quaternion::FromEulerAngles(
			std::atan2(lMat1[2], lMat2[2]),
			std::atan2(-lMat0[2], std::sqrt(lMat1[2] * lMat1[2] + lMat2[2] * lMat2[2])),
			std::atan2(lMat0[1], lMat0[0])
		);
	}

	template<typename T>
	Vector<3, T, float> GetScale(const Matrix<4, 4, T, float> &pMat)
	{
		Vector<3, T, float> lRes;
		lRes[0] = pMat[0].Length();
		lRes[1] = pMat[1].Length();
		lRes[2] = pMat[2].Length();
		return lRes;
	}

	template<typename T>
	void SetTranslation(Matrix<4, 4, T, float> &pMat, const Vector<3, T, float> &pTrans)
	{
		memcpy(&pMat[3], &pTrans, 3 * sizeof(T));
	}

	template<typename T>
	void SetRotation(Matrix<4, 4, T, float> &pMat, const Vector<3, T, float> &pRot)
	{
		const Vector3 lScale = GetScale(pMat);
		const Vector3 lTrans = GetTranslation(pMat);

		pMat = MakeRotation(Vector<3, T, float>::Right(), pRot.x) * MakeRotation(Vector<3, T, float>::Up(), pRot.y) * MakeRotation(Vector<3, T, float>::Forward(), pRot.z);
		pMat[0] *= lScale.x;
		pMat[1] *= lScale.y;
		pMat[2] *= lScale.z;
		pMat[3][0] = lTrans.x;
		pMat[3][1] = lTrans.y;
		pMat[3][2] = lTrans.z;
	}

	template<typename T>
	void SetRotation(Matrix<4, 4, T, float> &pMat, const QuaternionT<T, float> &pQuat)
	{
		SetRotation(pMat, pQuat.ToEulerAngles());
	}

	template<typename T>
	void SetScale(Matrix<4, 4, T, float> &pMat, const Vector<3, T, float> &pScale)
	{
		pMat[0] = pMat[0].Normalized() * pScale.x;
		pMat[1] = pMat[1].Normalized() * pScale.y;
		pMat[2] = pMat[2].Normalized() * pScale.z;
	}
}

