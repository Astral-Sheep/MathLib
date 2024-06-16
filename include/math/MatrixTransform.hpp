#pragma once

#include "Matrix.hpp"
#include "Vector.hpp"
#include "Vector3.hpp"
#include "Quaternion.hpp"
#include "Math.hpp"
#include <cfloat>
#include <cmath>
#include <limits>

namespace Math
{
	template<typename T>
	Matrix<4, 4, T, float> MakeTransform(const Vector<3, T, float> &trans, const Vector<3, T, float> &rot, const Vector<3, T, float> &scale)
	{
		Matrix<4, 4, T, float> transMat(
			1.0f, 0.0f, 0.0f, trans.x,
			0.0f, 1.0f, 0.0f, trans.y,
			0.0f, 0.0f, 1.0f, trans.z,
			0.0f, 0.0f, 0.0f, 1.0f
		);

		Matrix<4, 4, T, float> rotMat(1.0f);
		Rotate(
			Rotate(
				Rotate(rotMat, rot.z, Vector<3, T, float>::Up()),
				rot.y,
				Vector<3, T, float>::Right()
			),
			rot.x,
			Vector<3, T, float>::Forward()
		);

		Matrix<4, 4, T, float> scaleMat(
			scale.x, 0.0f, 0.0f, 0.0f,
			0.0f, scale.y, 0.0f, 0.0f,
			0.0f, 0.0f, scale.z, 0.0f,
			0.0f, 0.0f, 0.0f, 1.0f
		);

		return transMat * rotMat * scaleMat;
	}

	/* template<typename T> */
	/* Matrix<4, 4, T> MakeTransform(const Vector<3, T, float> &trans, const Vector<3, T, float> &rot, const Vector<3, T, float> &scale) */
	/* { */
	/* 	// Apply rotation */
	/* 	const T right[3] = { T(1), T(0), T(0) }; */
	/* 	const T up[3] = { T(0), T(1), T(0) }; */
	/* 	const T forward[3] = { T(0), T(0), T(1) }; */
	/* 	Matrix<4, 4, T> res = MakeRotation(Vector<3, T, float>(right), rot.x) * */
	/* 		MakeRotation(Vector<3, T, float>(up), rot.y) * */
	/* 		MakeRotation(Vector<3, T, float>(forward), rot.z); */

	/* 	// Apply scale */
	/* 	for (int i = 0; i < 3; i++) */
	/* 	{ */
	/* 		res[i] *= Abs(scale[i]) < std::numeric_limits<T>::epsilon() ? T(0.001) : scale[i]; */
	/* 	} */

	/* 	// Apply translation */
	/* 	memcpy(&res[3], &trans, 3 * sizeof(T)); */
	/* 	return res; */
	/* } */

	/* template<typename T> */
	/* Matrix<4, 4, T> MakeTransform(const Vector<3, T, float> &trans, const QuaternionT<T, float> &rot, const Vector<3, T, float> &scale) */
	/* { */
	/* 	// Apply rotation */
	/* 	const T right[3] = { T(1), T(0), T(0) }; */
	/* 	const T up[3] = { T(0), T(1), T(0) }; */
	/* 	const T forward[3] = { T(0), T(0), T(1) }; */
	/* 	Vector<3, T, float> euler = rot.ToEulerAngles(); */
	/* 	Matrix<4, 4, T> res = MakeRotation(Vector<3, T, float>(right), euler.x) * */
	/* 		MakeRotation(Vector<3, T, float>(up), euler.y) * */
	/* 		MakeRotation(Vector<3, T, float>(forward), euler.z); */

	/* 	// Apply scale */
	/* 	for (int i = 0; i < 3; i++) */
	/* 	{ */
	/* 		res[i] *= Abs(scale[i]) < std::numeric_limits<T>::epsilon() ? T(0.001) : scale[i]; */
	/* 	} */

	/* 	// Apply translation */
	/* 	memcpy(&res[3], &trans, 3 * sizeof(T)); */
	/* 	return res; */
	/* } */

	/* template<typename T> */
	/* Matrix<4, 4, T> MakeRotation(const Vector<3, T, float> &axis, const T angle) */
	/* { */
	/* 	float sqrlength = axis.LengthSquared(); */

	/* 	if (axis.LengthSquared() < FLT_EPSILON) */
	/* 	{ */
	/* 		return Matrix<4, 4, T>(1.0f); */
	/* 	} */

	/* 	Vector<3, T, float> norm = axis.IsNormalized() ? axis : axis / std::sqrt(sqrlength); */
	/* 	T sin = std::sin(angle); */
	/* 	T cos = std::cos(angle); */
	/* 	T k = 1.0f - cos; */

	/* 	T xx = norm.x * norm.x * k + cos; */
	/* 	T yy = norm.y * norm.y * k + cos; */
	/* 	T zz = norm.z * norm.z * k + cos; */
	/* 	T xy = norm.x * norm.y * k; */
	/* 	T yz = norm.y * norm.z * k; */
	/* 	T zx = norm.z * norm.x * k; */
	/* 	T xs = norm.x * sin; */
	/* 	T ys = norm.y * sin; */
	/* 	T zs = norm.z * sin; */

	/* 	Matrix<4, 4, T> res(1.0f); */
	/* 	res[0][0] = xx; */
	/* 	res[0][1] = xy + zs; */
	/* 	res[0][2] = zx - ys; */
	/* 	res[1][0] = xy - zs; */
	/* 	res[1][1] = yy; */
	/* 	res[1][2] = yz + xs; */
	/* 	res[0][2] = zx + ys; */
	/* 	res[1][2] = yz - xs; */
	/* 	res[2][2] = zz; */
	/* 	return res; */
	/* } */

	template<typename T>
	Matrix<4, 4, T, float> Translate(const Matrix<4, 4, T, float> &mat, const Vector<3, T, float> &vec)
	{
		Matrix<4, 4, T, float> res(mat);
		res[3] = mat[0] * vec[0] + mat[1] * vec[1] + mat[2] * vec[2] + mat[3];
		return res;
	}

	template<typename T>
	Matrix<4, 4, T, float> Rotate(const Matrix<4, 4, T, float> &mat, const T angle, const Vector<3, T, float> &vec)
	{
		const T cos = std::cos(angle);
		const T sin = std::sin(angle);

		Vector<3, T, float> axis = vec.Normalized();
		Vector<3, T, float> temp(axis * (T(1) - cos));

		Matrix<4, 4, T, float> rotate;
		rotate[0][0] = cos + temp[0] * axis[0];
		rotate[0][1] = temp[0] * axis[1] + sin * axis[2];
		rotate[0][2] = temp[0] * axis[2] - sin * axis[1];

		rotate[1][0] = temp[1] * axis[0] - sin * axis[2];
		rotate[1][1] = cos + temp[1] * axis[1];
		rotate[1][2] = temp[1] * axis[2] + sin * axis[0];

		rotate[2][0] = temp[2] * axis[0] + sin * axis[1];
		rotate[2][1] = temp[2] * axis[1] - sin * axis[0];
		rotate[2][2] = cos + temp[2] * axis[2];

		Matrix<4, 4, T, float> res;
		res[0] = mat[0] * rotate[0][0] + mat[1] * rotate[0][1] + mat[2] * rotate[0][2];
		res[1] = mat[0] * rotate[1][0] + mat[1] * rotate[1][1] + mat[2] * rotate[1][2];
		res[2] = mat[0] * rotate[2][0] + mat[1] * rotate[2][1] + mat[2] * rotate[2][2];
		res[3] = mat[3];
		return res;
	}

	template<typename T>
	Matrix<4, 4, T, float> Scale(const Matrix<4, 4, T, float> &mat, const Vector<3, T, float> &vec)
	{
		Matrix<4, 4, T, float> res;
		res[0] = mat[0] * vec[0];
		res[1] = mat[1] * vec[1];
		res[2] = mat[2] * vec[2];
		res[3] = mat[3];
		return res;
	}

	template<typename T>
	Vector<3, T, float> GetTranslation(const Matrix<4, 4, T, float> &mat)
	{
		Vector<3, T, float> res;
		res[0] = mat[3][0];
		res[1] = mat[3][1];
		res[2] = mat[3][2];
		return res;
	}

	template<typename T>
	Vector<3, T, float> GetEulerRotation(const Matrix<4, 4, T, float> &mat)
	{
		Vector<3, T, float> mat0 = mat[0].Normalized();
		Vector<3, T, float> mat1 = mat[1].Normalized();
		Vector<3, T, float> mat2 = mat[2].Normalized();

		Vector<3, T, float> res;
		res[0] = std::atan2(mat1[2], mat2[2]);
		res[1] = std::atan2(-mat0[2], std::sqrt(mat1[2] * mat1[2] + mat2[2] * mat2[2]));
		res[2] = std::atan2(mat0[1], mat0[0]);
		return res;
	}

	template<typename T>
	Quaternion GetQuatRotation(const Matrix<4, 4, T, float> &mat)
	{
		Vector<3, T, float> mat0 = mat[0].Normalized();
		Vector<3, T, float> mat1 = mat[1].Normalized();
		Vector<3, T, float> mat2 = mat[2].Normalized();

		return Quaternion::FromEulerAngles(
			std::atan2(mat1[2], mat2[2]),
			std::atan2(-mat0[2], std::sqrt(mat1[2] * mat1[2] + mat2[2] * mat2[2])),
			std::atan2(mat0[1], mat0[0])
		);
	}

	template<typename T>
	Vector<3, T, float> GetScale(const Matrix<4, 4, T, float> &mat)
	{
		Vector<3, T, float> res;
		res[0] = mat[0].Length();
		res[1] = mat[1].Length();
		res[2] = mat[2].Length();
		return res;
	}

	template<typename T>
	void SetTranslation(Matrix<4, 4, T, float> &mat, const Vector<3, T, float> &pos)
	{
		memcpy(&mat[3], &pos, 3 * sizeof(T));
	}

	template<typename T>
	void SetRotation(Matrix<4, 4, T, float> &mat, const Vector<3, T, float> &rot)
	{
		Vector3 scale = GetScale(mat);
		Vector3 trans = GetTranslation(mat);

		mat = MakeRotation(Vector<3, T, float>::Right(), rot.x) * MakeRotation(Vector<3, T, float>::Up(), rot.y) * MakeRotation(Vector<3, T, float>::Forward(), rot.z);
		mat[0] *= scale.x;
		mat[1] *= scale.y;
		mat[2] *= scale.z;
		mat[3][0] = trans.x;
		mat[3][1] = trans.y;
		mat[3][2] = trans.z;
	}

	template<typename T>
	void SetRotation(Matrix<4, 4, T, float> &mat, const QuaternionT<T, float> &quat)
	{
		SetRotation(mat, quat.ToEulerAngles());
	}

	template<typename T>
	void SetScale(Matrix<4, 4, T, float> &mat, const Vector<3, T, float> &scale)
	{
		mat[0] = mat[0].Normalized() * scale.x;
		mat[1] = mat[1].Normalized() * scale.y;
		mat[2] = mat[2].Normalized() * scale.z;
	}
}

