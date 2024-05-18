#pragma once

#include "Matrix.hpp"

namespace Math
{
	template<unsigned int N, typename T>
	struct Matrix<N, N, T>
	{
		MATRIX_BODY(N, N, T, float)
		SQUARE_MAT_FUNCTIONS(N, T)

	public:
		static MatrixNxN<N, T> Identity()
		{
			Vector<N, T, float> columns[N];

			for (int i = 0; i < N; i++)
			{
				for (int j = 0; j < N; j++)
				{
					columns[i][j] = i == j;
				}
			}
		}
	};
}

#include "MatrixNxN.cpp"

