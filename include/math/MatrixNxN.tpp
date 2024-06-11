/* #include "MatrixNxN.hpp" */

namespace Math
{
#define SQUARE_MAT_TEMPLATE template<unsigned int N, typename T>
#define SQUARE_MAT_GENERIC MatrixNxN<N, T>

	// -- Constructors --

	SQUARE_MAT_TEMPLATE SQUARE_MAT_GENERIC::Matrix()
	{
		for (int i = 0; i < N; i++)
		{
			values[i] = Vector<N, T, float>();
		}
	}

	SQUARE_MAT_TEMPLATE SQUARE_MAT_GENERIC::Matrix(const T value)
	{
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				values[i][j] = i == j ? value : T(0);
			}
		}
	}

	SQUARE_MAT_TEMPLATE SQUARE_MAT_GENERIC::Matrix(const Vector<N, T, float> columns[N])
	{
		for (int i = 0; i < N; i++)
		{
			values[i] = columns[i];
		}
	}

	SQUARE_MAT_TEMPLATE SQUARE_MAT_GENERIC::Matrix(const T values[N * N])
	{
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				this->values[i][j] = values[i * N + j];
			}
		}
	}

	SQUARE_MAT_TEMPLATE SQUARE_MAT_GENERIC::Matrix(const MatrixNxN<N, T> &matrix)
	{
		for (int i = 0; i < N; i++)
		{
			values[i] = matrix[i];
		}
	}

	// -- Accesses --

	SQUARE_MAT_TEMPLATE const Vector<N, T, float> &SQUARE_MAT_GENERIC::operator[](const int index) const noexcept
	{
		return values[index];
	}

	SQUARE_MAT_TEMPLATE Vector<N, T, float> &SQUARE_MAT_GENERIC::operator[](const int index) noexcept
	{
		return values[index];
	}

	// -- Unary arithmetic operators --

	SQUARE_MAT_TEMPLATE SQUARE_MAT_GENERIC &SQUARE_MAT_GENERIC::operator+=(const SQUARE_MAT_GENERIC &matrix)
	{
		for (int i = 0; i < N; i++)
		{
			values[i] += matrix[i];
		}

		return *this;
	}

	SQUARE_MAT_TEMPLATE SQUARE_MAT_GENERIC &SQUARE_MAT_GENERIC::operator-=(const SQUARE_MAT_GENERIC &matrix)
	{
		for (int i = 0; i < N; i++)
		{
			values[i] -= matrix[i];
		}

		return *this;
	}

	SQUARE_MAT_TEMPLATE SQUARE_MAT_GENERIC &SQUARE_MAT_GENERIC::operator*=(const SQUARE_MAT_GENERIC &matrix)
	{
		SQUARE_MAT_GENERIC temp(*this);
		Vector<N, T, float> column;

		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				for (int k = 0; k < N; k++)
				{
					column[j] += temp[k][i] * matrix[j][k];
				}
			}

			values[i] = column;
		}

		return *this;
	}

	SQUARE_MAT_TEMPLATE SQUARE_MAT_GENERIC &SQUARE_MAT_GENERIC::operator*=(const T scalar)
	{
		for (int i = 0; i < N; i++)
		{
			values[i] *= scalar;
		}

		return *this;
	}

	SQUARE_MAT_TEMPLATE SQUARE_MAT_GENERIC &SQUARE_MAT_GENERIC::operator/=(const SQUARE_MAT_GENERIC &matrix)
	{
		return *this *= matrix.Inverted();
	}

	SQUARE_MAT_TEMPLATE SQUARE_MAT_GENERIC &SQUARE_MAT_GENERIC::operator/=(const T scalar)
	{
		for (int i = 0; i < N; i++)
		{
			values[i] /= scalar;
		}

		return *this;
	}

	// -- Unary operators --

	SQUARE_MAT_TEMPLATE SQUARE_MAT_GENERIC SQUARE_MAT_GENERIC::operator-() const
	{
		SQUARE_MAT_GENERIC mat;

		for (int i = 0; i < N; i++)
		{
			mat[i] = -values[i];
		}

		return mat;
	}

	// -- Binary operators --

	SQUARE_MAT_TEMPLATE SQUARE_MAT_GENERIC SQUARE_MAT_GENERIC::operator+(const SQUARE_MAT_GENERIC &matrix) const
	{
		SQUARE_MAT_GENERIC mat;

		for (int i = 0; i < N; i++)
		{
			mat[i] = values[i] + matrix[i];
		}

		return mat;
	}

	SQUARE_MAT_TEMPLATE SQUARE_MAT_GENERIC SQUARE_MAT_GENERIC::operator-(const SQUARE_MAT_GENERIC &matrix) const
	{
		SQUARE_MAT_GENERIC mat;

		for (int i = 0; i < N; i++)
		{
			mat[i] = values[i] - matrix[i];
		}

		return mat;
	}

	SQUARE_MAT_TEMPLATE SQUARE_MAT_GENERIC SQUARE_MAT_GENERIC::operator*(const SQUARE_MAT_GENERIC &matrix) const
	{
		SQUARE_MAT_GENERIC mat;

		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				for (int k = 0; k < N; k++)
				{
					mat[i][j] += values[k][j] * matrix[i][k];
				}
			}
		}

		return mat;
	}

	SQUARE_MAT_TEMPLATE template<unsigned int S>
	Matrix<S, N, T> SQUARE_MAT_GENERIC::operator*(const Matrix<S, N, T> &matrix) const
	{
		Matrix<S, N, T> mat;

		for (int i = 0; i < S; i++)
		{
			for (int j = 0; j < N; j++)
			{
				for (int k = 0; k < N; k++)
				{
					mat[i][j] += values[k][j] * matrix[i][k];
				}
			}
		}

		return mat;
	}

	SQUARE_MAT_TEMPLATE SQUARE_MAT_GENERIC SQUARE_MAT_GENERIC::operator*(const T scalar) const
	{
		SQUARE_MAT_GENERIC mat;

		for (int i = 0; i < N; i++)
		{
			mat[i] = values[i] * scalar;
		}

		return mat;
	}

	SQUARE_MAT_TEMPLATE SQUARE_MAT_GENERIC SQUARE_MAT_GENERIC::operator/(const SQUARE_MAT_GENERIC &matrix) const
	{
		return *this * matrix.Inverted();
	}

	SQUARE_MAT_TEMPLATE SQUARE_MAT_GENERIC SQUARE_MAT_GENERIC::operator/(const T scalar) const
	{
		SQUARE_MAT_GENERIC mat;

		for (int i = 0; i < N; i++)
		{
			mat[i] = values[i] / scalar;
		}

		return mat;
	}

	// -- Boolean operators --

	SQUARE_MAT_TEMPLATE bool SQUARE_MAT_GENERIC::operator==(const SQUARE_MAT_GENERIC &matrix) const
	{
		for (int i = 0; i < N; i++)
		{
			if (values[i] != matrix[i])
			{
				return false;
			}
		}

		return true;
	}

	SQUARE_MAT_TEMPLATE bool SQUARE_MAT_GENERIC::operator!=(const SQUARE_MAT_GENERIC &matrix) const
	{
		for (int i = 0; i < N; i++)
		{
			if (values[i] == matrix[i])
			{
				return false;
			}
		}

		return true;
	}

	// -- Stream operators --

	SQUARE_MAT_TEMPLATE std::ostream &operator<<(std::ostream &ostream, const SQUARE_MAT_GENERIC &matrix)
	{
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				ostream << matrix[j][i] << (j < N - 1 ? ',' : '\n');
			}
		}

		return ostream;
	}

	// -- Getters --

	SQUARE_MAT_TEMPLATE SQUARE_MAT_GENERIC SQUARE_MAT_GENERIC::GetAdjugate() const
	{
		SQUARE_MAT_GENERIC mat;

		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				Matrix<N - 1, N - 1, T> submat;
				int x = 0;
				int y = 0;

				for (int k = 0; k < N; k++)
				{
					if (k == i)
					{
						continue;
					}

					for (int l = 0; k < N; k++)
					{
						if (l == j)
						{
							continue;
						}

						submat[x][y] = values[i][j];
						y++;
					}

					x++;
				}

				mat[j][i] = submat.GetDeterminant() * (((i + j) % 2) ? -1 : 1);
			}
		}

		return mat;
	}

	SQUARE_MAT_TEMPLATE SQUARE_MAT_GENERIC SQUARE_MAT_GENERIC::GetCofactor() const
	{
		SQUARE_MAT_GENERIC mat;

		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				Matrix<N - 1, N - 1, T> submat;
				int x = 0;
				int y = 0;

				for (int k = 0; k < N; k++)
				{
					if (k == i)
					{
						continue;
					}

					for (int l = 0; k < N; k++)
					{
						if (l == j)
						{
							continue;
						}

						submat[x][y] = values[i][j];
						y++;
					}

					x++;
				}

				mat[i][j] = submat.GetDeterminant() * (((i + j) % 2) ? -1 : 1);
			}
		}

		return mat;
	}

	SQUARE_MAT_TEMPLATE T SQUARE_MAT_GENERIC::GetDeterminant() const
	{
		if (N == 0)
		{
			return 0;
		}

		if (N == 1)
		{
			return values[0][0];
		}

		if (N == 2)
		{
			return values[0][0] * values[1][1] - values[1][0] * values[0][1];
		}

		// Matrix3x3 and above

		T determinant = T(0);
		int x;
		int y;

		for (int i = 0; i < N; i++)
		{
			Matrix<N - 1, N - 1, T> subMat;
			x = 0;

			for (int j = 0; j < N; j++)
			{
				if (j == i)
				{
					continue;
				}

				y = 0;

				for (int k = 0; k < N; k++)
				{
					if (k == 0)
					{
						continue;
					}

					subMat[x][y] = values[j][k];
					y++;
				}

				x++;
			}

			determinant += values[0][i] * subMat.GetDeterminant() * ((i % 2) ? -1 : 1);
		}

		return determinant;
	}

	SQUARE_MAT_TEMPLATE SQUARE_MAT_GENERIC SQUARE_MAT_GENERIC::Inverted() const
	{
		return GetAdjugate() * (T(1) / GetDeterminant());

	}

	SQUARE_MAT_TEMPLATE SQUARE_MAT_GENERIC SQUARE_MAT_GENERIC::Transposed() const
	{
		SQUARE_MAT_GENERIC transposed;

		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				transposed[i][j] = values[j][i];
			}
		}

		return transposed;
	}

	// -- Transformations --

	SQUARE_MAT_TEMPLATE SQUARE_MAT_GENERIC &SQUARE_MAT_GENERIC::Invert()
	{
		*this = GetAdjugate() * (T(1) / GetDeterminant());
		return *this;
	}

	SQUARE_MAT_TEMPLATE SQUARE_MAT_GENERIC &SQUARE_MAT_GENERIC::Transpose()
	{
		const SQUARE_MAT_GENERIC temp(*this);

		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				values[i][j] = temp[j][i];
			}
		}

		return *this;
	}

#undef SQUARE_MAT_TEMPLATE
#undef SQUARE_MAT_GENERIC
}
