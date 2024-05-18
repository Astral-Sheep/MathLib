namespace Math
{
	// -- Constructors --

	template<unsigned int N, typename T>
	MatrixNxN<N, T>::Matrix()
	{
		for (int i = 0; i < N; i++)
		{
			values[i] = Vector<N, T, float>();
		}
	}

	template<unsigned int N, typename T>
	MatrixNxN<N, T>::Matrix(const T value)
	{
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				values[i][j] = i == j ? value : T(0);
			}
		}
	}

	template<unsigned int N, typename T>
	MatrixNxN<N, T>::Matrix(const Vector<N, T, float> columns[N])
	{
		for (int i = 0; i < N; i++)
		{
			values[i] = columns[i];
		}
	}

	template<unsigned int N, typename T>
	MatrixNxN<N, T>::Matrix(const T values[N * N])
	{
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				this->values[i][j] = values[i * N + j];
			}
		}
	}

	template<unsigned int N, typename T>
	MatrixNxN<N, T>::Matrix(const MatrixNxN<N, T> &matrix)
	{
		for (int i = 0; i < N; i++)
		{
			values[i] = matrix[i];
		}
	}

	// -- Accesses --

	template<unsigned int N, typename T>
	const Vector<N, T, float> &MatrixNxN<N, T>::operator[](const int index) const noexcept
	{
		return values[index];
	}

	template<unsigned int N, typename T>
	Vector<N, T, float> &MatrixNxN<N, T>::operator[](const int index) noexcept
	{
		return values[index];
	}

	// -- Unary arithmetic operators --

	template<unsigned int N, typename T>
	MatrixNxN<N, T> &MatrixNxN<N, T>::operator+=(const MatrixNxN<N, T> &matrix)
	{
		for (int i = 0; i < N; i++)
		{
			values[i] += matrix[i];
		}

		return *this;
	}

	template<unsigned int N, typename T>
	MatrixNxN<N, T> &MatrixNxN<N, T>::operator-=(const MatrixNxN<N, T> &matrix)
	{
		for (int i = 0; i < N; i++)
		{
			values[i] -= matrix[i];
		}

		return *this;
	}

	template<unsigned int N, typename T>
	MatrixNxN<N, T> &MatrixNxN<N, T>::operator*=(const MatrixNxN<N, T> &matrix)
	{
		MatrixNxN<N, T> temp(*this);
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

	template<unsigned int N, typename T>
	MatrixNxN<N, T> &MatrixNxN<N, T>::operator*=(const T scalar)
	{
		for (int i = 0; i < N; i++)
		{
			values[i] *= scalar;
		}

		return *this;
	}

	template<unsigned int N, typename T>
	MatrixNxN<N, T> &MatrixNxN<N, T>::operator/=(const MatrixNxN<N, T> &matrix)
	{
		return *this *= matrix.Inverted();
	}

	template<unsigned int N, typename T>
	MatrixNxN<N, T> &MatrixNxN<N, T>::operator/=(const T scalar)
	{
		for (int i = 0; i < N; i++)
		{
			values[i] /= scalar;
		}

		return *this;
	}

	// -- Unary operators --

	template<unsigned int N, typename T>
	MatrixNxN<N, T> MatrixNxN<N, T>::operator-() const
	{
		MatrixNxN<N, T> mat;

		for (int i = 0; i < N; i++)
		{
			mat[i] = -values[i];
		}

		return mat;
	}

	// -- Binary operators --

	template<unsigned int N, typename T>
	MatrixNxN<N, T> MatrixNxN<N, T>::operator+(const MatrixNxN<N, T> &matrix) const
	{
		MatrixNxN<N, T> mat;

		for (int i = 0; i < N; i++)
		{
			mat[i] = values[i] + matrix[i];
		}

		return mat;
	}

	template<unsigned int N, typename T>
	MatrixNxN<N, T> MatrixNxN<N, T>::operator-(const MatrixNxN<N, T> &matrix) const
	{
		MatrixNxN<N, T> mat;

		for (int i = 0; i < N; i++)
		{
			mat[i] = values[i] - matrix[i];
		}

		return mat;
	}

	template<unsigned int N, typename T>
	MatrixNxN<N, T> MatrixNxN<N, T>::operator*(const MatrixNxN<N, T> &matrix) const
	{
		MatrixNxN<N, T> mat;

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

	template<unsigned int N, typename T>
	template<unsigned int S>
	Matrix<S, N, T> MatrixNxN<N, T>::operator*(const Matrix<S, N, T> &matrix) const
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

	template<unsigned int N, typename T>
	MatrixNxN<N, T> MatrixNxN<N, T>::operator*(const T scalar) const
	{
		MatrixNxN<N, T> mat;

		for (int i = 0; i < N; i++)
		{
			mat[i] = values[i] * scalar;
		}

		return mat;
	}

	template<unsigned int N, typename T>
	MatrixNxN<N, T> MatrixNxN<N, T>::operator/(const MatrixNxN<N, T> &matrix) const
	{
		return *this * matrix.Inverted();
	}

	template<unsigned int N, typename T>
	MatrixNxN<N, T> MatrixNxN<N, T>::operator/(const T scalar) const
	{
		MatrixNxN<N, T> mat;

		for (int i = 0; i < N; i++)
		{
			mat[i] = values[i] / scalar;
		}

		return mat;
	}

	// -- Boolean operators --

	template<unsigned int N, typename T>
	bool MatrixNxN<N, T>::operator==(const MatrixNxN<N, T> &matrix) const
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

	template<unsigned int N, typename T>
	bool MatrixNxN<N, T>::operator!=(const MatrixNxN<N, T> &matrix) const
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

	template<unsigned int N, typename T>
	std::ostream &operator<<(std::ostream &ostream, const MatrixNxN<N, T> &matrix)
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

	template<unsigned int N, typename T>
	MatrixNxN<N, T> MatrixNxN<N, T>::GetAdjugate() const
	{
		MatrixNxN<N, T> mat;

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

	template<unsigned int N, typename T>
	MatrixNxN<N, T> MatrixNxN<N, T>::GetCofactor() const
	{
		MatrixNxN<N, T> mat;

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

	template<unsigned int N, typename T>
	T MatrixNxN<N, T>::GetDeterminant() const
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

	template<unsigned int N, typename T>
	MatrixNxN<N, T> MatrixNxN<N, T>::Inverted() const
	{
		return GetAdjugate() * (T(1) / GetDeterminant());

	}

	template<unsigned int N, typename T>
	MatrixNxN<N, T> MatrixNxN<N, T>::Transposed() const
	{
		MatrixNxN<N, T> transposed;

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

	template<unsigned int N, typename T>
	MatrixNxN<N, T> &MatrixNxN<N, T>::Invert()
	{
		*this = GetAdjugate() * (T(1) / GetDeterminant());
		return *this;
	}

	template<unsigned int N, typename T>
	MatrixNxN<N, T> &MatrixNxN<N, T>::Transpose()
	{
		const MatrixNxN<N, T> temp(*this);

		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				values[i][j] = temp[j][i];
			}
		}

		return *this;
	}
}
