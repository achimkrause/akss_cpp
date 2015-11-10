#include <exception>
#include <functional>
#include <iostream>
#include <vector>

#include <gmpxx.h>


template <typename T> 
class Matrix {
public:
	Matrix(const std::size_t height, const std::size_t width);

	const T& operator()(const std::size_t i, const std::size_t j) const;
	T& operator()(const std::size_t i, const std::size_t j);

	static Matrix<T> identity(const std::size_t n);

	//	friend Matrix<T> operator*(const Matrix<T>& g, const Matrix<T>& f);
	//  friend std::ostream& operator<< (std::ostream& stream, const Matrix<T>& f);
	friend bool operator==(const Matrix<T>& f, const Matrix<T>& g) 
	{
		if (f.height != g.height)
			return false;

		if (f.width != g.width)
			return false;

		for (std::size_t i = 0; i < f.entries_.size(); ++i)
		{
			if (f.entries_[i] != g.entries_[i])
				return false;
		}

		return true;
	}

	const std::size_t height;
	const std::size_t width;

private:
	void row_add(const std::size_t i1, const std::size_t i2, const T& lambda);
	void row_mul(const std::size_t i, const T& lambda);
	void row_swap(const std::size_t i1, const std::size_t i2);

 	void col_add(const std::size_t j1, const std::size_t j2, const T& lambda);
	void col_mul(const std::size_t j,  const T& lambda);
	void col_swap(const std::size_t j1, const std::size_t j2);

 	std::vector<T> entries_;
};

template <typename T>
using MatrixRef= std::reference_wrapper<Matrix<T>>;
template <typename T>
using MatrixRefList = std::vector<MatrixRef<T>>; 

template <typename T>
void basis_vectors_add(MatrixRefList<T>& to_X, MatrixRefList<T>& from_X, const std::size_t i1, const std::size_t i2, const T& lambda);
template <typename T>
void basis_vectors_mul(MatrixRefList<T>& to_X, MatrixRefList<T>& from_X, const std::size_t i, const T& lambda);
template <typename T>
void basis_vectors_swap(MatrixRefList<T>& to_X, MatrixRefList<T>& from_X, const std::size_t i1, const std::size_t i2);


using MatrixQ = Matrix<mpq_class>;
using MatrixQList = MatrixRefList<mpq_class>;


template <typename T>
Matrix<T>::Matrix(const std::size_t height, const std::size_t width) :
	height(height), width(width), entries_(height * width) 
{}

template <typename T>
const T& Matrix<T>::operator()(const std::size_t i, const std::size_t j) const
{
	return entries_[i * width + j];
}

template <typename T>
T& Matrix<T>::operator()(const std::size_t i, const std::size_t j)
{
	return entries_[i * width + j];
}

template <typename T>
Matrix<T> Matrix<T>::identity(const std::size_t n)
{
	Matrix<T> id(n, n);

	for (std::size_t i = 0; i < n; ++i)
	{
		id(i, i) = 1;
	}

	return id;
}

template <typename T>
void Matrix<T>::row_add(const std::size_t i1, const std::size_t i2, const T& lambda)
{
	for (std::size_t j = 0; j < width; ++j)
	{
		entries_[i2*width+j] += lambda*entries_[i1*width+j];
	}
}

template <typename T>
void Matrix<T>::row_mul(const std::size_t i, const T& lambda)
{
	for (std::size_t j = 0; j < width; ++j)
	{
		entries_[i*width+j] *= lambda;
	}
}

template <typename T>
void Matrix<T>::row_swap(const std::size_t i1, const std::size_t i2)
{
	using std::swap;

	for (std::size_t j = 0; j < width; ++j)
	{
		swap(entries_[i1*width+j], entries_[i2*width+j]);
	}
}

template <typename T>
void Matrix<T>::col_add(const std::size_t j1, const std::size_t j2, const T& lambda)
{
	for (std::size_t i = 0; i < width; ++i)
	{
		entries_[i*width+j2] += lambda*entries_[i*width+j1];
	}
}

template <typename T>
void Matrix<T>::col_mul(const std::size_t j, const T& lambda)
{
	for (std::size_t i = 0; i < width; ++i)
	{
		entries_[i*width+j] *= lambda;
	}
}

template <typename T>
void Matrix<T>::col_swap(const std::size_t j1, const std::size_t j2)
{
	using std::swap;

	for (std::size_t i = 0; i < width; ++i)
	{
		swap(entries_[i*width+j1], entries_[i*width+j2]);
	}
}

template <typename T>
Matrix<T> operator*(const Matrix<T>& g, const Matrix<T>& f) 
{
	if (g.width != f.height)
		throw std::logic_error("Dimension mismatch");

	Matrix<T> gf(g.height, f.width);
	for (std::size_t i = 0; i < g.height; ++i)
	{
		for (std::size_t j = 0; j < f.width; ++j)
		{
			T acc;

			for (std::size_t k = 0; k < g.width; ++k)
			{
				acc += g(i, k) * f(k, j);
			}

			gf(i, j) = acc;
		}
	}

	return gf;
}

template <typename T>
std::ostream& operator<< (std::ostream& stream, const Matrix<T>& f)
{
	stream << "Matrix (" << f.height << "x" << f.width << ")\n";

	for (std::size_t i = 0; i < f.height; ++i)
	{
		for (std::size_t j = 0; j < f.width; ++j)
		{
			stream << f(i, j) << " ";
		}

		stream << "\n";
	}

	return stream;
}


template <typename T>
void basis_vectors_add(MatrixRefList<T>& to_X, MatrixRefList<T>& from_X, const std::size_t i1, const std::size_t i2, const T& lambda)
{
	for (MatrixRef<T> f : to_X)
	{
		f.row_add(i2, i1, -lambda);
	}

	for (MatrixRef<T> f : from_X)
	{
		f.col_add(i1, i2, lambda);
	}
}

template <typename T>
void basis_vectors_mul(MatrixRefList<T>& to_X, MatrixRefList<T>& from_X, const std::size_t i, const T& lambda)
{
	for (MatrixRef<T> f : to_X)
	{
		f.row_mul(i, 1 / lambda);
	}

	for (MatrixRef<T> f : from_X)
	{
		f.col_mul(i, lambda);
	}
}

template <typename T>
void basis_vectors_swap(MatrixRefList<T>& to_X, MatrixRefList<T>& from_X, const std::size_t i1, const std::size_t i2)
{
	for (MatrixRef<T> f : to_X)
	{
		f.row_swap(i2, i1);
	}

	for (MatrixRef<T> f : from_X)
	{
		f.col_swap(i1, i2);
	}
}
