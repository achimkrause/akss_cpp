#include <exception>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

#include <gmpxx.h>

template <typename T, template <typename> class E>
class MatrixExpression
{
 public:
  MatrixExpression()
  {
  }

  inline std::size_t height() const
  {
    return static_cast<const E<T>&>(*this).height();
  }

  inline std::size_t width() const
  {
    return static_cast<const E<T>&>(*this).width();
  }

  inline T operator()(const std::size_t i, const std::size_t j) const
  {
    return static_cast<const E<T>&> (*this)(i, j);
  }

  inline T& operator()(const std::size_t i, const std::size_t j)
  {
    return static_cast<E<T>&> (*this)(i, j);
  }
};

template <typename T>
class IdentityMatrix : public MatrixExpression<T, IdentityMatrix>
{
 public:
  IdentityMatrix(const std::size_t n);

  std::size_t height() const;
  std::size_t width() const;

  T operator()(const std::size_t i, const std::size_t j) const;

 private:
  std::size_t n_;
};

template <typename T>
class MatrixSlice;

template <typename T>
class Matrix : public MatrixExpression<T, Matrix>
{
 public:
  Matrix(const std::size_t height, const std::size_t width);
  Matrix(const Matrix<T>& other);

  template <template <typename> class E>
  Matrix(const MatrixExpression<T, E>& expr);

  inline std::size_t height() const
  {
    return height_;
  }
  inline std::size_t width() const
  {
    return width_;
  }

  T operator()(const std::size_t i, const std::size_t j) const;
  T& operator()(const std::size_t i, const std::size_t j);

  MatrixSlice<T> operator()(const std::size_t i, const std::size_t j,
                            const std::size_t height, const std::size_t width);

  static IdentityMatrix<T> identity(const std::size_t n);

  //	friend Matrix<T> operator*(const Matrix<T>& g, const Matrix<T>& f);
  //  friend std::ostream& operator<< (std::ostream& stream, const Matrix<T>&
  //  f);
  friend bool operator==(const Matrix<T>& f, const Matrix<T>& g)
  {
    if (f.height_ != g.height_) return false;

    if (f.width_ != g.width_) return false;

    for (std::size_t i = 0; i < f.entries_.size(); ++i) {
      if (f.entries_[i] != g.entries_[i]) return false;
    }

    return true;
  }

 private:
  void row_add(const std::size_t i1, const std::size_t i2, const T& lambda);
  void row_mul(const std::size_t i, const T& lambda);
  void row_swap(const std::size_t i1, const std::size_t i2);

  void col_add(const std::size_t j1, const std::size_t j2, const T& lambda);
  void col_mul(const std::size_t j, const T& lambda);
  void col_swap(const std::size_t j1, const std::size_t j2);

  const std::size_t height_;
  const std::size_t width_;
  std::vector<T> entries_;
};

template <typename T>
class MatrixSlice : public MatrixExpression<T, MatrixSlice>
{
 public:
  MatrixSlice(Matrix<T>& mat, const std::size_t i, const std::size_t j,
              const std::size_t height, const std::size_t width);

  MatrixSlice<T>& operator=(const MatrixSlice<T>& other);

  template <template <typename> class E>
  MatrixSlice<T>& operator=(const MatrixExpression<T, E>& expr);

  inline std::size_t height() const
  {
    return height_;
  }

  inline std::size_t width() const
  {
    return width_;
  }

  T operator()(const std::size_t i, const std::size_t j) const;
  T& operator()(const std::size_t i, const std::size_t j);

 private:
  Matrix<T>& mat_;
  const std::size_t i_;
  const std::size_t j_;
  const std::size_t height_;
  const std::size_t width_;
};

template <typename T>
using MatrixRef = std::reference_wrapper<Matrix<T>>;
template <typename T>
using MatrixRefList = std::vector<MatrixRef<T>>;

template <typename T>
void basis_vectors_add(MatrixRefList<T>& to_X, MatrixRefList<T>& from_X,
                       const std::size_t i1, const std::size_t i2,
                       const T& lambda);
template <typename T>
void basis_vectors_mul(MatrixRefList<T>& to_X, MatrixRefList<T>& from_X,
                       const std::size_t i, const T& lambda);
template <typename T>
void basis_vectors_swap(MatrixRefList<T>& to_X, MatrixRefList<T>& from_X,
                        const std::size_t i1, const std::size_t i2);

using MatrixQ = Matrix<mpq_class>;
using MatrixQList = MatrixRefList<mpq_class>;

#include "matrix_impl.h"
