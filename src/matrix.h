#include <exception>
#include <functional>
#include <iostream>
#include <vector>

#include <gmpxx.h>

template <typename T, template <typename> class E>
class MatrixExpression
{
 public:
  std::size_t height() const
  {
    return static_cast<const E<T>&>(*this).height();
  }

  std::size_t width() const { return static_cast<const E<T>&>(*this).width(); }
  const T& operator()(const std::size_t i, const std::size_t j) const
  {
    return static_cast<const E<T>&> (*this)(i, j);
  }

  T& operator()(const std::size_t i, const std::size_t j)
  {
    return static_cast<E<T>&> (*this)(i, j);
  }
};

template <typename T>
class IdentityMatrix : MatrixExpression<T, IdentityMatrix>
{
 public:
  typedef T value_type;

  IdentityMatrix(const std::size_t n);

  std::size_t height() const;
  std::size_t width() const;

  const T& operator()(const std::size_t i, const std::size_t j) const;

 private:
  std::size_t n_;
};

template <typename T>
class Matrix : MatrixExpression<T, Matrix>
{
 public:
  typedef T value_type;

  Matrix(const std::size_t height, const std::size_t width);

  std::size_t height() const { return height_; }
  std::size_t width() const { return width_; }
  const T& operator()(const std::size_t i, const std::size_t j) const;
  T& operator()(const std::size_t i, const std::size_t j);

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
IdentityMatrix<T>::IdentityMatrix(const std::size_t n)
    : n_(n)
{
}

template <typename T>
std::size_t IdentityMatrix<T>::height() const
{
  return n_;
}

template <typename T>
std::size_t IdentityMatrix<T>::width() const
{
  return n_;
}

template <typename T>
const T& IdentityMatrix<T>::operator()(const std::size_t i,
                                       const std::size_t j) const
{
  return i == j ? 1 : 0;
}

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

template <typename T>
Matrix<T>::Matrix(const std::size_t height, const std::size_t width)
    : height_(height), width_(width), entries_(height_ * width_)
{
}

template <typename T>
const T& Matrix<T>::operator()(const std::size_t i, const std::size_t j) const
{
  return entries_[i * width_ + j];
}

template <typename T>
T& Matrix<T>::operator()(const std::size_t i, const std::size_t j)
{
  return entries_[i * width_ + j];
}

template <typename T>
IdentityMatrix<T> Matrix<T>::identity(const std::size_t n)
{
  return IdentityMatrix<T>(n);
}

template <typename T>
void Matrix<T>::row_add(const std::size_t i1, const std::size_t i2,
                        const T& lambda)
{
  for (std::size_t j = 0; j < width_; ++j) {
    entries_[i2 * width_ + j] += lambda * entries_[i1 * width_ + j];
  }
}

template <typename T>
void Matrix<T>::row_mul(const std::size_t i, const T& lambda)
{
  for (std::size_t j = 0; j < width_; ++j) {
    entries_[i * width_ + j] *= lambda;
  }
}

template <typename T>
void Matrix<T>::row_swap(const std::size_t i1, const std::size_t i2)
{
  using std::swap;

  for (std::size_t j = 0; j < width_; ++j) {
    swap(entries_[i1 * width_ + j], entries_[i2 * width_ + j]);
  }
}

template <typename T>
void Matrix<T>::col_add(const std::size_t j1, const std::size_t j2,
                        const T& lambda)
{
  for (std::size_t i = 0; i < height_; ++i) {
    entries_[i * width_ + j2] += lambda * entries_[i * width_ + j1];
  }
}

template <typename T>
void Matrix<T>::col_mul(const std::size_t j, const T& lambda)
{
  for (std::size_t i = 0; i < height_; ++i) {
    entries_[i * width_ + j] *= lambda;
  }
}

template <typename T>
void Matrix<T>::col_swap(const std::size_t j1, const std::size_t j2)
{
  using std::swap;

  for (std::size_t i = 0; i < height_; ++i) {
    swap(entries_[i * width_ + j1], entries_[i * width_ + j2]);
  }
}

template <typename T>
Matrix<T> operator*(const Matrix<T>& g, const Matrix<T>& f)
{
  if (g.width() != f.height()) throw std::logic_error("Dimension mismatch");

  Matrix<T> gf(g.height(), f.width());
  for (std::size_t i = 0; i < g.height(); ++i) {
    for (std::size_t j = 0; j < f.width(); ++j) {
      T acc;

      for (std::size_t k = 0; k < g.width(); ++k) {
        acc += g(i, k) * f(k, j);
      }

      gf(i, j) = acc;
    }
  }

  return gf;
}

template <typename T>
std::ostream& operator<<(std::ostream& stream, const Matrix<T>& f)
{
  stream << "Matrix (" << f.height() << "x" << f.width() << ")\n";

  for (std::size_t i = 0; i < f.height(); ++i) {
    for (std::size_t j = 0; j < f.width(); ++j) {
      stream << f(i, j) << " ";
    }

    stream << "\n";
  }

  return stream;
}

template <typename T>
void basis_vectors_add(MatrixRefList<T>& to_X, MatrixRefList<T>& from_X,
                       const std::size_t i1, const std::size_t i2,
                       const T& lambda)
{
  for (MatrixRef<T> f : to_X) {
    f.row_add(i2, i1, -lambda);
  }

  for (MatrixRef<T> f : from_X) {
    f.col_add(i1, i2, lambda);
  }
}

template <typename T>
void basis_vectors_mul(MatrixRefList<T>& to_X, MatrixRefList<T>& from_X,
                       const std::size_t i, const T& lambda)
{
  for (MatrixRef<T> f : to_X) {
    f.row_mul(i, 1 / lambda);
  }

  for (MatrixRef<T> f : from_X) {
    f.col_mul(i, lambda);
  }
}

template <typename T>
void basis_vectors_swap(MatrixRefList<T>& to_X, MatrixRefList<T>& from_X,
                        const std::size_t i1, const std::size_t i2)
{
  for (MatrixRef<T> f : to_X) {
    f.row_swap(i2, i1);
  }

  for (MatrixRef<T> f : from_X) {
    f.col_swap(i1, i2);
  }
}
