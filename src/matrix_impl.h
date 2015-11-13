template <typename T, template <typename> class E>
MatrixExpression<T, E>::MatrixExpression(const MatrixExpression<T, E>&)
{
}

template <typename T, template <typename> class E1,
          template <typename> class E2>
bool operator==(const MatrixExpression<T, E1>& f,
                const MatrixExpression<T, E2>& g)
{
  if (f.height() != g.height()) return false;

  if (f.width() != g.width()) return false;

  for (std::size_t i = 0; i < f.height(); ++i) {
    for (std::size_t j = 0; j < f.width(); ++j) {
      if (f(i, j) != g(i, j)) return false;
    }
  }

  return true;
}

template <typename T, template <typename> class E1,
          template <typename> class E2>
bool operator!=(const MatrixExpression<T, E1>& f,
                const MatrixExpression<T, E2>& g)
{
  return !(f == g);
}

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
T IdentityMatrix<T>::operator()(const std::size_t i, const std::size_t j) const
{
  return i == j ? 1 : 0;
}

template <typename T>
Matrix<T>::Matrix(const std::size_t height, const std::size_t width)
    : height_(height), width_(width), entries_(height_ * width_)
{
}

template <typename T>
Matrix<T>::Matrix(std::initializer_list<std::initializer_list<T>> lst)
    : height_(lst.size()), width_(height_ ? lst.begin()->size() : 0)
{
  entries_.reserve(height_ * width_);

  for (const auto& row : lst) {
    for (const T& value : row) {
      entries_.push_back(value);
    }
  }
}

template <typename T>
template <template <typename> class E>
Matrix<T>::Matrix(const MatrixExpression<T, E>&& expr)
    : height_(expr.height()), width_(expr.width())
{
  entries_.reserve(height_ * width_);

  for (std::size_t i = 0; i < height_; ++i) {
    for (std::size_t j = 0; j < width_; ++j) {
      entries_.push_back(expr(i, j));
    }
  }
}

template <typename T>
T Matrix<T>::operator()(const std::size_t i, const std::size_t j) const
{
  return entries_[i * width_ + j];
}

template <typename T>
T& Matrix<T>::operator()(const std::size_t i, const std::size_t j)
{
  return entries_[i * width_ + j];
}

template <typename T>
MatrixSlice<T> Matrix<T>::operator()(const std::size_t i, const std::size_t j,
                                     const std::size_t height,
                                     const std::size_t width)
{
  return MatrixSlice<T>(*this, i, j, height, width);
}

template <typename T>
IdentityMatrix<T> Matrix<T>::identity(const std::size_t n)
{
  return IdentityMatrix<T>(n);
}

template <typename T>
Matrix<T>& Matrix<T>::row_add(const std::size_t i1, const std::size_t i2,
                              const T& lambda)
{
  for (std::size_t j = 0; j < width_; ++j) {
    entries_[i2 * width_ + j] += lambda * entries_[i1 * width_ + j];
  }
  return *this;
}

template <typename T>
Matrix<T>& Matrix<T>::row_mul(const std::size_t i, const T& lambda)
{
  for (std::size_t j = 0; j < width_; ++j) {
    entries_[i * width_ + j] *= lambda;
  }
  return *this;
}

template <typename T>
Matrix<T>& Matrix<T>::row_swap(const std::size_t i1, const std::size_t i2)
{
  using std::swap;

  for (std::size_t j = 0; j < width_; ++j) {
    swap(entries_[i1 * width_ + j], entries_[i2 * width_ + j]);
  }
  return *this;
}

template <typename T>
Matrix<T>& Matrix<T>::col_add(const std::size_t j1, const std::size_t j2,
                              const T& lambda)
{
  for (std::size_t i = 0; i < height_; ++i) {
    entries_[i * width_ + j2] += lambda * entries_[i * width_ + j1];
  }
  return *this;
}

template <typename T>
Matrix<T>& Matrix<T>::col_mul(const std::size_t j, const T& lambda)
{
  for (std::size_t i = 0; i < height_; ++i) {
    entries_[i * width_ + j] *= lambda;
  }
  return *this;
}

template <typename T>
Matrix<T>& Matrix<T>::col_swap(const std::size_t j1, const std::size_t j2)
{
  using std::swap;

  for (std::size_t i = 0; i < height_; ++i) {
    swap(entries_[i * width_ + j1], entries_[i * width_ + j2]);
  }
  return *this;
}

template <typename T>
MatrixSlice<T>::MatrixSlice(Matrix<T>& mat, const std::size_t i,
                            const std::size_t j, const std::size_t height,
                            const std::size_t width)
    : mat_(mat), i_(i), j_(j), height_(height), width_(width)
{
}

template <typename T>
MatrixSlice<T>& MatrixSlice<T>::operator=(const MatrixSlice<T>& other)
{
  if (&mat_ == &other.mat_) {
    if (i_ < other.i_ + other.height_ && i_ + height_ > other.i_ &&
        j_ < other.j_ + other.width_ && j_ + width_ > other.j_)
      throw std::logic_error("MatrixSlice::operator=: Matrix slices overlap");
  }

  return (*this =
              static_cast<const MatrixExpression<T, ::MatrixSlice>&>(other));
}

template <typename T>
template <template <typename> class E>
MatrixSlice<T>& MatrixSlice<T>::operator=(const MatrixExpression<T, E>& expr)
{
  if (height() != expr.height())
    throw std::logic_error("MatrixSlice::operator=: Dimension mismatch: " +
                           std::to_string(height()) + " != " +
                           std::to_string(expr.height()));
  if (width() != expr.width())
    throw std::logic_error("MatrixSlice::operator=: Dimension mismatch:  " +
                           std::to_string(width()) + " != " +
                           std::to_string(expr.width()));

  for (std::size_t i = 0; i < height(); ++i) {
    for (std::size_t j = 0; j < width(); ++j) {
      (*this)(i, j) = expr(i, j);
    }
  }

  return *this;
}

template <typename T>
T MatrixSlice<T>::operator()(const std::size_t i, const std::size_t j) const
{
  return mat_(i_ + i, j_ + j);
}

template <typename T>
T& MatrixSlice<T>::operator()(const std::size_t i, const std::size_t j)
{
  return mat_(i_ + i, j_ + j);
}

template <typename T>
Matrix<T> operator*(const Matrix<T>& g, const Matrix<T>& f)
{
  if (g.width() != f.height())
    throw std::logic_error("Matrix<T>::operator*: Dimension mismatch" +
                           std::to_string(g.width()) + " != " +
                           std::to_string(f.height()));

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
  for (Matrix<T>& f : to_X) {
    f.row_add(i2, i1, -lambda);
  }

  for (Matrix<T>& f : from_X) {
    f.col_add(i1, i2, lambda);
  }
}

template <typename T>
void basis_vectors_mul(MatrixRefList<T>& to_X, MatrixRefList<T>& from_X,
                       const std::size_t i, const T& lambda)
{
  for (Matrix<T>& f : to_X) {
    f.row_mul(i, 1 / lambda);
  }

  for (Matrix<T>& f : from_X) {
    f.col_mul(i, lambda);
  }
}

template <typename T>
void basis_vectors_swap(MatrixRefList<T>& to_X, MatrixRefList<T>& from_X,
                        const std::size_t i1, const std::size_t i2)
{
  for (Matrix<T>& f : to_X) {
    f.row_swap(i2, i1);
  }

  for (Matrix<T>& f : from_X) {
    f.col_swap(i1, i2);
  }
}

template <typename T>
MatrixList<T> deref(const MatrixRefList<T>& ref_list)
{
  MatrixList<T> list;
  list.reserve(ref_list.size());

  for (Matrix<T> mat : ref_list) {
    list.push_back(mat);
  }

  return list;
}

template <typename T>
MatrixRefList<T> ref(MatrixList<T>& list)
{
  MatrixRefList<T> ref_list;
  ref_list.reserve(list.size());

  for (Matrix<T>& mat : list) {
    ref_list.emplace_back(mat);
  }

  return ref_list;
}
