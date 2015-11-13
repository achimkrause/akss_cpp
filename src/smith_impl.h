#include <exception>
#include <iostream>

#include "p_local.h"

template <typename T>
void smith_reduce_p(const std::size_t p, Matrix<T>& f, MatrixRefList<T>& to_X,
                    MatrixRefList<T>& from_X, MatrixRefList<T>& to_Y,
                    MatrixRefList<T>& from_Y)
{
  from_X.emplace_back(f);
  to_Y.emplace_back(f);

  mpq_class lambda;
  for (std::size_t diagonal_block_size = 0;
       diagonal_block_size < std::min(f.height(), f.width());
       ++diagonal_block_size) {
    std::size_t i_min = 0;
    std::size_t j_min = 0;
    T min_value;
    int min_valuation = 0;

    for (std::size_t i = diagonal_block_size; i < f.height(); ++i) {
      for (std::size_t j = diagonal_block_size; j < f.width(); ++j) {
        if (f(i, j)) {
          if (!min_value || p_valuation(p, f(i, j)) < min_valuation) {
            i_min = i;
            j_min = j;
            min_value = f(i, j);
            min_valuation = p_valuation(p, min_value);
          }
        }
      }
    }

    if (!min_value) break;
    if (min_valuation < 0)
      throw std::logic_error(
          "smith_reduce_p: matrix entry has negative valuation");

    for (std::size_t i = diagonal_block_size; i < f.height(); ++i) {
      if (i == i_min) continue;
      lambda = f(i, j_min) / min_value;
      basis_vectors_add(to_Y, from_Y, i, i_min, lambda);
    }

    for (std::size_t j = diagonal_block_size; j < f.width(); ++j) {
      if (j == j_min) continue;
      lambda = -f(i_min, j) / min_value;
      basis_vectors_add(to_X, from_X, j_min, j, lambda);
    }

    basis_vectors_swap(to_Y, from_Y, i_min, diagonal_block_size);
    basis_vectors_swap(to_X, from_X, j_min, diagonal_block_size);

    lambda = p_pow(p, static_cast<std::size_t>(min_valuation)) / min_value;
    basis_vectors_mul(to_X, from_X, diagonal_block_size, lambda);
  }

  from_X.pop_back();
  to_Y.pop_back();
}
