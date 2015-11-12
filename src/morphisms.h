#pragma once

#include "abelian_group.h"
#include "matrix.h"

template <typename T>
struct Morphism {
  AbelianGroup domain;
  Matrix<T> map;
  AbelianGroup codomain;
};

template <typename T>
AbelianGroup compute_cokernel(const std::size_t p, const Morphism& f,
                              const MatrixRefList<T>& to_Y,
                              const MatrixRefList<T>& from_Y,
                              MatrixRefList<T>& to_C, MatrixRefList<T>& from_C);

#include "morphisms_impl.h"
