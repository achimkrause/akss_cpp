#pragma once

#include "abelian_group.h"
#include "matrix.h"

AbelianGroup compute_cokernel(const std::size_t p, const MatrixQ& f,
                              const AbelianGroup& Y,
                              const MatrixQRefList& to_Y_ref,
                              const MatrixQRefList& from_Y_ref,
                              MatrixQList& to_C, MatrixQList& from_C);

AbelianGroup compute_kernel(const std::size_t p, const MatrixQ& f,
                              const AbelianGroup& X,
                              const AbelianGroup& Y, 
                              const MatrixQRefList& to_X_ref,
                              const MatrixQRefList& from_X_ref,
                              MatrixQList& to_K, MatrixQList& from_K);
