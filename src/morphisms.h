#pragma once

#include "abelian_group.h"
#include "matrix.h"

AbelianGroup compute_cokernel(const std::size_t p, const MatrixQ& f,
                              const AbelianGroup& Y,
                              const MatrixQRefList& to_Y_ref,
                              const MatrixQRefList& from_Y_ref,
                              MatrixQList& to_C, MatrixQList& from_C);

