#pragma once

#include "abelian_group.h"
#include "matrix.h"

struct GroupWithMorphisms {
 public:
  GroupWithMorphisms(const std::size_t free_rank, const std::size_t tor_rank);

  AbelianGroup group;
  MatrixQList maps_to;
  MatrixQList maps_from;
};

GroupWithMorphisms compute_cokernel(const std::size_t p, const MatrixQ& f,
                                    const AbelianGroup& Y,
                                    const MatrixQRefList& to_Y_ref,
                                    const MatrixQRefList& from_Y_ref);

GroupWithMorphisms compute_kernel(const std::size_t p, const MatrixQ& f,
                                  const AbelianGroup& X, const AbelianGroup& Y,
                                  const MatrixQRefList& to_X_ref,
                                  const MatrixQRefList& from_X_ref);

GroupWithMorphisms compute_epi_mono(const std::size_t p, const MatrixQ& f,
                                    const AbelianGroup& X,
                                    const AbelianGroup& Y);
