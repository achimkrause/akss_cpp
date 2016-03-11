#pragma once

#include "abelian_group.h"
#include "matrix.h"

struct GroupWithMorphisms {
 public:
  GroupWithMorphisms() = default;
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

GroupWithMorphisms compute_image(const std::size_t p, const MatrixQ& f,
                                 const AbelianGroup& X, const AbelianGroup& Y);

MatrixQ lift_from_free(const std::size_t p, const MatrixQ& f,
                       const MatrixQ& map, const AbelianGroup& Y);

bool morphism_equal(std::size_t p, const MatrixQ& f, const MatrixQ& g,
                    const AbelianGroup& Y);
bool morphism_zero(std::size_t p, const MatrixQ& f, const AbelianGroup& Y);
