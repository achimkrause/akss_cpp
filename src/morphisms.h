#pragma once

#include "abelian_group.h"
#include "matrix.h"

template <unsigned int P>
struct GroupWithMorphisms<P> {
 public:
  GroupWithMorphisms(const std::size_t free_rank, const std::size_t tor_rank);

  AbelianGroup<P> group;
  MatrixQList maps_to;
  MatrixQList maps_from;
};

template <unsigned int P>
GroupWithMorphisms<P> compute_cokernel(const MatrixQ& f,
                                    const AbelianGroup<P>& Y,
                                    const MatrixQRefList& to_Y_ref,
                                    const MatrixQRefList& from_Y_ref);

template <unsigned int P>
GroupWithMorphisms<P> compute_kernel(const MatrixQ& f,
                                  const AbelianGroup<P>& X, const AbelianGroup<P>& Y,
                                  const MatrixQRefList& to_X_ref,
                                  const MatrixQRefList& from_X_ref);

template <unsigned int P>
GroupWithMorphisms<P> compute_image(const MatrixQ& f,
                                    const AbelianGroup<P>& X,
                                    const AbelianGroup<P>& Y);

template <unsigned int P>
bool morphism_equal(const MatrixQ& f, const MatrixQ& g, const AbelianGroup<P>& Y);

template <unsigned int P>
bool morphism_zero(const MatrixQ& f, const AbelianGroup<P>& Y);


#include <morphisms_impl.h>
