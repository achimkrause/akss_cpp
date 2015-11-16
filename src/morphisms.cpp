#include "morphisms.h"

#include <iostream>

#include "abelian_group.h"
#include "matrix.h"
#include "p_local.h"
#include "smith.h"

GroupWithMorphisms::GroupWithMorphisms(const std::size_t free_rank, const
    std::size_t tor_rank)
    : group(free_rank, tor_rank)
{
}

GroupWithMorphisms compute_cokernel(const std::size_t p, const MatrixQ& f,
                                    const AbelianGroup& Y,
                                    const MatrixQRefList& to_Y_ref,
                                    const MatrixQRefList& from_Y_ref)
{
  MatrixQ f_rel_Y(f.height(), f.width() + Y.tor_rank());
  f_rel_Y(0, 0, f.height(), f.width()) = f;
  f_rel_Y(0, f.width(), Y.tor_rank(), Y.tor_rank()) = Y.torsion_matrix(p);

  MatrixQList to_Y_copy = deref(to_Y_ref);
  MatrixQList from_Y_copy = deref(from_Y_ref);

  MatrixQRefList to_X;
  MatrixQRefList from_X;
  MatrixQRefList to_Y_copy_ref = ref(to_Y_copy);
  MatrixQRefList from_Y_copy_ref = ref(from_Y_copy);

  smith_reduce_p(p, f_rel_Y, to_X, from_X, to_Y_copy_ref, from_Y_copy_ref);

  std::size_t rank_diff = 0;
  std::size_t torsion_rank = 0;

  for (std::size_t i = 0; i < std::min(f_rel_Y.height(), f_rel_Y.width());
       ++i) {
    if (f_rel_Y(i, i) == 1)
      ++rank_diff;
    else if (f_rel_Y(i, i) != 0)
      ++torsion_rank;
    else
      break;
  }

  GroupWithMorphisms C(f_rel_Y.height() - rank_diff - torsion_rank, torsion_rank);
  for (std::size_t i = rank_diff; i < rank_diff + torsion_rank; ++i)
    C.group(i - rank_diff) =
        static_cast<std::size_t>(p_val_q(p, f_rel_Y(i, i)));

  for (MatrixQ& g_to_Y : to_Y_copy)
    C.maps_to.emplace_back(
        g_to_Y(rank_diff, 0, g_to_Y.height() - rank_diff, g_to_Y.width()));

  for (MatrixQ& g_from_Y : from_Y_copy)
    C.maps_from.emplace_back(g_from_Y(0, rank_diff, g_from_Y.height(),
                                      g_from_Y.width() - rank_diff));

  return C;
}

GroupWithMorphisms compute_kernel(const std::size_t p, const MatrixQ& f,
                                  const AbelianGroup& X, const AbelianGroup& Y,
                                  const MatrixQRefList& to_X_ref,
                                  const MatrixQRefList& from_X_ref)
{
  MatrixQ f_rel_Y(f.height(), f.width() + Y.tor_rank());
  f_rel_Y(0, 0, f.height(), f.width()) = f;
  f_rel_Y(0, f.width(), Y.tor_rank(), Y.tor_rank()) = Y.torsion_matrix(p);

  MatrixQ rel_x_lift(f.width() + Y.tor_rank(), X.tor_rank());
  rel_x_lift(0, 0, X.tor_rank(), X.tor_rank()) = X.torsion_matrix(p);
  // rel_x_lift(f.width(), 0, Y.tor_rank(), X.tor_rank()) = -lift of f\circ
  // rel_x over rel_Y.
  //      Can be computed by multiplying the columns of f with the orders of X,
  //      and dividing the rows by the orders of Y.

  for (std::size_t i = 0; i < Y.tor_rank(); ++i) {
    for (std::size_t j = 0; j < X.tor_rank(); ++j) {
      long order_diff = static_cast<long>(X(j) - Y(i));
      rel_x_lift(f.width() + i, j) = -f(i, j) * p_pow_q(p, order_diff);
    }
  }

  // build to_X_rel_Y, from_X_rel_Y.
  MatrixQList to_X_rel_Y;
  for (MatrixQ& g_to_X : to_X_ref) {
    to_X_rel_Y.emplace_back(f.width() + Y.tor_rank(), g_to_X.width());
    to_X_rel_Y.back()(0, 0, f.width(), g_to_X.width()) = g_to_X;

    MatrixQ fg = f * g_to_X;
    for (std::size_t i = 0; i < Y.tor_rank(); ++i) {
      for (std::size_t j = 0; j < g_to_X.width(); ++j) {
        to_X_rel_Y.back()(f.width() + i, j) = -fg(i, j) / p_pow_z(p, Y(i));
      }
    }
  }

  MatrixQList from_X_rel_Y;
  for (MatrixQ& g_from_X : from_X_ref) {
    from_X_rel_Y.emplace_back(g_from_X.height(), f.width() + Y.tor_rank());
    from_X_rel_Y.back()(0, 0, g_from_X.height(), f.width()) = g_from_X;
  }

  MatrixQRefList to_X_rel_Y_ref = ref(to_X_rel_Y);
  MatrixQRefList from_X_rel_Y_ref = ref(from_X_rel_Y);

  to_X_rel_Y_ref.emplace_back(rel_x_lift);

  MatrixQRefList to_Y;
  MatrixQRefList from_Y;
  smith_reduce_p(p, f_rel_Y, to_X_rel_Y_ref, from_X_rel_Y_ref, to_Y, from_Y);

  std::size_t rank_diff;
  for (rank_diff = 0; rank_diff < std::min(f_rel_Y.height(), f_rel_Y.width());
       ++rank_diff) {
    if (f_rel_Y(rank_diff, rank_diff) == 0) break;
  }

  // next, restrict attention to the entries corresponding to zero columns of
  // f_rel_Y:
  // for rel_x_lift as well as the entries of to_X_rel_Y, take the submatrices
  // formed
  // by the corresponding rows. For from_X_rel_Y, take the submatrices formed
  // by
  // the corresponding columns.
  MatrixQ rel_K = rel_x_lift(rank_diff, 0, rel_x_lift.height() - rank_diff,
                             rel_x_lift.width());
  MatrixQList to_free_K;
  for (MatrixQ& g_to_X_rel_Y : to_X_rel_Y) {
    to_free_K.emplace_back(g_to_X_rel_Y(
        rank_diff, 0, g_to_X_rel_Y.height() - rank_diff, g_to_X_rel_Y.width()));
  }
  MatrixQList from_free_K;

  for(MatrixQ& g_from_X_rel_Y : from_X_rel_Y) {
    from_free_K.emplace_back(g_from_X_rel_Y(0,rank_diff,
                     g_from_X_rel_Y.height(), g_from_X_rel_Y.width()-rank_diff));

  }

  MatrixQRefList to_free_K_ref = ref(to_free_K);
  MatrixQRefList from_free_K_ref = ref(from_free_K);

  AbelianGroup free_K(rel_K.height(), 0);
  // then, compute the cokernel of the new rel_x_lift with the respective
  // to_Y, from_Y.
  return compute_cokernel(p, rel_K, free_K, to_free_K_ref, from_free_K_ref);
}



int morphism_equal(std::size_t p, const MatrixQ& f, const MatrixQ& g, const AbelianGroup& Y) {

	for(std::size_t i=0; i<f.height(); i++) {
		for(std::size_t j=0; j<f.width(); j++) {
			if(f(i,j)!=g(i,j)) {
				if(i<Y.tor_rank()) {
					if(p_val_q(p, f(i,j)-g(i,j)) < Y(i)) {
						return 0;
					}
				}
				else {
					return 0;
				}
			}
		}
	}
	return 1;
}

int morphism_zero(std::size_t p, const MatrixQ& f, const AbelianGroup& Y) {
	MatrixQ g(f.height(), f.width());

	return morphism_equal(p, f, g, Y);
}

GroupWithMorphisms compute_epi_mono(const std::size_t p, const MatrixQ& f,
                                    const AbelianGroup& X,
                                    const AbelianGroup& Y)
{
  MatrixQRefList to_X;
  MatrixQList from_X = {MatrixQ::identity(f.width())};
  GroupWithMorphisms K = compute_kernel(p, f, X, Y, to_X, ref(from_X));

  MatrixQList to_X_2 = {MatrixQ::identity(f.width())};
  MatrixQList from_X_2 = {f};
  GroupWithMorphisms img =
      compute_cokernel(p, f, Y, ref(to_X_2), ref(from_X_2));

  return img;
}
