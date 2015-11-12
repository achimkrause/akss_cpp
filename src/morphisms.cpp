#include "morphisms.h"

#include <iostream>

#include "abelian_group.h"
#include "matrix.h"
#include "p_local.h"
#include "smith.h"

AbelianGroup compute_cokernel(const std::size_t p, const MatrixQ& f,
                              const AbelianGroup& Y,
                              const MatrixQRefList& to_Y_ref,
                              const MatrixQRefList& from_Y_ref,
                              MatrixQList& to_C, MatrixQList& from_C)
{
  MatrixQ f_rel_y(f.height(), f.width() + Y.tor_rank());
  f_rel_y(0, 0, f.height(), f.width()) = f;
  f_rel_y(0, f.width(), Y.tor_rank(), Y.tor_rank()) = Y.torsion_matrix(p);

  MatrixQList to_Y_copy = deref(to_Y_ref);
  MatrixQList from_Y_copy = deref(from_Y_ref);

  MatrixQRefList to_X;
  MatrixQRefList from_X;
  MatrixQRefList to_Y_copy_ref = ref(to_Y_copy);
  MatrixQRefList from_Y_copy_ref = ref(from_Y_copy);

  smith_reduce_p(p, f_rel_y, to_X, from_X, to_Y_copy_ref, from_Y_copy_ref);

  std::size_t rank_diff = 0;
  std::size_t torsion_rank = 0;

  for (std::size_t i = 0; i < std::min(f_rel_y.height(), f_rel_y.width());
       ++i) {
    if (f_rel_y(i, i) == 1)
      ++rank_diff;
    else if (f_rel_y(i, i))
      ++torsion_rank;
    else
      break;
  }

  AbelianGroup C(f_rel_y.height() - rank_diff - torsion_rank, torsion_rank);
  for (std::size_t i = rank_diff; i < rank_diff + torsion_rank; ++i)
    C(i - rank_diff) = p_valuation(p, f_rel_y(i, i));

  for (MatrixQ& g_to_Y : to_Y_copy)
    to_C.emplace_back(
        g_to_Y(rank_diff, 0, g_to_Y.height() - rank_diff, g_to_Y.width()));

  for (MatrixQ& g_from_Y : from_Y_copy)
    from_C.emplace_back(g_from_Y(0, rank_diff, g_from_Y.height(),
                                 g_from_Y.width() - rank_diff));

  return C;
}
