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
    else if (f_rel_y(i, i) != 0)
      ++torsion_rank;
    else
      break;
  }

  AbelianGroup C(f_rel_y.height() - rank_diff - torsion_rank, torsion_rank);
  for (std::size_t i = rank_diff; i < rank_diff + torsion_rank; ++i)
    C(i - rank_diff) = static_cast<std::size_t>(p_valuation(p, f_rel_y(i, i)));

  for (MatrixQ& g_to_Y : to_Y_copy)
    to_C.emplace_back(
        g_to_Y(rank_diff, 0, g_to_Y.height() - rank_diff, g_to_Y.width()));

  for (MatrixQ& g_from_Y : from_Y_copy)
    from_C.emplace_back(g_from_Y(0, rank_diff, g_from_Y.height(),
                                 g_from_Y.width() - rank_diff));

  return C;
}


AbelianGroup compute_kernel(const std::size_t p, const MatrixQ& f,
                              const AbelianGroup& X,
                              const AbelianGroup& Y, 
                              const MatrixQRefList& to_X_ref,
                              const MatrixQRefList& from_X_ref,
                              MatrixQList& to_K, MatrixQList& from_K)
{
  MatrixQ f_rel_y(f.height(), f.width() + Y.tor_rank());
  f_rel_y(0, 0, f.height(), f.width()) = f;
  f_rel_y(0, f.width(), Y.tor_rank(), Y.tor_rank()) = Y.torsion_matrix(p);

  MatrixQ rel_x_lift(f.width() + Y.tor_rank(), X.tor_rank());
  rel_x_lift(0, 0, X.tor_rank(), X.tor_rank()) = X.torsion_matrix(p);
  //rel_x_lift(f.width(), 0, Y.tor_rank(), X.tor_rank()) = -lift of f\circ rel_x over rel_y. 
  //      Can be computed by multiplying the columns of f with the orders of X, and dividing the rows by the orders of Y. 

  for(std::size_t i=0; i<Y.tor_rank(); ++i) {
    for(std::size_t j=0; j<X.tor_rank(); ++j) {
      rel_x_lift(f.width() + i, j) = -f(i,j)*p_pow(p,X(j))/p_pow(p,Y(i));   
        //if p_pow could produce rationals of negative val, we could write p_pow(p, X(j)-Y(i)) here.
    }
  }

  //build to_X_rel_y, from_X_rel_y.
  MatrixQList to_X_rel_y;
  for (MatrixQ& g_to_X : to_X_ref) {
    MatrixQ fg = f*g_to_X;
    MatrixQ g_lift(f.width() + Y.tor_rank(), g_to_X.width());
    g_lift(0, 0, f.width(), g_to_X.width()) = g_to_X;
    for(std::size_t i=0; i<Y.tor_rank(); ++i) {
      for(std::size_t j=0; j<g_to_X.width(); ++j) {
        g_lift(f.width()+i, j) = -fg(i, j)/p_pow(p,Y(i));
      }
    }
    to_X_rel_y.emplace_back(g_lift);
  }
  
  MatrixQList from_X_rel_y;
  for (MatrixQ& g_from_X : from_X_ref) {
    MatrixQ g_rel_y(g_from_X.height(), f.width() + Y.tor_rank());
    g_rel_y(0, 0, g_from_X.height(), f.width()) = g_from_X;
    from_X_rel_y.emplace_back(g_rel_y);
    //is this copying around? should I instead do 
    //to_X_rel_y.emplace_back(g_from_X.height(), f.width() + Y.tor_rank());
    //then modify the entry in to_X_rel_y?
    //same above.
  }

  MatrixQRefList to_X_rel_y_ref = ref(to_X_rel_y);
  MatrixQRefList from_X_rel_y_ref = ref(from_X_rel_y);

  to_X_rel_y_ref.emplace_back(rel_x_lift);

  MatrixQRefList to_Y;
  MatrixQRefList from_Y;
  smith_reduce_p(p, f_rel_y, to_X_rel_y_ref, from_X_rel_y_ref, 
                              to_Y, from_Y);


  std::size_t rank_diff = 0;
  for(std::size_t i = 0; i<std::min(f_rel_y.height(), f_rel_y.width()); ++i) {
    if(f_rel_y(i,i) != 0) {
      ++rank_diff;
    }
    else {
      break;
    }
  }


  //next, restrict attention to the entries corresponding to zero columns of f_rel_y:
  // for rel_x_lift as well as the entries of to_X_rel_y, take the submatrices formed 
  // by the corresponding rows. For from_X_rel_y, take the submatrices formed by the corresponding columns.
  MatrixQ rel_K = rel_x_lift(rank_diff,0,rel_x_lift.height()-rank_diff, rel_x_lift.width());
  MatrixQList to_free_K;
  for(MatrixQ& g_to_X_rel_y : to_X_rel_y) {
    to_free_K.emplace_back(g_to_X_rel_y(rank_diff,0,
                     g_to_X_rel_y.height()-rank_diff, g_to_X_rel_y.width()));

  }
  MatrixQList from_free_K;
  for(MatrixQ& g_from_X_rel_y : from_X_rel_y) {
    to_free_K.emplace_back(g_from_X_rel_y(0,rank_diff,
                     g_from_X_rel_y.height(), g_from_X_rel_y.width()-rank_diff));

  }

  MatrixQRefList to_free_K_ref = ref(to_free_K);
  MatrixQRefList from_free_K_ref = ref(from_free_K);

  AbelianGroup free_K(rel_K.height(),0);
  //then, compute the cokernel of the new rel_x_lift with the respective to_Y, from_Y.
  return (compute_cokernel(p,rel_K,free_K,to_free_K_ref, from_free_K_ref, to_K, from_K));
}


