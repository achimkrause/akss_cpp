#pragma once

#include <vector>

#include "matrix.h"
#include "types.h"

class AbelianGroup
{
  template <typename T>
  class TorsionMatrix : public MatrixExpression<T, TorsionMatrix>
  {
   public:
    TorsionMatrix(const AbelianGroup& group, const mod_t p);

    dim_t height() const;
    dim_t width() const;

    T operator()(const dim_t i, const dim_t j) const;

   private:
    const AbelianGroup& group_;
    const mod_t p_;
  };

 public:
  AbelianGroup() = default;
  AbelianGroup(const dim_t free_rank, const dim_t tor_rank);

  inline dim_t operator()(const dim_t i) const
  {
    return orders_[i];
  }

  inline dim_t& operator()(const dim_t i)
  {
    return orders_[i];
  }

  inline dim_t free_rank() const
  {
    return free_rank_;
  }

  inline dim_t tor_rank() const
  {
    return orders_.size();
  }

  inline dim_t rank() const
  {
    return free_rank() + tor_rank();
  }

  TorsionMatrix<mpq_class> torsion_matrix(const mod_t p) const;
  void print(std::ostream& stream, mod_t p);
 private:
  dim_t free_rank_;
  std::vector<dim_t> orders_;
};
