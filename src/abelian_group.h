#pragma once

#include <vector>

#include "matrix.h"
#include "types.h"

typedef std::size_t OrderExponent;

class AbelianGroup
{
  template <typename T>
  class TorsionMatrix : public MatrixExpression<T, TorsionMatrix>
  {
   public:
    TorsionMatrix(const AbelianGroup& group, const std::size_t p);

    std::size_t height() const;
    std::size_t width() const;

    T operator()(const std::size_t i, const std::size_t j) const;

   private:
    const AbelianGroup& group_;
    const std::size_t p_;
  };

 public:
  AbelianGroup();
  AbelianGroup(const std::size_t free_rank, const std::size_t tor_rank);

  inline OrderExponent operator()(const std::size_t i) const
  {
    return orders_[i];
  }

  inline OrderExponent& operator()(const std::size_t i)
  {
    return orders_[i];
  }

  inline std::size_t free_rank() const
  {
    return free_rank_;
  }

  inline std::size_t tor_rank() const
  {
    return orders_.size();
  }

  inline std::size_t rank() const
  {
	  return free_rank() + tor_rank();
  }

  TorsionMatrix<mpq_class> torsion_matrix(const std::size_t p) const;

 private:
  const std::size_t free_rank_;
  std::vector<OrderExponent> orders_;
};
