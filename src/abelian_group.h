#pragma once

#include <vector>
#include "matrix.h"
#include "types.h"

typedef std::size_t OrderExponent;

template <unsigned int P>
class AbelianGroup
{
  template <typename T>
  class TorsionMatrix : public MatrixExpression<T, TorsionMatrix>
  {
   public:
    TorsionMatrix(const AbelianGroup<P>& group);

    std::size_t height() const;
    std::size_t width() const;

    T operator()(const std::size_t i, const std::size_t j) const;

   private:
    const AbelianGroup<P>& group_;
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

  inline std::size_t p() const
  {
	return p_;
  }

  inline std::size_t free_rank() const
  {
    return free_rank_;
  }

  inline std::size_t tor_rank() const
  {
    return orders_.size();
  }

  inline std::size_t total_rank() const
  {
	  return free_rank_ + orders_.size();
  }

  TorsionMatrix<mpq_class> torsion_matrix() const;

 private:
  const std::size_t p_;
  const std::size_t free_rank_;
  std::vector<OrderExponent> orders_;
};

template <unsigned int P>
std::ostream& operator<<(std::ostream& stream, const AbelianGroup<P>& ab);


#include "abelian_group_impl.h"
