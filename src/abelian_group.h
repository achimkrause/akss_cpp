#pragma once

#include <vector>

#include "types.h"

typedef std::size_t OrderExponent;

class AbelianGroup
{
 public:
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

 private:
  const std::size_t free_rank_;
  std::vector<OrderExponent> orders_;
};
