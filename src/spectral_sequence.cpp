#include "spectral_sequence.h"

#include <tuple>

TrigradedIndex::TrigradedIndex(const int p, const int q,
                             const int s)
    : p_(p), q_(q), s_(s)
{
}

bool operator==(const TrigradedIndex& a, const TrigradedIndex& b)
{
  return (a.p_ == b.p_) && (a.q_ == b.q_) && (a.s_ == b.s_);
}

bool operator<(const TrigradedIndex& a, const TrigradedIndex& b)
{
  int a_deg = a.p_ + a.q_;
  int b_deg = b.p_ + b.q_;

  return std::tie(a_deg, a.p_, a.s_) < std::tie(b_deg, b.p_, b.s_);
}
