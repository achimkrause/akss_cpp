#pragma once

#include <map>
#include <vector>

class BigradedIndex
{
 public:
  BigradedIndex(const std::size_t p, const std::size_t q, const std::size_t s);

  inline std::size_t p() const
  {
    return p_;
  }
  inline std::size_t q() const
  {
    return q_;
  }
  inline std::size_t s() const
  {
    return s_;
  }

  friend bool operator==(const BigradedIndex& a, const BigradedIndex& b);
  friend bool operator<(const BigradedIndex& a, const BigradedIndex& b);

 private:
  const std::size_t p_;
  const std::size_t q_;
  const std::size_t s_;
};

bool operator==(const BigradedIndex& a, const BigradedIndex& b);
bool operator<(const BigradedIndex& a, const BigradedIndex& b);
