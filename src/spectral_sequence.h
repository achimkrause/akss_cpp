#pragma once

#include <map>
#include <vector>

class BigradedIndex
{
 public:
  BigradedIndex(const int p, const int q, const int s);

  inline int p() const
  {
    return p_;
  }
  inline int q() const
  {
    return q_;
  }
  inline int s() const
  {
    return s_;
  }

  friend bool operator==(const BigradedIndex& a, const BigradedIndex& b);
  friend bool operator<(const BigradedIndex& a, const BigradedIndex& b);

 private:
  const int p_;
  const int q_;
  const int s_;
};

bool operator==(const BigradedIndex& a, const BigradedIndex& b);
bool operator<(const BigradedIndex& a, const BigradedIndex& b);
