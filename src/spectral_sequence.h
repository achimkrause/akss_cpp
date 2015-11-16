#pragma once

#include <map>
#include <vector>

class TrigradedIndex
{
 public:
  TrigradedIndex(const int p, const int q, const int s);

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

  friend bool operator==(const TrigradedIndex& a, const TrigradedIndex& b);
  friend bool operator<(const TrigradedIndex& a, const TrigradedIndex& b);

 private:
  const int p_;
  const int q_;
  const int s_;
};

bool operator==(const TrigradedIndex& a, const TrigradedIndex& b);
bool operator<(const TrigradedIndex& a, const TrigradedIndex& b);

class GroupSequence {
 public:
	GroupSequence(const AbelianGroup grp);

	const AbelianGroup& getGroup(std::size_t index);
	const MatrixQ& getMatrix(std::size_t index);

 private:
	int current;
	std::map<std::size_t, AbelianGroup> groups;
	std::map<std::size_t, MatrixQ> matrices;
	//The n-th matrix represents the map between the 0-th group and the n+1-st group.
	//if number n is not set explicitly, but smaller than current (current+1 in the case of matrices),
	//its value is given by the next smaller index.



};
