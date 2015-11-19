#pragma once

#include <map>
#include <vector>
#include "abelian_group.h"
#include "morphisms.h"

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
  friend TrigradedIndex operator+(const TrigradedIndex&a, const TrigradedIndex& b);
 private:
  const int p_;
  const int q_;
  const int s_;
};

bool operator==(const TrigradedIndex& a, const TrigradedIndex& b);
bool operator<(const TrigradedIndex& a, const TrigradedIndex& b);
TrigradedIndex operator+(const TrigradedIndex&a, const TrigradedIndex& b);

class GroupSequence {
 public:
	GroupSequence(const std::size_t index_min, const AbelianGroup& grp);
	const AbelianGroup& get_group(const std::size_t index);
	const MatrixQ& get_matrix(const std::size_t index);
	void append(const std::size_t index, const AbelianGroup& grp, const MatrixQ& map);
	void done();
	std::size_t get_current();
	void inc();

 private:
	bool done_;
	std::map<std::size_t, std::tuple<AbelianGroup,MatrixQ>> entries_;
	std::size_t current_;

	//The matrix nr n represents the map between the group nr index_min and the n-th group.
	//if number n is not set explicitly, but smaller than current,
	//its value is given by the next smaller index.
	//if done is set, arbitrarily large indices are allowed,
	//and all higher things are treated as equal to the highest one that has been set explicitly.
};

class SpectralSequence {

public:
	void set_diff(TrigradedIndex pqs, std::size_t r, MatrixQ matrix);
	const AbelianGroup& get_e_ab(TrigradedIndex pqs, std::size_t a, std::size_t b);
private:
	std::map<TrigradedIndex, GroupSequence> kernels_;
	std::map<TrigradedIndex, GroupSequence> cokernels_;
	std::map<TrigradedIndex, MatrixQList> differentials_;
	//const TrigradedIndex diff_offset_; oops, depends on r. Do we want a function object for that?
	std::size_t prime_;
};
