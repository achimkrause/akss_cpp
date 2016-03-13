#pragma once

#include <map>
#include <vector>
#include "abelian_group.h"
#include "morphisms.h"

class TrigradedIndex
{
 public:
  TrigradedIndex(const deg_t p, const deg_t q, const deg_t s);

  inline deg_t p() const
  {
    return p_;
  }
  inline deg_t q() const
  {
    return q_;
  }
  inline deg_t s() const
  {
    return s_;
  }

  friend bool operator==(const TrigradedIndex& a, const TrigradedIndex& b);
  friend bool operator<(const TrigradedIndex& a, const TrigradedIndex& b);
  friend TrigradedIndex operator+(const TrigradedIndex& a,
                                  const TrigradedIndex& b);
  friend TrigradedIndex operator-(const TrigradedIndex& a,
                                  const TrigradedIndex& b);

 private:
  const deg_t p_;
  const deg_t q_;
  const deg_t s_;
};

bool operator==(const TrigradedIndex& a, const TrigradedIndex& b);
bool operator<(const TrigradedIndex& a, const TrigradedIndex& b);
TrigradedIndex operator+(const TrigradedIndex& a, const TrigradedIndex& b);
TrigradedIndex source(const TrigradedIndex& pqs, dim_t r);
TrigradedIndex target(const TrigradedIndex& pqs, dim_t r);

class GroupSequence
{
 public:
  GroupSequence(const dim_t index_min, const AbelianGroup& grp);
  const AbelianGroup& get_group(const dim_t index) const;
  const MatrixQ& get_matrix(const dim_t index) const;
  void append(const dim_t index, const AbelianGroup& grp, const MatrixQ& map);
  void done();
  dim_t get_current() const;
  void inc();

 private:
  bool done_;
  std::map<dim_t, std::pair<AbelianGroup, MatrixQ>> entries_;
  dim_t current_;

  // The matrix nr n represents the map between the group nr index_min and the
  // n-th group.
  // if number n is not set explicitly, but smaller than current,
  // its value is given by the next smaller index.
  // if done is set, arbitrarily large indices are allowed,
  // and all higher things are treated as equal to the highest one that has been
  // set explicitly.
};

class SpectralSequence
{
 public:
  void set_diff_zero(TrigradedIndex pqs, dim_t r);
  void set_diff(TrigradedIndex pqs, dim_t r, MatrixQ matrix);
  MatrixQ get_diff_from(TrigradedIndex pqs, dim_t r) const;
  MatrixQ get_diff_to(TrigradedIndex pqs, dim_t r) const;
  GroupWithMorphisms get_e_ab(TrigradedIndex pqs, dim_t a, dim_t b) const;
  AbelianGroup get_e_2(TrigradedIndex pqs) const;
  AbelianGroup get_kernel(TrigradedIndex pqs, dim_t r) const;
  AbelianGroup get_cokernel(TrigradedIndex pqs, dim_t r) const;
  MatrixQ get_inclusion(TrigradedIndex pqs, dim_t r) const;
  MatrixQ get_projection(TrigradedIndex pqs, dim_t r) const;
  void set_e2(TrigradedIndex pqs, AbelianGroup grp);
  mod_t get_prime() const;
  std::pair<deg_t, deg_t> get_bounds(deg_t q) const;
  void set_bounds(deg_t q, deg_t min_s, deg_t max_s);

 private:
  std::map<TrigradedIndex, GroupSequence> kernels_;
  std::map<TrigradedIndex, GroupSequence> cokernels_;
  std::map<TrigradedIndex, std::map<dim_t, MatrixQ>> differentials_;
  std::map<deg_t, std::pair<deg_t, deg_t>> bounds_;
  // const TrigradedIndex diff_offset_; oops, depends on r. Do we want a
  // function object for that?
  mod_t prime_;
};
