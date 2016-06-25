#include "spectral_sequence.h"

#include <tuple>
#include <sstream>

TrigradedIndex::TrigradedIndex(const deg_t p, const deg_t q, const deg_t s)
    : p_(p), q_(q), s_(s)
{
}

bool operator==(const TrigradedIndex& a, const TrigradedIndex& b)
{
  return (a.p_ == b.p_) && (a.q_ == b.q_) && (a.s_ == b.s_);
}

bool operator<(const TrigradedIndex& a, const TrigradedIndex& b)
{
  deg_t a_deg = a.p_ + a.q_;
  deg_t b_deg = b.p_ + b.q_;

  return std::tie(a_deg, a.p_, a.s_) < std::tie(b_deg, b.p_, b.s_);
}

std::ostream& operator<<(std::ostream& stream, const TrigradedIndex& pqs)
{
  stream << "(" << pqs.p() << ", " << pqs.q() << ", " << pqs.s() << ")";
  return stream;
}

TrigradedIndex operator+(const TrigradedIndex& a, const TrigradedIndex& b)
{
  return TrigradedIndex(a.p_ + b.p_, a.q_ + b.q_, a.s_ + b.s_);
}

TrigradedIndex operator-(const TrigradedIndex& a, const TrigradedIndex& b)
{
  return TrigradedIndex(a.p_ - b.p_, a.q_ - b.q_, a.s_ - b.s_);
}

TrigradedIndex source(const TrigradedIndex& pqs, dim_t r)
{
  const TrigradedIndex diff_offset(-static_cast<deg_t>(r),
                                   static_cast<deg_t>(r) - 1, 1);
  return pqs - diff_offset;
}

TrigradedIndex target(const TrigradedIndex& pqs, dim_t r)
{
  const TrigradedIndex diff_offset(-static_cast<deg_t>(r),
                                   static_cast<deg_t>(r) - 1, 1);
  return pqs + diff_offset;
}

GroupSequence::GroupSequence(const dim_t index_min, const AbelianGroup& grp)
    : done_(false), current_(index_min)
{
  entries_.emplace(index_min,
                   std::make_pair(grp, MatrixQ::identity(grp.rank())));
}

void GroupSequence::append(const dim_t index, const AbelianGroup& grp,
                           const MatrixQ& map)
{
  if (index <= current_) {
    throw std::logic_error("GroupSequence::append: Index is already set");
  }

  entries_.emplace(index, std::make_pair(grp, map));
  current_ = index;
}

void GroupSequence::done()
{
  done_ = true;
}

const AbelianGroup& GroupSequence::get_group(const dim_t index) const
{
  if (index < entries_.begin()->first) {
    throw std::logic_error(
        "GroupSequence::get_group: Index is less than min_index");
  }
  if (!done_ && index > current_) {
    throw std::logic_error("GroupSequence::get_group: Index is not yet set");
  }
  auto entries_it = entries_.upper_bound(index);
  --entries_it;

  return entries_it->second.first;
}

const MatrixQ& GroupSequence::get_matrix(const dim_t index) const
{
  if (index < entries_.begin()->first) {
    throw std::logic_error(
        "GroupSequence::get_matrix: Index is less than min_index");
  }
  if (!done_ && index > current_) {
    throw std::logic_error("GroupSequence::get_matrix: Index is not yet set");
  }
  auto entries_it = entries_.upper_bound(index);
  --entries_it;

  return entries_it->second.second;
}

dim_t GroupSequence::get_current() const
{
  return current_;
}

void GroupSequence::inc()
{
  ++current_;
}

void SpectralSequence::set_diff_zero(TrigradedIndex pqs, dim_t r)
{
    std::pair<deg_t, deg_t> bounds_ker = get_bounds(pqs.q());
    std::pair<deg_t, deg_t> bounds_coker = get_bounds(pqs.q() + r -1);

    if(bounds_ker.first <= pqs.s() && bounds_ker.second >= pqs.s()){
        auto kers_it = kernels_.find(pqs);
        if (kers_it == kernels_.end()) {
            std::stringstream str;
            str << "SpectralSequence::set_diff_zero: Kernel at " << pqs << " is not set.";
            throw std::logic_error(str.str());
        }
        if (kers_it->second.get_current() != r) {
            throw std::logic_error(
                    "SpectralSequence::set_diff_zero: Kernel is at wrong r.");
        }
        kers_it->second.inc();
    }
    if(bounds_coker.first <= pqs.s()+1 && bounds_coker.second >= pqs.s()+1){
        auto cokers_it = cokernels_.find(target(pqs, r));
        if (cokers_it == cokernels_.end()) {
            throw std::logic_error(
                    "SpectralSequence::set_diff_zero: Cokernel is not set.");
        }
        if (cokers_it->second.get_current() != r) {
            throw std::logic_error(
                    "SpectralSequence::set_diff_zero: Cokernel is at wrong r.");
        }
        cokers_it->second.inc();
    }
}

void SpectralSequence::set_diff(TrigradedIndex pqs, dim_t r, MatrixQ matrix)
{
  auto kers_it = kernels_.find(pqs);
  auto cokers_it = cokernels_.find(target(pqs, r));
  if (kers_it == kernels_.end()) {
    throw std::logic_error("SpectralSequence::set_diff: Kernel is not set.");
  }
  if (cokers_it == cokernels_.end()) {
    throw std::logic_error("SpectralSequence::set_diff: Cokernel is not set.");
  }
  if (kers_it->second.get_current() != r) {
    throw std::logic_error("SpectralSequence::set_diff: Kernel is at wrong r.");
  }
  if (cokers_it->second.get_current() != r) {
    throw std::logic_error(
        "SpectralSequence::set_diff: Cokernel is at wrong r.");
  }

  AbelianGroup X = kers_it->second.get_group(r);
  AbelianGroup Y = cokers_it->second.get_group(r);

  if (morphism_zero(prime_, matrix, Y)) {
    kers_it->second.inc();
    cokers_it->second.inc();
    return;
  }

  MatrixQ inc_X = kers_it->second.get_matrix(r);
  MatrixQ proj_Y = cokers_it->second.get_matrix(r);

  MatrixQList from_X, to_Y;
  from_X.emplace_back(inc_X);
  to_Y.emplace_back(proj_Y);

  GroupWithMorphisms new_kernel =
      compute_kernel(prime_, matrix, X, Y, MatrixQRefList(), ref(from_X));
  GroupWithMorphisms new_cokernel =
      compute_cokernel(prime_, matrix, Y, ref(to_Y), MatrixQRefList());

  kers_it->second.append(r + 1, new_kernel.group, new_kernel.maps_from[0]);
  cokers_it->second.append(r + 1, new_cokernel.group, new_cokernel.maps_to[0]);

  auto diffmap_it = differentials_.find(pqs);

  if (diffmap_it == differentials_.end()) {
    std::map<dim_t, MatrixQ> diffs_at_pqs;
    diffs_at_pqs.insert(std::pair<dim_t, MatrixQ>(r, matrix));
    differentials_.emplace(pqs, diffs_at_pqs);
  } else {
    diffmap_it->second.emplace(r, matrix);
  }
}

std::pair<deg_t, deg_t> SpectralSequence::get_bounds(deg_t q) const
{
  auto bounds_it = bounds_.find(q);
  if (bounds_it == bounds_.end()) {
    throw std::logic_error("SpectralSequence::get_bounds: Bounds not set.");
  }
  return bounds_it->second;
}

void SpectralSequence::set_bounds(deg_t q, deg_t min_s, deg_t max_s)
{
  auto bounds_it = bounds_.find(q);
  if (bounds_it == bounds_.end()) {
    bounds_.emplace(q, std::make_pair(min_s, max_s));
  } else {
    throw std::logic_error("SpectralSequence::set_bounds: Bounds already set.");
  }
}

MatrixQ SpectralSequence::get_diff_from(TrigradedIndex pqs, dim_t r) const
{
  TrigradedIndex pqs_target = target(pqs, r);

  if (pqs_target.p() < 0 || pqs_target.q() < 0) {
    return MatrixQ(0, 0);
  }

  std::pair<deg_t, deg_t> bounds = get_bounds(pqs.q());
  if (pqs.s() < bounds.first || pqs.s() > bounds.second) {
    return MatrixQ(0, 0);
  }

  std::pair<deg_t, deg_t> bounds_target = get_bounds(pqs_target.q());
  if (pqs_target.s() < bounds_target.first ||
      pqs_target.s() > bounds_target.second) {
    return MatrixQ(0, 0);
  }

  auto kers_it = kernels_.find(pqs);
  auto cokers_it = cokernels_.find(pqs_target);
  if (kers_it == kernels_.end()) {
    throw std::logic_error(
        "SpectralSequence::get_diff_from: Domain is not set.");
  }
  if (cokers_it == kernels_.end()) {
    throw std::logic_error(
        "SpectralSequence::get_diff_from: Codomain is not set.");
  }
  if (kers_it->second.get_current() <= r) {
    throw std::logic_error(
        "SpectralSequence::get_diff_from: Differential isn't set yet.");
  }
  if (cokers_it->second.get_current() <= r) {
    throw std::logic_error(
        "SpectralSequence::get_diff_from: Differential isn't set yet.");
  }

  auto diffmap_it = differentials_.find(pqs);
  if (diffmap_it == differentials_.end()) {
    dim_t height = cokers_it->second.get_group(r).rank();
    dim_t width = kers_it->second.get_group(r).rank();
    MatrixQ result(height, width);
    return result;
  }
  auto diff_it = diffmap_it->second.find(r);
  if (diff_it == diffmap_it->second.end()) {
    dim_t height = cokers_it->second.get_group(r).rank();
    dim_t width = kers_it->second.get_group(r).rank();
    MatrixQ result(height, width);
    return result;
  }
  return diff_it->second;
}

MatrixQ SpectralSequence::get_diff_to(TrigradedIndex pqs, std::size_t r) const
{
  return get_diff_from(source(pqs, r), r);
}

// computes the kernel of the differentials up to d_{a-1} mod the image of the
// differentials up to d_{b-1}.
// for example, get_e_ab(pqs, r, r) computes the E_r page at pqs.
GroupWithMorphisms SpectralSequence::get_e_ab(TrigradedIndex pqs, dim_t a,
                                              dim_t b) const
{
  std::pair<deg_t, deg_t> bounds = get_bounds(pqs.q());
  if (pqs.s() < bounds.first || pqs.s() > bounds.second) {
    return GroupWithMorphisms(0, 0);
  }
  auto kers_it = kernels_.find(pqs);
  auto cokers_it = cokernels_.find(pqs);
  if (kers_it == kernels_.end()) {
    throw std::logic_error("SpectralSequence::get_e_ab: Kernel is not set.");
  }
  if (cokers_it == cokernels_.end()) {
    throw std::logic_error("SpectralSequence::get_e_ab: Cokernel is not set.");
  }
  if (kers_it->second.get_current() < a) {
    throw std::logic_error("SpectralSequence::get_e_ab: Kernel is at wrong r.");
  }
  if (cokers_it->second.get_current() < b) {
    throw std::logic_error(
        "SpectralSequence::get_e_ab: Cokernel is at wrong r.");
  }

  AbelianGroup K = kers_it->second.get_group(a);
  AbelianGroup C = cokers_it->second.get_group(b);

  MatrixQ map =
      (cokers_it->second.get_matrix(b)) * (kers_it->second.get_matrix(a));

  GroupWithMorphisms I = compute_image(prime_, map, K, C);
  return I;
}

AbelianGroup SpectralSequence::get_e_2(TrigradedIndex pqs) const
{
  std::pair<deg_t, deg_t> bounds = get_bounds(pqs.q());
  if (pqs.s() < bounds.first || pqs.s() > bounds.second) {
    return AbelianGroup(0, 0);
  }
  auto kers = kernels_.find(pqs);
  if (kers == kernels_.end()) {
    throw std::logic_error("SpectralSequence::get_e_2: Group is not set.");
  }

  return kers->second.get_group(2);
}

AbelianGroup SpectralSequence::get_kernel(TrigradedIndex pqs, dim_t r) const
{
  std::pair<deg_t, deg_t> bounds = get_bounds(pqs.q());
  if (pqs.s() < bounds.first || pqs.s() > bounds.second) {
    return AbelianGroup(0, 0);
  }
  auto kers_it = kernels_.find(pqs);
  if (kers_it == kernels_.end()) {
    std::stringstream str;
    str << "SpectralSequence::get_kernel: Group at " << pqs << " is not set.";
    throw std::logic_error(str.str());
  }
  if (kers_it->second.get_current() < r) {
    throw std::logic_error(
        "SpectralSequence::get_kernel: Kernel is at wrong r.");
  }
  return kers_it->second.get_group(r);
}

AbelianGroup SpectralSequence::get_cokernel(TrigradedIndex pqs, dim_t r) const
{
  std::pair<deg_t, deg_t> bounds = get_bounds(pqs.q());
  if (pqs.s() < bounds.first || pqs.s() > bounds.second) {
    return AbelianGroup(0, 0);
  }
  auto cokers_it = cokernels_.find(pqs);
  if (cokers_it == cokernels_.end()) {
    std::stringstream str;
    str << "SpectralSequence::get_cokernel: Group at " << pqs << " is not set.";
    throw std::logic_error(str.str());
  }
  if (cokers_it->second.get_current() < r) {
    throw std::logic_error(
        "SpectralSequence::get_cokernel: Cokernel is at wrong r.");
  }
  return cokers_it->second.get_group(r);
}

MatrixQ SpectralSequence::get_inclusion(TrigradedIndex pqs, dim_t r) const
{
  std::pair<deg_t, deg_t> bounds = get_bounds(pqs.q());
  if (pqs.s() < bounds.first || pqs.s() > bounds.second) {
    return MatrixQ(0, 0);
  }
  auto kers_it = kernels_.find(pqs);
  if (kers_it == kernels_.end()) {
    throw std::logic_error(
        "SpectralSequence::get_inclusion: Group is not set.");
  }
  if (kers_it->second.get_current() < r) {
    throw std::logic_error(
        "SpectralSequence::get_inclusion: Kernel is at wrong r.");
  }
  return kers_it->second.get_matrix(r);
}

MatrixQ SpectralSequence::get_projection(TrigradedIndex pqs, dim_t r) const
{
  std::pair<deg_t, deg_t> bounds = get_bounds(pqs.q());
  if (pqs.s() < bounds.first || pqs.s() > bounds.second) {
    return MatrixQ(0, 0);
  }
  auto cokers_it = cokernels_.find(pqs);
  if (cokers_it == cokernels_.end()) {
    throw std::logic_error(
        "SpectralSequence::get_projection: Group is not set.");
  }
  if (cokers_it->second.get_current() < r) {
    throw std::logic_error(
        "SpectralSequence::get_projection: Cokernel is at wrong r.");
  }
  return cokers_it->second.get_matrix(r);
}

void SpectralSequence::set_e2(TrigradedIndex pqs, AbelianGroup grp)
{
  std::pair<deg_t, deg_t> bounds = get_bounds(pqs.q());
  if (pqs.s() < bounds.first || pqs.s() > bounds.second) {
    throw std::logic_error("SpectralSequence::set_e2: Group is already set.");
  }
  if (kernels_.find(pqs) != kernels_.end()) {
    throw std::logic_error("SpectralSequence::set_e2: Group is already set.");
  }
  GroupSequence ker2(2, grp);
  GroupSequence coker2(2, grp);
  kernels_.emplace(pqs, ker2);
  cokernels_.emplace(pqs, ker2);
}

mod_t SpectralSequence::get_prime() const
{
  return prime_;
}
