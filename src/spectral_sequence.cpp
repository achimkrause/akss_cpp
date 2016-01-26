#include "spectral_sequence.h"

#include <tuple>

TrigradedIndex::TrigradedIndex(const int p, const int q, const int s)
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

TrigradedIndex operator+(const TrigradedIndex& a, const TrigradedIndex& b)
{
  return TrigradedIndex(a.p_ + b.p_, a.q_ + b.q_, a.s_ + b.s_);
}

GroupSequence::GroupSequence(const std::size_t index_min,
                             const AbelianGroup& grp)
    : done_(false), current_(index_min)
{
  entries_.emplace(index_min,
                   std::make_pair(grp, MatrixQ::identity(grp.rank())));
}

void GroupSequence::append(const std::size_t index, const AbelianGroup& grp,
                           const MatrixQ& map)
{
  if (index <= current_) {
    throw std::logic_error("GroupSequence::append: Index is already set");
  }

  entries_.emplace(index, std::make_tuple(grp, map));
  current_ = index;
}

void GroupSequence::done()
{
  done_ = true;
}

const AbelianGroup& GroupSequence::get_group(const std::size_t index)
{
  if (index < entries_.begin()->first) {
    throw std::logic_error(
        "GroupSequence::get_group: Index is less than min_index");
  }
  if (!done_ && index > current_) {
    throw std::logic_error("GroupSequence::get_group: Index is not yet set");
  }
  std::map<std::size_t, std::tuple<AbelianGroup, MatrixQ>>::iterator pos =
      entries_.upper_bound(index);

  return (std::get<0>(pos->second));
}

const MatrixQ& GroupSequence::get_matrix(const std::size_t index)
{
  if (index < entries_.begin()->first) {
    throw std::logic_error(
        "GroupSequence::get_matrix: Index is less than min_index");
  }
  if (!done_ && index > current_) {
    throw std::logic_error("GroupSequence::get_matrix: Index is not yet set");
  }
  std::map<std::size_t, std::tuple<AbelianGroup, MatrixQ>>::iterator pos =
      entries_.upper_bound(index);

  return (std::get<1>(pos->second));
}

std::size_t GroupSequence::get_current()
{
  return current_;
}

void GroupSequence::inc()
{
  ++current_;
}
void SpectralSequence::set_diff_zero(TrigradedIndex pqs, std::size_t r){
	TrigradedIndex diff_offset(-r,r-1,1); //is that cast okay?

	std::map<TrigradedIndex,GroupSequence>::iterator
	           kers = kernels_.find(pqs);
	std::map<TrigradedIndex,GroupSequence>::iterator
		       cokers = cokernels_.find(pqs+diff_offset);
	if(kers == kernels_.end()) {
		throw std::logic_error("SpectralSequence::set_diff_zero: Kernel is not set.");
	}
	if(cokers == cokernels_.end()) {
		throw std::logic_error("SpectralSequence::set_diff_zero: Cokernel is not set.");
	}
	if(kers -> second . get_current() != r) {
		throw std::logic_error("SpectralSequence::set_diff_zero: Kernel is at wrong r.");
	}
	if(cokers -> second . get_current() != r) {
		throw std::logic_error("SpectralSequence::set_diff_zero: Cokernel is at wrong r.");
	}

	kers->second.inc();
	cokers->second.inc();
}

void SpectralSequence::set_diff(TrigradedIndex pqs, std::size_t r, MatrixQ matrix){
	TrigradedIndex diff_offset(-r,r-1,1); //is that cast okay?


	std::map<TrigradedIndex,GroupSequence>::iterator
	           kers = kernels_.find(pqs);
	std::map<TrigradedIndex,GroupSequence>::iterator
		       cokers = cokernels_.find(pqs+diff_offset);
	if(kers == kernels_.end()) {
		throw std::logic_error("SpectralSequence::set_diff: Kernel is not set.");
	}
	if(cokers == cokernels_.end()) {
		throw std::logic_error("SpectralSequence::set_diff: Cokernel is not set.");
	}
	if(kers -> second . get_current() != r) {
		throw std::logic_error("SpectralSequence::set_diff: Kernel is at wrong r.");
	}
	if(cokers -> second . get_current() != r) {
		throw std::logic_error("SpectralSequence::set_diff: Cokernel is at wrong r.");
	}

	AbelianGroup X = kers->second.get_group(r);
	AbelianGroup Y = cokers->second.get_group(r);

	if(morphism_zero(prime_,matrix, Y)){
		kers->second.inc();
		cokers->second.inc();
		return;
	}

	MatrixQ inc_X = kers->second.get_matrix(r);
	MatrixQ proj_Y = cokers->second.get_matrix(r);

	MatrixQList from_X, to_Y;
	from_X.emplace_back(inc_X);
	to_Y.emplace_back(proj_Y);

	GroupWithMorphisms new_kernel = compute_kernel(prime_, matrix, X, Y, MatrixQRefList(), ref(from_X));
	GroupWithMorphisms new_cokernel = compute_cokernel(prime_, matrix, Y, ref(to_Y),MatrixQRefList());

	kers->second.append(r+1,new_kernel.group, new_kernel.maps_from[0]);
	cokers->second.append(r+1,new_cokernel.group, new_cokernel.maps_to[0]);

	std::map<TrigradedIndex,std::map<std::size_t,MatrixQ>>::iterator
		           diffmap = differentials_.find(pqs);

	if(diffmap==differentials_.end()) {
		std::map<std::size_t,MatrixQ> diffs_at_pqs;
		diffs_at_pqs.insert(std::pair<std::size_t,MatrixQ>(r,matrix));
		differentials_.emplace(pqs,diffs_at_pqs);
	}
	else {
		diffmap->second.emplace(r,matrix);
	}
}

const MatrixQ SpectralSequence::get_diff_from(TrigradedIndex pqs, std::size_t r){
	TrigradedIndex diff_offset(-r,r-1,1);
	std::map<TrigradedIndex,GroupSequence>::iterator
	     kers = kernels_.find(pqs);
	std::map<TrigradedIndex,GroupSequence>::iterator
		 cokers = cokernels_.find(pqs+diff_offset);
	if(kers == kernels_.end()) {
		throw std::logic_error("SpectralSequence::get_diff_from: Domain is not set.");
	}
	if(cokers == kernels_.end()) {
		throw std::logic_error("SpectralSequence::get_diff_from: Codomain is not set.");
	}
	if(kers -> second . get_current() <=r) {
		throw std::logic_error("SpectralSequence::get_diff_from: Differential isn't set yet.");
	}
	if(cokers -> second . get_current() <=r) {
		throw std::logic_error("SpectralSequence::get_diff_from: Differential isn't set yet.");
	}
	std::map<TrigradedIndex,std::map<std::size_t,MatrixQ>>::iterator
	      diffmap = differentials_.find(pqs);
	if(diffmap == differentials_.end()){
		int height = cokers->second.get_group(r);
		int width = kers->second.get_group(r);
		MatrixQ result(height, width);
		return result;
	}
	std::map<std::size_t,MatrixQ>::iterator diff = diffmap->second.find(r);
	if(diff==diffmap->second.end()){
		int height = cokers->second.get_group(r);
		int width = kers->second.get_group(r);
		MatrixQ result(height, width);
		return result;
	}
	return diff->second;
}

const MatrixQ SpectralSequence::get_diff_to(TrigradedIndex pqs, std::size_t r){
	TrigradedIndex diff_offset(-r,r-1,1);
	return get_diff_from(pqs-diff_offset, r);
}

//computes the kernel of the differentials up to d_{a-1} mod the image of the differentials up to d_{b-1}.
//for example, get_e_ab(pqs, r, r) computes the E_r page at pqs.
const AbelianGroup SpectralSequence::get_e_ab(TrigradedIndex pqs, std::size_t a, std::size_t b) {
	std::map<TrigradedIndex,GroupSequence>::iterator
        kers = kernels_.find(pqs);
	std::map<TrigradedIndex,GroupSequence>::iterator
		cokers = cokernels_.find(pqs);
	if(kers == kernels_.end()) {
		throw std::logic_error("SpectralSequence::get_e_ab: Kernel is not set.");
	}
	if(cokers == cokernels_.end()) {
		throw std::logic_error("SpectralSequence::get_e_ab: Cokernel is not set.");
	}
	if(kers -> second . get_current() <a) {
		throw std::logic_error("SpectralSequence::get_e_ab: Kernel is at wrong r.");
	}
	if(cokers -> second . get_current() < b) {
		throw std::logic_error("SpectralSequence::get_e_ab: Cokernel is at wrong r.");
	}

	AbelianGroup K = kers->second.get_group(a);
	AbelianGroup C = cokers->second.get_group(b);

	MatrixQ map = (cokers->second.get_matrix(b)) * (kers->second.get_matrix(a));

	GroupWithMorphisms I = compute_image(prime_,map, K, C);
	return I.group;
}

const AbelianGroup SpectralSequence::get_e_2(TrigradedIndex pqs){
	std::map<TrigradedIndex,GroupSequence>::iterator
	        kers = kernels_.find(pqs);
	if(kers == kernels_.end()) {
		throw std::logic_error("SpectralSequence::get_e_2: Group is not set.");
	}

	return kers->second.get_group(2);
}

const AbelianGroup SpectralSequence::get_kernel(TrigradedIndex pqs, std::size_t r){
	std::map<TrigradedIndex,GroupSequence>::iterator
			kers = kernels_.find(pqs);
	if(kers == kernels_.end()) {
		throw std::logic_error("SpectralSequence::get_kernel: Group is not set.");
	}
	if(kers -> second . get_current() <r) {
		throw std::logic_error("SpectralSequence::get_kernel: Kernel is at wrong r.");
	}
	return kers->second.get_group(r);
}

const AbelianGroup SpectralSequence::get_cokernel(TrigradedIndex pqs, std::size_t r){
	std::map<TrigradedIndex,GroupSequence>::iterator
			cokers = cokernels_.find(pqs);
	if(cokers == cokernels_.end()) {
		throw std::logic_error("SpectralSequence::get_cokernel: Group is not set.");
	}
	if(cokers -> second . get_current() <r) {
		throw std::logic_error("SpectralSequence::get_cokernel: Cokernel is at wrong r.");
	}
	return cokers->second.get_group(r);
}

const MatrixQ SpectralSequence::get_inclusion(TrigradedIndex pqs, std::size_t r){
	std::map<TrigradedIndex,GroupSequence>::iterator
		kers = kernels_.find(pqs);
	if(kers == kernels_.end()) {
		throw std::logic_error("SpectralSequence::get_inclusion: Group is not set.");
	}
	if(kers -> second . get_current() <r) {
		throw std::logic_error("SpectralSequence::get_inclusion: Kernel is at wrong r.");
	}
	return kers->second.get_matrix(r);
}


const MatrixQ SpectralSequence::get_projection(TrigradedIndex pqs, std::size_t r){
	std::map<TrigradedIndex,GroupSequence>::iterator
			cokers = cokernels_.find(pqs);
	if(cokers == cokernels_.end()) {
		throw std::logic_error("SpectralSequence::get_projection: Group is not set.");
	}
	if(cokers -> second . get_current() <r) {
		throw std::logic_error("SpectralSequence::get_projection: Cokernel is at wrong r.");
	}
	return cokers->second.get_matrix(r);
}

void SpectralSequence::set_e2(TrigradedIndex pqs, AbelianGroup grp){
	std::map<TrigradedIndex,GroupSequence>::iterator
        kers = kernels_.find(pqs);
	if(kers != kernels_.end()) {
		throw std::logic_error("SpectralSequence::set_e2: Group is already set.");
	}
	kernels_.emplace(pqs,2,grp);
	cokernels_.emplace(pqs,2,grp);
}

const std::size_t SpectralSequence::get_prime(){
	return prime_;
}
