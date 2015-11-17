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


GroupSequence::GroupSequence(const std::size_t index_min, const AbelianGroup& grp)
	: done_(false)
{
	entries_.emplace(index_min,std::make_tuple(grp, MatrixQ::identity(grp.rank())));
}

void GroupSequence::append(const std::size_t index, const AbelianGroup& grp, const MatrixQ& map){
	if(index <= entries_.rbegin()->first) {
		throw std::logic_error("GroupSequence::append: Index is already set");
	}

	entries_.emplace(index, std::make_tuple(grp,map));
}

void GroupSequence::done() {
	done_ = true;
}

const AbelianGroup& GroupSequence::get_group(const std::size_t index){
	if(index > entries_.rbegin()->first) {
		throw std::logic_error("GroupSequence::get_group: Index is not yet set");
	}
	if(index < entries_.begin()->first) {
		throw std::logic_error("GroupSequence::get_group: Index is less than min_index");
	}
	std::map<std::size_t, std::tuple<AbelianGroup,MatrixQ>>::iterator
	    pos = entries_.upper_bound(index);

	return(std::get<0>(pos->second));
}

const MatrixQ& GroupSequence::get_matrix(const std::size_t index){
	if(index > entries_.rbegin()->first) {
		throw std::logic_error("GroupSequence::get_group: Index is not yet set");
	}
	if(index < entries_.begin()->first) {
		throw std::logic_error("GroupSequence::get_group: Index is less than min_index");
	}
	std::map<std::size_t, std::tuple<AbelianGroup,MatrixQ>>::iterator
	    pos = entries_.upper_bound(index);

	return(std::get<1>(pos->second));
}
