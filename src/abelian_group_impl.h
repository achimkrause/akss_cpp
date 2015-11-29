#include <gmpxx.h>

#include "p_local.h"

template <unsigned int P>
AbelianGroup<P>::AbelianGroup() : free_rank_(0), orders_(0)
{
}

template <unsigned int P>
AbelianGroup<P>::AbelianGroup(const std::size_t free_rank,
                           const std::size_t tor_rank)
    : free_rank_(free_rank), orders_(tor_rank)
{
}

template <unsigned int P>
template <>
AbelianGroup<P>::TorsionMatrix<mpq_class>::TorsionMatrix(const AbelianGroup<P>& group)
    : group_(group)
{
}

template <unsigned int P>
template <>
std::size_t AbelianGroup<P>::TorsionMatrix<mpq_class>::height() const
{
  return group_.tor_rank();
}

template <unsigned int P>
template <>
std::size_t AbelianGroup<P>::TorsionMatrix<mpq_class>::width() const
{
  return group_.tor_rank();
}

template <unsigned int P>
template <>
mpq_class AbelianGroup<P>::TorsionMatrix<mpq_class>::operator()(
    const std::size_t i, const std::size_t j) const
{
  return i == j ? p_pow_z(P, group_.orders_[i]) : 0;
}

template <unsigned int P>
AbelianGroup<P>::TorsionMatrix<mpq_class> AbelianGroup<P>::torsion_matrix() const
{
  return AbelianGroup<P>::TorsionMatrix<mpq_class>(*this);
}

template <unsigned int P>
std::ostream& operator<<(std::ostream& stream, const AbelianGroup<P>& ab)
{

	std::size_t run=1;
	OrderExponent ord = ab(1);
	for(int i=1; i<ab.tor_rank(); i++){
		if(ord == ab(i)){
			run++;
		}
		else {
			stream << "Z/";
			stream << p_pow_z(P, ab(i-1));
			stream << "^" << run << " + ";
			run=1;
			ord=ab(i);
		}
	}
	stream << "Z/";
	stream << "^" << run << " + ";

	stream << "Z^" << ab.free_rank();

	return stream;
}

