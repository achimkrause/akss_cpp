#include <gmpxx.h>

#include "abelian_group.h"
#include "p_local.h"

AbelianGroup::AbelianGroup() : free_rank_(0), orders_(0)
{
}

AbelianGroup::AbelianGroup(const std::size_t free_rank,
                           const std::size_t tor_rank)
    : free_rank_(free_rank), orders_(tor_rank)
{
}

template <>
AbelianGroup::TorsionMatrix<mpq_class>::TorsionMatrix(const AbelianGroup& group,
                                                      const std::size_t p)
    : group_(group), p_(p)
{
}

template <>
std::size_t AbelianGroup::TorsionMatrix<mpq_class>::height() const
{
  return group_.tor_rank();
}

template <>
std::size_t AbelianGroup::TorsionMatrix<mpq_class>::width() const
{
  return group_.tor_rank();
}

template <>
mpq_class AbelianGroup::TorsionMatrix<mpq_class>::operator()(
    const std::size_t i, const std::size_t j) const
{
  return i == j ? p_pow_z(p_, group_.orders_[i]) : 0;
}

AbelianGroup::TorsionMatrix<mpq_class> AbelianGroup::torsion_matrix(
    const std::size_t p) const
{
  return AbelianGroup::TorsionMatrix<mpq_class>(*this, p);
}
