#include <gmpxx.h>

#include "abelian_group.h"
#include "p_local.h"

AbelianGroup::AbelianGroup(const dim_t free_rank,
                           const dim_t tor_rank)
    : free_rank_(free_rank), orders_(tor_rank)
{
}

template <>
AbelianGroup::TorsionMatrix<mpq_class>::TorsionMatrix(const AbelianGroup& group,
                                                      const mod_t p)
    : group_(group), p_(p)
{
}

template <>
dim_t AbelianGroup::TorsionMatrix<mpq_class>::height() const
{
  return group_.tor_rank();
}

template <>
dim_t AbelianGroup::TorsionMatrix<mpq_class>::width() const
{
  return group_.tor_rank();
}

template <>
mpq_class AbelianGroup::TorsionMatrix<mpq_class>::operator()(
    const dim_t i, const dim_t j) const
{
  return i == j ? p_pow_z(p_, group_.orders_[i]) : 0;
}

AbelianGroup::TorsionMatrix<mpq_class> AbelianGroup::torsion_matrix(
    const mod_t p) const
{
  return AbelianGroup::TorsionMatrix<mpq_class>(*this, p);
}
