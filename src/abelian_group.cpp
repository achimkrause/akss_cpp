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

void AbelianGroup::print(std::ostream& stream, mod_t p){
  bool plus = false;
  for (dim_t i = 0; i < free_rank(); i++) {
    if (plus) std::cout << " + ";
    std::cout << "Z";
    plus = true;
  }
  for (dim_t i = 0; i < tor_rank(); i++) {
    if (plus) std::cout << " + ";
    mpz_class order =
            p_pow_z(p, orders_[i]);
    std::cout << "Z/";
    std::cout << order.get_ui();  // unsafe placeholder. FABIAAAAN :(
    plus = true;
  }
  if(!plus){
    std::cout << "0";
  }
}
