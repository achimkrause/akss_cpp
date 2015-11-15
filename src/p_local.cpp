#include <exception>

#include "p_local.h"

std::size_t p_val_z(const std::size_t p, const mpz_class& x)
{
  if (x == 0) throw std::logic_error("p_valuation: x=0");

  std::size_t val = 0;
  mpz_class remainder = x;

  while (mpz_divisible_ui_p(remainder.get_mpz_t(), p)) {
    val++;
    mpz_divexact_ui(remainder.get_mpz_t(), remainder.get_mpz_t(), p);
  }

  return val;
}

long p_val_q(const std::size_t p, const mpq_class& x)
{
  std::size_t num_valuation = p_val_z(p, x.get_num());

  if (num_valuation > 0)
    return static_cast<long>(num_valuation);
  else
    return -static_cast<long>(p_val_z(p, x.get_den()));
}

mpz_class p_pow_z(const std::size_t p, const std::size_t exp)
{
  mpz_class pow;
  mpz_ui_pow_ui(pow.get_mpz_t(), p, exp);
  return pow;
}

mpq_class p_pow_q(const std::size_t p, const long exp)
{
  if (exp >= 0)
    return p_pow_z(p, static_cast<std::size_t>(exp));
  else {
    mpq_class inv_pow = p_pow_z(p, static_cast<std::size_t>(-exp));
    return 1 / inv_pow;
  }
}
