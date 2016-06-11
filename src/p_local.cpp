#include <limits>

#include "p_local.h"

val_t p_val_z(const mod_t p, const mpz_class& x)
{
  if (x == 0) return std::numeric_limits<val_t>::max();

  val_t val = 0;
  mpz_class remainder = x;

  while (mpz_divisible_ui_p(remainder.get_mpz_t(), p)) {
    val++;
    mpz_divexact_ui(remainder.get_mpz_t(), remainder.get_mpz_t(), p);
  }

  return val;
}

val_t p_val_q(const mod_t p, const mpq_class& x)
{
  val_t num_valuation = p_val_z(p, x.get_num());

  if (num_valuation > 0)
    return num_valuation;
  else
    return -p_val_z(p, x.get_den());
}

mpz_class p_pow_z(const mod_t p, const u_val_t exp)
{
  if (exp == std::numeric_limits<val_t>::max())
    return 0_mpz;

  mpz_class pow;
  mpz_ui_pow_ui(pow.get_mpz_t(), p, exp);
  return pow;
}

mpq_class p_pow_q(const mod_t p, const val_t exp)
{
  if (exp >= 0)
    return p_pow_z(p, static_cast<u_val_t>(exp));
  else {
    mpq_class inv_pow = p_pow_z(p, static_cast<u_val_t>(-exp));
    return 1 / inv_pow;
  }
}
