#include <exception>

#include "p_local.h"

int p_valuation_int(const std::size_t p, const mpz_class& x)
{
  if (x == 0) throw std::logic_error("p_valuation: x=0");

  int valuation = 0;
  mpz_class remainder = x;

  while (mpz_divisible_ui_p(remainder.get_mpz_t(), p)) {
    valuation++;
    mpz_divexact_ui(remainder.get_mpz_t(), remainder.get_mpz_t(), p);
  }

  return valuation;
}

int p_valuation(const std::size_t p, const mpq_class& x)
{
  int num_valuation = p_valuation_int(p, x.get_num());

  if (num_valuation > 0)
    return num_valuation;
  else
    return -p_valuation_int(p, x.get_den());
}

mpz_class p_pow(const std::size_t p, const std::size_t exp) {
  mpz_class pow;
  mpz_ui_pow_ui(pow.get_mpz_t(), p, exp);
  return pow;
}
