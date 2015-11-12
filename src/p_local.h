#pragma once

#include <gmpxx.h>

#include "types.h"

int p_valuation_int(const std::size_t p, const mpz_class& x);
int p_valuation(const std::size_t p, const mpq_class& x);

mpz_class p_pow(const std::size_t p, const std::size_t exp);
