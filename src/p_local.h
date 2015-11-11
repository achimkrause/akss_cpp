#pragma once

#include <gmpxx.h>

#include "types.h"

int p_valuation(const std::size_t p, const mpz_class& x);
int p_valuation(const std::size_t p, const mpq_class& x);
