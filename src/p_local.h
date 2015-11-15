#pragma once

#include <gmpxx.h>

std::size_t p_val_z(const std::size_t p, const mpz_class& x);
long p_val_q(const std::size_t p, const mpq_class& x);

mpz_class p_pow_z(const std::size_t p, const std::size_t exp);
mpq_class p_pow_q(const std::size_t p, const long exp);
