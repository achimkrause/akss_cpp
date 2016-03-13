#pragma once

#include <gmpxx.h>

#include "types.h"


val_t p_val_z(const mod_t p, const mpz_class& x);
val_t p_val_q(const mod_t p, const mpq_class& x);

mpz_class p_pow_z(const mod_t p, const unsigned long int exp);
mpq_class p_pow_q(const mod_t p, const val_t exp);
