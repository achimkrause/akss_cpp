#include <gmpxx.h>

#include "gtest/gtest.h"

#include "../src/p_local.h"

TEST(PLocal, ValuationInt)
{
  ASSERT_EQ(0, p_valuation_int(3, 1_mpz));
  ASSERT_THROW(p_valuation_int(3, 0_mpz), std::logic_error);
  ASSERT_EQ(2, p_valuation_int(5, 25_mpz));
  ASSERT_EQ(1, p_valuation_int(2, -6_mpz));
}

TEST(PLocal, ValuationRational)
{
  ASSERT_EQ(0, p_valuation(2, 1_mpq/3));
  ASSERT_EQ(2, p_valuation(5, -25_mpq/2));
  ASSERT_EQ(-3, p_valuation(2, 4_mpq/32));
}
