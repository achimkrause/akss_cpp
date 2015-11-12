#include <gmpxx.h>

#include "gtest/gtest.h"

#include "../src/p_local.h"

TEST(PLocal, ValuationInt)
{
  EXPECT_EQ(0, p_valuation_int(3, 1_mpz));
  EXPECT_THROW(p_valuation_int(3, 0_mpz), std::logic_error);
  EXPECT_EQ(2, p_valuation_int(5, 25_mpz));
  EXPECT_EQ(1, p_valuation_int(2, -6_mpz));
}

TEST(PLocal, ValuationRational)
{
  EXPECT_EQ(0, p_valuation(2, 1_mpq/3));
  EXPECT_EQ(2, p_valuation(5, -25_mpq/2));
  EXPECT_EQ(-3, p_valuation(2, 4_mpq/32));
}
