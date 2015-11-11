#include <gmpxx.h>

#include "gtest/gtest.h"

#include "../src/p_local.h"

TEST(PLocal, ValuationInt)
{
  ASSERT_EQ(0, p_valuation(3, 1_mpz));
  ASSERT_THROW(p_valuation(3, 0_mpz), std::logic_error);
  ASSERT_EQ(2, p_valuation(5, 25_mpz));
}
