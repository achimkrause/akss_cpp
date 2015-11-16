#include <gmpxx.h>

#include "gtest/gtest.h"

#include "../src/p_local.h"

TEST(PLocal, ValuationInt)
{
  EXPECT_EQ(0, p_val_z(3, 1_mpz));
  EXPECT_THROW(p_val_z(3, 0_mpz), std::logic_error);
  EXPECT_EQ(2, p_val_z(5, 25_mpz));
  EXPECT_EQ(1, p_val_z(2, -6_mpz));
}

TEST(PLocal, ValuationRational)
{
  EXPECT_EQ(0, p_val_q(2, 1_mpq / 3));
  EXPECT_EQ(2, p_val_q(5, -25_mpq / 2));
  EXPECT_EQ(-3, p_val_q(2, 4_mpq / 32));
}

TEST(PLocal, PowInt)
{
  EXPECT_EQ(1, p_pow_z(13, 0));
  EXPECT_EQ(25, p_pow_z(5, 2));
}

TEST(PLocal, PowRational)
{
  EXPECT_EQ(25, p_pow_q(5, 2));
  EXPECT_EQ(1/9_mpq, p_pow_q(3, -2));
}
