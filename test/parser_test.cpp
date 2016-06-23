#include <gmpxx.h>
#include <string>
#include <sstream>
#include "gtest/gtest.h"

#include "../src/parser.h"

TEST(Parser, Read)
{
  std::stringstream read_in("Test");
  EXPECT_EQ("Test", read(read_in, 4));
}

TEST(Parser, ParseMpzClass)
{
  std::stringstream read_in("-12412312");
  mpz_class val;
  EXPECT_TRUE(parse_mpz_class(read_in, val));
  EXPECT_EQ(-12412312_mpz, val);
}

TEST(Parser, ParseMpqClass)
{
  {
    std::stringstream read_in("-45/90");
    mpq_class val;
    EXPECT_TRUE(parse_mpq_class(read_in, val));
    EXPECT_EQ(-1/2_mpq, val);
  }

  {
    std::stringstream read_in("67");
    mpq_class val;
    EXPECT_TRUE(parse_mpq_class(read_in, val));
    EXPECT_EQ(67_mpq, val);
  }
}

TEST(Parser, EatWhitespace)
{
  std::stringstream read_in("  \r\n\t  abc");
  eat_whitespace(read_in);
  std::string abc;
  read_in >> abc;
  EXPECT_EQ("abc", abc);
}

TEST(Parser, ParseMatrix)
{
  std::stringstream read_in("2 3 1 2/5 3 -4 5 -6/2");
  MatrixQ result;
  EXPECT_TRUE(parse_matrix(read_in, result));
  MatrixQ expected = {{1_mpq, 2/5_mpq, 3_mpq}, {-4_mpq, 5_mpq, -6/2_mpq}};
  EXPECT_EQ(expected, result);
}
