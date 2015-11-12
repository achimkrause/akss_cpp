#include <gmpxx.h>

#include "gtest/gtest.h"

#include "../src/matrix.h"
#include "../src/smith.h"

TEST(SmithReduceP, Empty)
{
  MatrixQ f = {};

  auto to_X = MatrixQRefList();
  auto from_X = MatrixQRefList();
  auto to_Y = MatrixQRefList();
  auto from_Y = MatrixQRefList();

  smith_reduce_p(2, f, to_X, from_X, to_Y, from_Y);

  EXPECT_EQ(MatrixQ({}), f);
}

TEST(SmithReduceP, Diagonal)
{
  MatrixQ f = MatrixQ::identity(3);

  auto to_X = MatrixQRefList();
  auto from_X = MatrixQRefList();
  auto to_Y = MatrixQRefList();
  auto from_Y = MatrixQRefList();

  smith_reduce_p(2, f, to_X, from_X, to_Y, from_Y);

  EXPECT_EQ(MatrixQ::identity(3), f);
}

TEST(SmithReduceP, AntiDiagonal)
{
  MatrixQ f = {{0, 0, 1}, {0, 1, 0}, {1, 0, 0}};

  auto to_X = MatrixQRefList();
  auto from_X = MatrixQRefList();
  auto to_Y = MatrixQRefList();
  auto from_Y = MatrixQRefList();

  smith_reduce_p(2, f, to_X, from_X, to_Y, from_Y);

  EXPECT_EQ(MatrixQ({{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}), f);
}
