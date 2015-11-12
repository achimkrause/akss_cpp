#include <gmpxx.h>

#include "gtest/gtest.h"

#include "../src/matrix.h"
#include "../src/smith.h"

TEST(SmithReduceP, Empty)
{
  MatrixQ f = {};

  auto to_X = MatrixQList();
  auto from_X = MatrixQList();
  auto to_Y = MatrixQList();
  auto from_Y = MatrixQList();

  smith_reduce_p(2, f, to_X, from_X, to_Y, from_Y);

  ASSERT_EQ(MatrixQ({}), f);
}

TEST(SmithReduceP, Diagonal)
{
  MatrixQ f = MatrixQ::identity(3);

  auto to_X = MatrixQList();
  auto from_X = MatrixQList();
  auto to_Y = MatrixQList();
  auto from_Y = MatrixQList();

  smith_reduce_p(2, f, to_X, from_X, to_Y, from_Y);

  ASSERT_EQ(MatrixQ::identity(3), f);
}
