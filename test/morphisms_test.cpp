#include "gtest/gtest.h"

#include "../src/matrix.h"
#include "../src/morphisms.h"

TEST(Cokernel, Diagonal)
{
  AbelianGroup Y(3, 0);
  MatrixQ f = {{0, 0, 0}, {0, 2, 0}, {0, 0, 3}};

  MatrixQList to_C;
  MatrixQList from_C;
  AbelianGroup C = compute_cokernel(3, f, Y, MatrixQRefList(), MatrixQRefList(),
                                    to_C, from_C);

  EXPECT_EQ(1, C.free_rank());
  ASSERT_EQ(1, C.tor_rank());
  EXPECT_EQ(1, C(0));
}
