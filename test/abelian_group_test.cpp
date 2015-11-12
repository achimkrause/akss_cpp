#include "gtest/gtest.h"

#include "../src/abelian_group.h"
#include "../src/matrix.h"

TEST(AbelianGroup, TorsionMatrix)
{
  AbelianGroup X(2, 3);
  X(0) = 2;
  X(1) = 1;
  X(2) = 3;

  EXPECT_EQ(MatrixQ({{9, 0, 0}, {0, 3, 0}, {0, 0, 27}}), X.torsion_matrix(3));
}
