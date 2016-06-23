#include "gtest/gtest.h"

#include "../src/matrix.h"
#include "../src/morphisms.h"

TEST(Cokernel, Diagonal)
{
  AbelianGroup Y(3, 0);
  MatrixQ f = {{0, 0, 0}, {0, 2, 0}, {0, 0, 3}};

  MatrixQList to_C;
  MatrixQList from_C;
  GroupWithMorphisms C =
      compute_cokernel(3, f, Y, MatrixQRefList(), MatrixQRefList());

  EXPECT_EQ(1, C.group.free_rank());
  ASSERT_EQ(1, C.group.tor_rank());
  EXPECT_EQ(1, C.group(0));
}

TEST(Kernel, Diagonal)
{
  AbelianGroup X(0, 2);
  X(0) = 1;
  X(1) = 2;

  AbelianGroup Y(0, 1);
  Y(0) = 2;

  MatrixQ f = {{5, 2}};

  MatrixQList to_K;
  MatrixQList from_K;

  MatrixQList from_X;
  MatrixQ id = {{1, 0}, {0, 1}};
  from_X.emplace_back(id);

  GroupWithMorphisms K =
      compute_kernel(5, f, X, Y, MatrixQRefList(), ref(from_X));

  EXPECT_EQ(0, K.group.free_rank());
  ASSERT_EQ(1, K.group.tor_rank());
  EXPECT_EQ(1, K.group(0));
}

TEST(Image, Diagonal)
{
  AbelianGroup X(0, 2);
  X(0) = 2;
  X(1) = 2;

  AbelianGroup Y(0, 2);
  Y(0) = 2;
  Y(1) = 2;

  MatrixQ f = {{0, 6}, {0, 0}};

  // MatrixQList from_X = {MatrixQ::identity(f.width())};
  // GroupWithMorphisms K = compute_kernel(3, f, X, Y, MatrixQRefList(), ref(from_X));
  // EXPECT_EQ(0, K.group.free_rank());
  // EXPECT_EQ(2, K.group.tor_rank());
  // EXPECT_EQ(1 + 2, K.group(0) + K.group(1));
  // EXPECT_EQ(1 * 2, K.group(0) * K.group(1));
  // MatrixQ expected = {{1, 0}, {0, 3}};
  // EXPECT_EQ(expected, K.maps_from[0]);

  GroupWithMorphisms I = compute_image(3, f, X, Y);

  EXPECT_EQ(0, I.group.free_rank());
  ASSERT_EQ(1, I.group.tor_rank());
  EXPECT_EQ(1, I.group(0));
}

TEST(Morphism, LiftFromFree)
{
  AbelianGroup Y(0, 2);
  Y(0) = 1;
  Y(1) = 3; //Y = Z/2 + Z/8

  MatrixQ f = {{1,4},{3,0}};  //f: Z+Z -> Z/2 + Z/8 is lifted by {{3,0}} against map.
  MatrixQ map = {{1},{1}};  //map: X -> Y is diagonal n -> (n,n)

  MatrixQ lift = lift_from_free(2, f, map, Y);

  MatrixQ expected = {{3,0}};
  EXPECT_EQ(expected, lift);
}

TEST(Morphism, Equal)
{
  MatrixQ f = {{2, 1}, {5, 2}};
  MatrixQ g = {{0, 1}, {1, 2}};
  AbelianGroup Y(0, 2);
  Y(0) = 1;
  Y(1) = 2;

  EXPECT_TRUE(morphism_equal(2, f, g, Y));
}

TEST(Morphism, Zero)
{
  AbelianGroup Y(0, 2);
  Y(0) = 1;
  Y(1) = 2;

  MatrixQ h = {{2, 4, 0}, {12, 0, 8}};
  EXPECT_TRUE(morphism_zero(2, h, Y));
}
