#include <exception>

#include <gmpxx.h>
#include "gtest/gtest.h"

#include "../src/matrix.h"

TEST(Matrix, Properties)
{
  MatrixQ A(3, 2);

  ASSERT_EQ(3, A.height());
  ASSERT_EQ(2, A.width());
}

TEST(Matrix, InitializerList)
{
  MatrixQ A = {{1, 2, 3}, {4, 5, 6}};

  ASSERT_EQ(2, A(0, 1));
  ASSERT_EQ(2, A.height());
  ASSERT_EQ(3, A.width());
}

TEST(Matrix, Comparison)
{
  MatrixQ A = {{1, 2, 3}, {4, 5, 6}};
  MatrixQ B = {{1, 2, 3}, {4, 5, 6}};
  MatrixQ C = {{1, 2, 3}, {4, 6, 6}};
  MatrixQ D = {{1, 2, 3, 3}, {4, 5, 6, 6}};

  ASSERT_EQ(A, B);
  ASSERT_NE(A, C);
  ASSERT_NE(A, D);
}

TEST(Matrix, RowAdd)
{
  MatrixQ A = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
  MatrixQ B = {{1, 2, 3}, {4, 5, 6}, {9, 12, 15}};
  A.row_add(0, 2, 2_mpq);
  ASSERT_EQ(B, A);
}

TEST(Matrix, Composition)
{
  MatrixQ A = {{1, 0, 1}, {0, 1, 1}};

  MatrixQ B = {{1, 0}, {0, 1}, {1, 1}};

  MatrixQ C = A * B;

  MatrixQ C_ref = {{2, 1}, {1, 2}};

  EXPECT_EQ(C_ref, C);
}

TEST(MatrixSlice, DimensionMismatch)
{
  MatrixQ A(2, 2);
  MatrixQ B(3, 3);

  EXPECT_THROW(A(0, 0, 2, 2) = B(0, 0, 1, 1), std::logic_error);
}

TEST(MatrixSlice, Aliasing)
{
  MatrixQ A(5, 5);

  EXPECT_THROW(A(0, 0, 2, 2) = A(1, 1, 2, 2), std::logic_error);
  EXPECT_NO_THROW(A(0, 0, 2, 2) = A(2, 2, 2, 2));
  EXPECT_THROW(A(3, 3, 2, 2) = A(3, 2, 2, 2), std::logic_error);
  EXPECT_NO_THROW(A(1, 3, 2, 2) = A(1, 1, 2, 2));
}
