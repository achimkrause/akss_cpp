#include <exception>

#include "gtest/gtest.h"

#include "matrix.h"

TEST(Matrix, Composition)
{
  MatrixQ A(2, 3);

  A(0, 0) = 1;
  A(0, 1) = 0;
  A(0, 2) = 1;
  A(1, 0) = 0;
  A(1, 1) = 1;
  A(1, 2) = 1;

  MatrixQ B(3, 2);

  B(0, 0) = 1;
  B(0, 1) = 0;
  B(1, 0) = 0;
  B(1, 1) = 1;
  B(2, 0) = 1;
  B(2, 1) = 1;

  MatrixQ C = A * B;

  MatrixQ C_ref(2, 2);

  C_ref(0, 0) = 2;
  C_ref(0, 1) = 1;
  C_ref(1, 0) = 1;
  C_ref(1, 1) = 2;

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
