#include <exception>

#include <gmpxx.h>
#include "gtest/gtest.h"

#include "../src/matrix.h"

TEST(Matrix, Properties)
{
  MatrixQ A(3, 2);

  EXPECT_EQ(3, A.height());
  EXPECT_EQ(2, A.width());
}

TEST(Matrix, InitializerList)
{
  MatrixQ A = {{1, 2, 3}, {4, 5, 6}};

  EXPECT_EQ(2, A(0, 1));
  EXPECT_EQ(2, A.height());
  EXPECT_EQ(3, A.width());
}

TEST(Matrix, Comparison)
{
  MatrixQ A = {{1, 2, 3}, {4, 5, 6}};
  MatrixQ B = {{1, 2, 3}, {4, 5, 6}};
  MatrixQ C = {{1, 2, 3}, {4, 6, 6}};
  MatrixQ D = {{1, 2, 3, 3}, {4, 5, 6, 6}};

  EXPECT_EQ(A, B);
  EXPECT_NE(A, C);
  EXPECT_NE(A, D);
}

TEST(Matrix, RowAdd)
{
  MatrixQ A = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
  MatrixQ B = {{1, 2, 3}, {4, 5, 6}, {9, 12, 15}};
  EXPECT_EQ(B, A.row_add(0, 2, 2_mpq));
}
TEST(Matrix, RowMul)
{
  MatrixQ A = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
  MatrixQ B = {{1_mpq / 2, 1_mpq, 3_mpq / 2}, {4, 5, 6}, {7, 8, 9}};
  EXPECT_EQ(B, A.row_mul(0, 1 / 2_mpq));
}

TEST(Matrix, RowSwap)
{
  MatrixQ A = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
  MatrixQ B = {{1, 2, 3}, {7, 8, 9}, {4, 5, 6}};
  EXPECT_EQ(B, A.row_swap(1, 2));

  MatrixQ C = A;
  EXPECT_EQ(C, A.row_swap(2, 2));
}
TEST(Matrix, ColAdd)
{
  MatrixQ A = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
  MatrixQ B = {{-1, 2, 3}, {-1, 5, 6}, {-1, 8, 9}};
  EXPECT_EQ(B, A.col_add(1, 0, -1_mpq));
}
TEST(Matrix, ColMul)
{
  MatrixQ A = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
  MatrixQ B = {{1_mpq / 2, 2, 3}, {2, 5, 6}, {7 / 2_mpq, 8, 9}};
  EXPECT_EQ(B, A.col_mul(0, 1 / 2_mpq));
}

TEST(Matrix, ColSwap)
{
  MatrixQ A = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
  MatrixQ B = {{2, 1, 3}, {5, 4, 6}, {8, 7, 9}};
  EXPECT_EQ(B, A.col_swap(0, 1));

  MatrixQ C = A;
  EXPECT_EQ(C, A.col_swap(0, 0));
}

TEST(Matrix, BasisVectorsAdd)
{
  MatrixQ A = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
  MatrixQ B(A);

  MatrixQRefList to_X = {A};
  MatrixQRefList from_X = {B};

  basis_vectors_add(to_X, from_X, 0, 2, 2_mpq);

  EXPECT_EQ(MatrixQ({{-13, -14, -15}, {4, 5, 6}, {7, 8, 9}}), A);
  EXPECT_EQ(MatrixQ({{1, 2, 5}, {4, 5, 14}, {7, 8, 23}}), B);
}

TEST(Matrix, BasisVectorsMul)
{
  MatrixQ A = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
  MatrixQ B(A);

  MatrixQRefList to_X = {A};
  MatrixQRefList from_X = {B};

  basis_vectors_mul(to_X, from_X, 1, 2_mpq);

  EXPECT_EQ(MatrixQ({{1, 2, 3}, {2, 5 / 2_mpq, 3}, {7, 8, 9}}), A);
  EXPECT_EQ(MatrixQ({{1, 4, 3}, {4, 10, 6}, {7, 16, 9}}), B);
}

TEST(Matrix, BasisVectorsSwap)
{
  MatrixQ A = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
  MatrixQ B(A);

  MatrixQRefList to_X = {A};
  MatrixQRefList from_X = {B};

  basis_vectors_swap(to_X, from_X, 1, 2);

  EXPECT_EQ(MatrixQ({{1, 2, 3}, {7, 8, 9}, {4, 5, 6}}), A);
  EXPECT_EQ(MatrixQ({{1, 3, 2}, {4, 6, 5}, {7, 9, 8}}), B);
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
