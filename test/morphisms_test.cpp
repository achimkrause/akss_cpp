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

TEST(Kernel, Diagonal)
{
	AbelianGroup X(0,2);
	X(0) = 1;
	X(1) = 2;

	AbelianGroup Y(0,1);
	Y(0) = 2;

	MatrixQ f = {{5,2}};

	MatrixQList to_K;
	MatrixQList from_K;

	MatrixQList from_X;
	MatrixQ id = {{1,0},{0,1}};
	from_X.emplace_back(id);

	AbelianGroup K = compute_kernel(5, f, X, Y, MatrixQRefList(), ref(from_X), to_K, from_K);

	EXPECT_EQ(0, K.free_rank());
	ASSERT_EQ(1, K.tor_rank());
	EXPECT_EQ(1, K(0));

	//EXPECT_EQ(0, f*from_K[0]);
}
