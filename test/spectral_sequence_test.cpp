#include "gtest/gtest.h"

#include "../src/spectral_sequence.h"

TEST(TrigradedIndex, Equality)
{
  TrigradedIndex ind_1(3, 2, 1);
  TrigradedIndex ind_2(3, 2, 1);

  EXPECT_EQ(ind_1, ind_2);
}

TEST(TrigradedIndex, LessThan)
{
  TrigradedIndex ind_1(4, -1, 11);
  TrigradedIndex ind_2(2, 2, 10);
  TrigradedIndex ind_3(2, 3, 4);
  TrigradedIndex ind_4(2, 5, 5);
  TrigradedIndex ind_5(2, 5, 6);

  EXPECT_LT(ind_1, ind_2);
  EXPECT_LT(ind_2, ind_3);
  EXPECT_LT(ind_3, ind_4);
  EXPECT_LT(ind_4, ind_5);
}

TEST(GroupSequence, GetGroup)
{
  GroupSequence seq(2, AbelianGroup(2, 0));
  EXPECT_THROW(seq.get_group(1), std::logic_error);
  EXPECT_THROW(seq.get_group(3), std::logic_error);
  EXPECT_EQ(2, seq.get_group(2).free_rank());
  EXPECT_EQ(0, seq.get_group(2).tor_rank());
}

TEST(GroupSequence, GetMatrix)
{
  GroupSequence seq(2, AbelianGroup(2, 0));
  EXPECT_THROW(seq.get_matrix(1), std::logic_error);
  EXPECT_THROW(seq.get_matrix(3), std::logic_error);
  EXPECT_EQ(MatrixQ::identity(2), seq.get_matrix(2));
}
