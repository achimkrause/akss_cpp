#include "gtest/gtest.h"

#include <string>

#include "common.h"

std::string TEST_DATA_DIR;

int main(int argc, char **argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  assert(argc == 2);
  TEST_DATA_DIR = argv[1];
  int ret = RUN_ALL_TESTS();
  return ret;
}
