#include "gtest/gtest.h"

#include <string>

#include "common.h"

#include "../src/session.h"

TEST(SessionInit, Parse)
{
  std::string TEST_DATA_PATH = TEST_DATA_DIR + "SessionInitParse/";
  Session session(TEST_DATA_PATH + "ranks.dat",
                  TEST_DATA_PATH + "v_inclusions.dat",
                  TEST_DATA_PATH + "r_operations.dat.",
                  10);
}
