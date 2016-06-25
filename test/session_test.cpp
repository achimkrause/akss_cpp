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

TEST(SessionInit, Step)
{
  std::string TEST_DATA_PATH = TEST_DATA_DIR + "SessionInitParse/";
  Session session(TEST_DATA_PATH + "ranks.dat",
                  TEST_DATA_PATH + "v_inclusions.dat",
                  TEST_DATA_PATH + "r_operations.dat.",
                  10);
  session.step();
}

TEST(SessionInit, TwoSteps)
{
  std::string TEST_DATA_PATH = TEST_DATA_DIR + "SessionInitParse/";
  Session session(TEST_DATA_PATH + "ranks.dat",
                  TEST_DATA_PATH + "v_inclusions.dat",
                  TEST_DATA_PATH + "r_operations.dat.",
                  10);
  session.step();
  //std::cerr << session.get_sequence().get_e_2(TrigradedIndex(0,1,1)).tor_rank() << "\n";
  session.step();
}
