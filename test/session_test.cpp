#include "gtest/gtest.h"

#include <string>

#include "common.h"

#include "../src/session.h"

TEST(SessionInit, Parse)
{
  std::string TEST_DATA_PATH = TEST_DATA_DIR + "SessionInitParse/";
  Session session(2, TEST_DATA_PATH + "ranks.dat",
                  TEST_DATA_PATH + "v_inclusions.dat",
                  TEST_DATA_PATH + "r_operations.dat.",
                  10);
}

TEST(SessionInit, Step)
{
  std::string TEST_DATA_PATH = TEST_DATA_DIR + "SessionInitParse/";
  Session session(2, TEST_DATA_PATH + "ranks.dat",
                  TEST_DATA_PATH + "v_inclusions.dat",
                  TEST_DATA_PATH + "r_operations.dat.",
                  10);
  session.step();
}

TEST(SessionInit, TwoSteps)
{
  std::string TEST_DATA_PATH = TEST_DATA_DIR + "SessionInitParse/";
  Session session(2, TEST_DATA_PATH + "ranks.dat",
                  TEST_DATA_PATH + "v_inclusions.dat",
                  TEST_DATA_PATH + "r_operations.dat.",
                  10);
  session.step();
  session.step();
}

TEST(SessionInit, ThreeSteps)
{
  std::string TEST_DATA_PATH = TEST_DATA_DIR + "SessionInitParse/";
  Session session(2, TEST_DATA_PATH + "ranks.dat",
                  TEST_DATA_PATH + "v_inclusions.dat",
                  TEST_DATA_PATH + "r_operations.dat.",
                  10);
  session.step();
  session.step();
  session.step();
}

TEST(SessionInit, TenSteps)
{
  std::string TEST_DATA_PATH = TEST_DATA_DIR + "SessionInitParse/";
  Session session(2, TEST_DATA_PATH + "ranks.dat",
                  TEST_DATA_PATH + "v_inclusions.dat",
                  TEST_DATA_PATH + "r_operations.dat.",
                  10);
  for(int i=0; i<10; i++){
    session.step();
  };
}

