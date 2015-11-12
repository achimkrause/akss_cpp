cmake_minimum_required(VERSION 2.8.9)
project(akss)

set(CMAKE_CXX_COMPILER "clang++")
set(CMAKE_CXX_FLAGS "-Wall -Werror -std=c++11")
set(CMAKE_EXPORT_COMPILE_COMMANDS "ON")
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/modules/")

include_directories(src)
file(GLOB SOURCES "src/*.cpp")
list(REMOVE_ITEM SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp)
file(GLOB TEST_SOURCES "test/*.cpp")
find_package(GMP REQUIRED)

enable_testing()
find_package(Threads REQUIRED)
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

add_executable(akss ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp ${SOURCES})
target_link_libraries(akss gmp gmpxx)

add_executable(akss_test ${SOURCES} ${TEST_SOURCES})
target_link_libraries(akss_test gmp gmpxx gtest gtest_main pthread)
add_test(NAME akss_test COMMAND akss_test)
