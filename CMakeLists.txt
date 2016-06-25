cmake_minimum_required(VERSION 2.8.9)
project(akss)

set(CMAKE_CXX_COMPILER "clang++")
set(CMAKE_CXX_FLAGS "-O0 -Weverything -Wno-exit-time-destructors -Wno-weak-vtables -Wno-padded -Wno-c++98-compat-bind-to-temporary-copy -Wno-global-constructors -Wno-c++98-compat -std=c++11 -g")
set(CMAKE_EXPORT_COMPILE_COMMANDS "ON")
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/modules/")

include_directories(src)
find_package(GMP REQUIRED)

enable_testing()
add_test(NAME akss_test COMMAND akss_test "${CMAKE_CURRENT_LIST_DIR}/test_data/")

add_subdirectory(src)
add_subdirectory(test)
