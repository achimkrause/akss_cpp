cmake_minimum_required(VERSION 2.8.9)
project(akss)

set(CMAKE_CXX_COMPILER "clang++")
set(CMAKE_CXX_FLAGS "-Wall -Werror -std=c++11")

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/modules/")

include_directories(src)
file(GLOB SOURCES "src/*.cpp")

find_package(GMP REQUIRED)
include_directories(${GOBJECT_INCLUDE_DIR})

add_executable(akss ${SOURCES})
target_link_libraries(akss gmp gmpxx)

