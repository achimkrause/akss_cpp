file(GLOB SOURCES "*.cpp")

find_package(GTest REQUIRED)

add_executable(akss_test ${SOURCES})
target_link_libraries(akss_test gmp gmpxx gtest gtest_main pthread akss_lib)
