file(GLOB LIB_SOURCES "*.cpp")
list(REMOVE_ITEM LIB_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)

add_library(akss_lib ${LIB_SOURCES} shell.cpp shell.h)
add_executable(akss_main main.cpp)
target_link_libraries(akss_main gmp gmpxx akss_lib)
