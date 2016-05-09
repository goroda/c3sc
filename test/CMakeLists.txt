cmake_minimum_required (VERSION 2.6)

#include_directories(/opt/local/include)

#add_executable(allTests CuTest.c vec_test.c mat_test.c 
#                        linalg_test.c tensor_test.c tt_test.c cross_test.c AllMyTests.c)
#TARGET_LINK_LIBRARIES(allTests -g array -lm)


include_directories(CuTest)
add_subdirectory(CuTest)

add_subdirectory(mca_test)
add_subdirectory(hashgrid_test)



