cmake_minimum_required (VERSION 2.6)

file(GLOB SRC_FILES *.c)
add_executable(hashgrid_test ${SRC_FILES})
target_link_libraries(hashgrid_test
    -g
    c3sc
    cutest 
    ${blas_lapack} 
    -lm)


