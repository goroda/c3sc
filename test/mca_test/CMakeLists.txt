cmake_minimum_required (VERSION 2.6)

file(GLOB SRC_FILES *.c)
add_executable(mca_test ${SRC_FILES})
target_link_libraries(mca_test
    -g
    c3sc
    cutest 
    ${blas_lapack} 
    -lm)


