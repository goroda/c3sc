cmake_minimum_required (VERSION 2.6)

add_executable(lqg2dr lqg2dr.c)
target_link_libraries(lqg2dr
    -g 
    c3sc
    ${blas_lapack} 
    -lm)

