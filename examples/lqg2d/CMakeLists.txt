cmake_minimum_required (VERSION 2.6)


add_executable(lqg2d lqg2d.c)
target_link_libraries(lqg2d
    -g 
    c3sc
    ${blas_lapack} 
    -lm)

