cmake_minimum_required (VERSION 2.6)


add_executable(dp2d dp2d.c)
target_link_libraries(dp2d
    -g 
    c3sc
    ${blas_lapack} 
    -lm)

