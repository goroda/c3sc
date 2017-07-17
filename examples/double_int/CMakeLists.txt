cmake_minimum_required (VERSION 2.6)


add_executable(double_int double_int.c)
target_link_libraries(double_int
    -g 
    c3sc
    ${blas_lapack} 
    -lm
    c3 )

