cmake_minimum_required (VERSION 2.6)


add_executable(co co.c)
target_link_libraries(co
    -g 
    c3sc
    ${blas_lapack} 
    -lm)

