cmake_minimum_required (VERSION 2.6)

add_executable(vtol vtol.c)
target_link_libraries(vtol
    -g 
    c3sc
    ${blas_lapack} 
    -lm)

