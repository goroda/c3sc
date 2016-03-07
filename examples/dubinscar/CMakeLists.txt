cmake_minimum_required (VERSION 2.6)

add_executable(dubinscar dubinscar.c)
target_link_libraries(dubinscar
    -g 
    c3sc
    ${blas_lapack} 
    -lm)

