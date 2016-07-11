cmake_minimum_required (VERSION 2.6)

add_executable(dubinscar_new dubinscar.c)
target_link_libraries(dubinscar_new
    -g 
    c3sc
    ${blas_lapack} 
    -lm)

