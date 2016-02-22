cmake_minimum_required (VERSION 2.6)


add_executable(sim2d sim2d.c)
target_link_libraries(sim2d
    -g 
    c3sc
    ${blas_lapack} 
    -lm)

