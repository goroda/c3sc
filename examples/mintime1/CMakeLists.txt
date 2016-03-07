cmake_minimum_required (VERSION 2.6)


add_executable(mintime1 mintime1.c)
target_link_libraries(mintime1
    -g 
    c3sc
    ${blas_lapack} 
    -lm)

