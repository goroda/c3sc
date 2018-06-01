cmake_minimum_required (VERSION 2.6)


add_executable(rossler rossler.c)
target_link_libraries(rossler
    -g 
    c3sc
    ${blas_lapack} 
    -lm)

