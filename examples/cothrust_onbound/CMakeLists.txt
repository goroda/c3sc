cmake_minimum_required (VERSION 2.6)


add_executable(cobound co.c)
target_link_libraries(cobound
    -g 
    c3sc
    ${blas_lapack} 
    -lm)

