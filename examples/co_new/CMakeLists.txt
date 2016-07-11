cmake_minimum_required (VERSION 2.6)


add_executable(co_new co.c)
target_link_libraries(co_new
    -g 
    c3sc
    ${blas_lapack} 
    -lm)

