cmake_minimum_required (VERSION 2.6)


add_executable(scar5d scar.c)
target_link_libraries(scar5d
    -g 
    c3sc
    ${blas_lapack} 
    -lm)

