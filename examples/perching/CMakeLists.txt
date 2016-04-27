cmake_minimum_required (VERSION 2.6)


add_executable(perch perch.c)
target_link_libraries(perch
    -g 
    c3sc
    ${blas_lapack} 
    -lm)

