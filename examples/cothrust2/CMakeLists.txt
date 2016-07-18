cmake_minimum_required (VERSION 2.6)


add_executable(cothrust_nostar copterposethrust.c)
target_link_libraries(cothrust_nostar
    -g 
    c3sc
    ${PROJECT_LINK_LIBS}
    ${blas_lapack} 
    -lm)
