cmake_minimum_required (VERSION 2.6)


add_executable(cothrust copterposethruststar1_test2.c)
target_link_libraries(cothrust
    -g 
    c3sc
    ${PROJECT_LINK_LIBS}
    ${blas_lapack} 
    -lm)

