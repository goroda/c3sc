cmake_minimum_required (VERSION 2.6)


add_executable(lqg_d lqgnd.c)
target_link_libraries(lqg_d
    -g 
    c3sc
    ${blas_lapack} 
    -lm)
# target_link_libraries(lqg2d_new PUBLIC c3)
