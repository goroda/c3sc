cmake_minimum_required (VERSION 2.6)


add_executable(lqg2d_new lqg2d.c)
target_link_libraries(lqg2d_new
    -g 
    c3sc
    ${blas_lapack} 
    -lm
    c3 )
# target_link_libraries(lqg2d_new PUBLIC c3)

