cmake_minimum_required (VERSION 2.6)


add_executable(markov2d markov2d.c)
target_link_libraries(markov2d
    -g 
    c3sc
    ${blas_lapack} 
    -lm)

