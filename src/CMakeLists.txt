cmake_minimum_required (VERSION 2.6)

#create a libray containing everything!
file(GLOB SRC_FILES *.c)
if (${static_libs} )
    add_library(c3sc ${SRC_FILES})
else ()
    add_library(c3sc SHARED ${SRC_FILES})
    target_link_libraries(c3sc
        ${blas_lapack}
        c3
        -lm
        )
        
endif()

