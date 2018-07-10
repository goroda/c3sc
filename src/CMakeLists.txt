cmake_minimum_required (VERSION 2.6)

set(header bellman.h
           boundary.h
           c3sc.h
           dynamics.h
           hashgrid.h
           nodeutil.h
           util.h
           valuefunc.h)

set(src bellman.c
        boundary.c
        dynamics.c
        hashgrid.c
        nodeutil.c
        util.c
        valuefunc.c)

add_library(c3sc SHARED ${header} ${src})
target_link_libraries(c3sc
  ${blas_lapack}
  -lm
  ${c3lib}
  ${cdynlib}
  )
target_include_directories(c3sc PUBLIC 
                     $<BUILD_INTERFACE:${C3SC_SOURCE_DIR}/include>
                     $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
                     $<INSTALL_INTERFACE:${include_dest}>
                     $<INSTALL_INTERFACE:${lib_dest}>)
                   
# install(TARGETS c3sc EXPORT c3sc DESTINATION "${lib_dest}")
# install(FILES ${header} DESTINATION "${include_dest}")

install(TARGETS c3sc DESTINATION lib)
install(FILES ${header} DESTINATION include/c3sc)

