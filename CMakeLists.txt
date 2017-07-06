cmake_minimum_required (VERSION 2.6)
project (C3SC)

#create my libray
set(C3SC_VERSION_MAJOR 1 CACHE STRING "major version" FORCE)
set(C3SC_VERSION_MINOR 0 CACHE STRING "minor version" FORCE)
set(C3SC_VERSION "${C3SC_VERSION_MAJOR}.${C3SC_VERSION_MINOR}" CACHE STRING "version" FORCE)

option(C3SC_DEBUG_MODE "whether or not debug mode is activated" OFF)

set(include_dest "include/c3sc")
set(main_lib_dest "lib/c3sc")
set(lib_dest "${main_lib_dest}/${CMAKE_BUILD_TYPE}")

if (APPLE)
    option(USE_CLANG "build application with clang" ON)
    include_directories(/opt/local/include)
    find_library(macblas
        NAMES ACCELERATE
        PATHS ${CMAKE_OSX_SYSROOT}/System/Library
        PATH_SUFFIXES Frameworks
        NO_DEFAULT_PATH)

    set(blas_lapack ${macblas}/Accelerate)
    set(compiler_flags "-std=c99 -O2 -Wall -Weverything -Wno-shorten-64-to-32 -Wno-padded -Wno-missing-prototypes -Wno-format-nonliteral -Wno-disabled-macro-expansion -Wno-missing-noreturn -fopenmp=libomp -Wno-reserved-id-macro")
    set(CMAKE_C_COMPILER "/opt/local/bin/clang-mp-4.0")
    #set(CMAKE_C_COMPILER "clang")

    set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} /opt/local/lib)
else ()
    set(compiler_flags "-g -O2 -ftree-vectorize -Wall -Wextra -std=c99 -fno-omit-frame-pointer")

    set(CMAKE_C_COMPILER "gcc")

    find_package(BLAS REQUIRED)
    find_package(LAPACK REQUIRED)

    include_directories(/usr/include)
    set(blas_lapack "-lblas -llapack")

endif()

if ("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
   # using GCC
    find_package(OpenMP)
    if (OPENMP_FOUND)
       set(compiler_flags "${compiler_flags} ${OpenMP_C_FLAGS}" )
    endif()
endif()



include_directories(${BLAS_INCLUDE_DIR}, ${LAPACK_INCLUDE_DIR})
include_directories(src)
include_directories(..)

set(CMAKE_C_FLAGS "${compiler_flags}")

SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)

add_subdirectory(src)
add_subdirectory(examples)
add_subdirectory(miscellaneous)
add_subdirectory(test)

install(FILES cmake/c3sc-config.cmake DESTINATION ${main_lib_dest})
install(EXPORT c3sc DESTINATION ${lib_dest})

