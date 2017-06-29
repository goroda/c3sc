cmake_minimum_required (VERSION 2.6)
project (C3SC)

#create my librray
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

    find_package(c3 REQUIRED)
    #find_package(cdyn REQUIRED)

    set(blas_lapack ${macblas}/Accelerate)
    set(compiler_flags "-std=c99 -O3 -Wall -Weverything -Wno-shorten-64-to-32 -Wno-padded -Wno-missing-prototypes -Wno-format-nonliteral -Wno-disabled-macro-expansion -Wno-missing-noreturn")
    set(CMAKE_C_COMPILER "clang")
else ()
    set(compiler_flags "-g -O2 -Wall -Wextra -std=c99 ")
    # test coverage
    #set(compiler_flags "-g -O0 -W -Wall -Wextra -std=c99 -fprofile-arcs -ftest-coverage")
    set(CMAKE_C_COMPILER "gcc")

    find_package(BLAS REQUIRED)
    find_package(LAPACK REQUIRED)

    #find_package(cdyn REQUIRED)
    #find_package(c3 REQUIRED)


    include_directories(/usr/include)
    set(blas_lapack "-lblas -llapack")
                       
    #find_package(SWIG REQUIRED)
    #include(${SWIG_USE_FILE})
endif()


find_package(OpenMP)
set(compiler_flags "${compiler_flags} ${OpenMP_C_FLAGS}" )

# set(CMAKE_CXX_COMPILER "gcc")

#set(static_libs 1)
set(static_libs 0)

include_directories(${BLAS_INCLUDE_DIR}, ${LAPACK_INCLUDE_DIR})
include_directories(..)
include_directories("src")

set(CMAKE_C_FLAGS "${compiler_flags}")

SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)

include_directories(src/include)

add_subdirectory(src)
add_subdirectory(examples)
add_subdirectory(miscellaneous)
add_subdirectory(test)

install(FILES cmake/c3sc-config.cmake DESTINATION ${main_lib_dest})
install(EXPORT c3sc DESTINATION ${lib_dest})

