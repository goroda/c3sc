cmake_minimum_required (VERSION 3.0)
project (C3SC)

#create my libray
set(C3SC_VERSION_MAJOR 1 CACHE STRING "major version" FORCE)
set(C3SC_VERSION_MINOR 0 CACHE STRING "minor version" FORCE)
set(C3SC_VERSION "${C3SC_VERSION_MAJOR}.${C3SC_VERSION_MINOR}" CACHE STRING "version" FORCE)

option(C3SC_DEBUG_MODE "whether or not debug mode is activated" OFF)

include(c3sc.cmake)
c3_always_full_rpath()

# set(include_dest "include/c3sc")
# set(main_lib_dest "lib/c3sc")
# set(lib_dest "${main_lib_dest}/${CMAKE_BUILD_TYPE}")

if (APPLE)
  set(CMAKE_MACOSX_RPATH 1)

  option(USE_CLANG "build application with clang" ON)
  # set(compiler_flags "-std=c99 -O2 -Wall -Weverything -Wno-shorten-64-to-32 -Wno-padded -Wno-missing-prototypes -Wno-format-nonliteral -Wno-disabled-macro-expansion -Wno-missing-noreturn -fopenmp=libomp -Wno-reserved-id-macro")
    # set(CMAKE_C_COMPILER "/opt/local/bin/clang-mp-4.0")

    set(compiler_flags "-std=c99 -O2 -Wall -Weverything -Wno-shorten-64-to-32 -Wno-padded -Wno-missing-prototypes -Wno-format-nonliteral -Wno-disabled-macro-expansion -Wno-missing-noreturn -Wno-reserved-id-macro")

    include_directories(/opt/local/include)
    find_library(macblas
      NAMES ACCELERATE
      PATHS ${CMAKE_OSX_SYSROOT}/System/Library
      PATH_SUFFIXES Frameworks
      NO_DEFAULT_PATH)

    set(blas_lapack ${macblas}/Accelerate)

else ()
    set(compiler_flags "-g -O2 -ftree-vectorize -Wall -Wextra -std=c99 -fno-omit-frame-pointer")

    set(CMAKE_C_COMPILER "gcc")
    include_directories(/usr/include)
    
    find_package(BLAS REQUIRED)
    find_package(LAPACK REQUIRED)

    include_directories(${BLAS_INCLUDE_DIR}, ${LAPACK_INCLUDE_DIR})
    set(blas_lapack ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

endif()

option(USE_OPEN_MP "Use OPEN_MP for parallelization" ON)
if (USE_OPEN_MP)
  if ("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
    # using GCC
    find_package(OpenMP)
    if (OPENMP_FOUND)
      set(compiler_flags "${compiler_flags} ${OpenMP_C_FLAGS}" )
    endif()
  endif()
endif()


include_directories(src)

set(C3_INCLUDE_DIR "" CACHE FILEPATH "Path to Compressed Continuous Computation (C3) includes")

if(NOT EXISTS "${C3_INCLUDE_DIR}/c3/lib_linalg.h")
  message(SEND_ERROR "Can't find c3/lib_linalg.h in ${C3_INCLUDE_DIR}")
endif()
include_directories(${C3_INCLUDE_DIR})

set(CDYN_INCLUDE_DIR "" CACHE FILEPATH "Path to CDYN includes")

if(NOT EXISTS "${CDYN_INCLUDE_DIR}/cdyn/simulate.h")
  message(SEND_ERROR "Can't find cdyn/simulate.h in ${CDYN_INCLUDE_DIR}")
endif()
include_directories(${CDYN_INCLUDE_DIR})


# Set c3 library
set(C3_LIB_PATH "" CACHE FILEPATH "Path to Compressed Continuous Computation (C3) library")
find_library(c3lib c3 PATHS ${C3_LIB_PATH})
if(NOT c3lib)
  message(FATAL_ERROR "c3 library not found in ${C3_LIB_PATH}")
endif()

set(CDYN_LIB_PATH "" CACHE FILEPATH "Path to CDYN library")
find_library(cdynlib cdyn PATHS ${CDYN_LIB_PATH})
if(NOT cdynlib)
  message(FATAL_ERROR "cdyn library not found in ${CDYN_LIB_PATH}")
endif()


set(CMAKE_C_FLAGS "${compiler_flags}")

# SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)

add_subdirectory(src)
add_subdirectory(examples)
add_subdirectory(miscellaneous)
add_subdirectory(test)

# install(FILES cmake/c3sc-config.cmake DESTINATION ${main_lib_dest})
# install(EXPORT c3sc DESTINATION ${lib_dest})

