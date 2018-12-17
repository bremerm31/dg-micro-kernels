cmake_minimum_required (VERSION 3.0 FATAL_ERROR)
project (DGBench CXX C)

find_program(CMAKE_CXX_COMPILER NAMES $ENV{CXX} g++ PATHS ENV PATH NO_DEFAULT_PATH)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # require at least intel 17.0
  if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 17.0)
    message(FATAL_ERROR "Intel version must be at least 17.0")
  endif()
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # require at least 6.0
  if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 6.0)
    message(FATAL_ERROR "GCC version must be at least 6.0")
  endif()
else()
  message(WARNING "Using untested compiler")
endif()

add_definitions(-std=c++14 -g -ftemplate-backtrace-limit=0 -march=native)

add_definitions(-Wall)

set(CMAKE_VERBOSE_MAKEFILE ON)

get_filename_component (default_prefix "../install" ABSOLUTE)
set (CMAKE_INSTALL_PREFIX ${default_prefix} CACHE STRING
      "Choose the installation directory; by default it installs in install."
FORCE)

find_package(blaze REQUIRED)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
link_libraries(blaze::blaze ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

find_package(benchmark REQUIRED)

add_executable(
  dg-micro-kernels
  src/main.cpp
)

target_include_directories(dg-micro-kernels PRIVATE ${benchmark_INCLUDE_DIRS})
target_link_libraries(dg-micro-kernels ${benchmark_LIBRARIES})