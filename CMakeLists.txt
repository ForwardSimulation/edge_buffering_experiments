cmake_minimum_required(VERSION 2.8.12)
project(wfbuffered)

find_package(GSL REQUIRED)
message(STATUS "GSL headers in ${GSL_INCLUDE_DIRS}")
include_directories(BEFORE ${GSL_INCLUDE_DIRS})
include_directories(BEFORE ${wfbuffered_SOURCE_DIR}/tskit/c)
include_directories(BEFORE ${wfbuffered_SOURCE_DIR}/tskit/c/subprojects/kastore)

set(WFBUFFERED_SOURCES wfbuffered.cc
    rng.cc)


add_executable(wfbuffered ${WFBUFFERED_SOURCES})
target_link_libraries(wfbuffered PRIVATE GSL::gsl GSL::gslcblas)
