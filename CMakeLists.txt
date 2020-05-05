cmake_minimum_required(VERSION 2.8.12)
project(wfbuffered)

find_package(GSL REQUIRED)
message(STATUS "GSL headers in ${GSL_INCLUDE_DIRS}")
include_directories(BEFORE ${GSL_INCLUDE_DIRS})
include_directories(BEFORE ${wfbuffered_SOURCE_DIR}/tskit/c)
include_directories(BEFORE ${wfbuffered_SOURCE_DIR}/tskit/c/subprojects/kastore)

set(WFBUFFERED_SOURCES wfbuffered.cc
    rng.cc
    tskit_tools.cc)

file(GLOB TSKIT_SOURCES ${wfbuffered_SOURCE_DIR}/tskit/c/tskit/*.c)
file(GLOB KASTORE_SOURCES ${wfbuffered_SOURCE_DIR}/tskit/c/subprojects/kastore/kastore.c)

set(ALL_SOURCES ${WFBUFFERED_SOURCES}
    ${TSKIT_SOURCES}
    ${KASTORE_SOURCES})


add_executable(wfbuffered ${ALL_SOURCES})
target_link_libraries(wfbuffered PRIVATE GSL::gsl GSL::gslcblas)
