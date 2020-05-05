cmake_minimum_required(VERSION 2.8.12)
project(wfbuffered)

include_directories(BEFORE ${wfbuffered_SOURCE_DIR}/tskit/c)
include_directories(BEFORE ${wfbuffered_SOURCE_DIR}/tskit/c/subprojects/kastore)

set(WFBUFFERED_SOURCES wfbuffered.cc)

add_executable(wfbuffered ${WFBUFFERED_SOURCES})

