#
# Copyright 2013 by Idiap Research Institute, http://www.idiap.ch
#
# See the file COPYING for the licence associated with this software.
#
# Author(s):
#   Phil Garner, December 2013
#

cmake_minimum_required(VERSION 2.8)

project(libssp)
set(VERSION 0.1)

set(CMAKE_CXX_FLAGS "-Wall -Werror -std=c++11")
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
include_directories(${CMAKE_SOURCE_DIR})

find_package(LibUBE REQUIRED)
include_directories(${LIBUBE_INCLUDE_DIR})
set(TARGET_LIBS ${LIBUBE_LIBRARIES})

add_subdirectory(ssp)

add_executable(waveplot waveplot.cpp)
target_link_libraries(waveplot shared-lib)

add_executable(specplot specplot.cpp)
target_link_libraries(specplot shared-lib)

add_executable(varcoder varcoder.cpp)
target_link_libraries(varcoder shared-lib)

enable_testing()
add_subdirectory(test)
