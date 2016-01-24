#
# Copyright 2015 by Idiap Research Institute, http://www.idiap.ch
#
# See the file COPYING for the licence associated with this software.
#
# Author(s):
#   Phil Garner, December 2013
#

# Testing
#
# The tests get built here, but they are run from a different cmake
# file.  This lets us compare output with a reference output rather
# than just test that something runs.

set(TEST_DIR ${CMAKE_SOURCE_DIR}/test)
add_definitions(-DTEST_DIR="${TEST_DIR}")

add_executable(test-ssp test-ssp.cpp)
target_link_libraries(test-ssp ssp-shared)
add_test(
  NAME ssp
  COMMAND ${CMAKE_COMMAND} -D TEST_DIR=${TEST_DIR} -P ${CMAKE_CURRENT_SOURCE_DIR}/test-ssp.cmake
  )

add_executable(test-cochlea test-cochlea.cpp)
target_link_libraries(test-cochlea ssp-shared)
add_test(
  NAME cochlea
  COMMAND ${CMAKE_COMMAND} -D TEST_DIR=${TEST_DIR} -P ${CMAKE_CURRENT_SOURCE_DIR}/test-cochlea.cmake
  )
