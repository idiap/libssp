#
# Copyright 2013 by Idiap Research Institute, http://www.idiap.ch
#
# See the file COPYING for the licence associated with this software.
#
# Author(s):
#   Phil Garner, December 2013
#

#
# Try to find Libube
# Once done this will define
#  LIBUBE_FOUND          - System has libube
#  LIBUBE_INCLUDE_DIR    - The libube include directory
#  LIBUBE_LIBRARIES      - The libraries needed to use Libube
#  LIBUBE_DEFINITIONS    - Compiler switches required for using libube
#  LIBUBE_VERSION_STRING - the version of libube found
#

find_package(PkgConfig)
pkg_check_modules(PC_LIBUBE QUIET libube)

set(LIBUBE_DEFINITIONS ${PC_LIBUBE_CFLAGS_OTHER})
set(LIBUBE_VERSION_STRING ${PC_LIBUBE_VERSION})

find_path(
  LIBUBE_INCLUDE_DIR lube.h
  HINTS ${PC_LIBUBE_INCLUDEDIR} ${PC_LIBUBE_INCLUDE_DIRS}
  PATH_SUFFIXES libube
)

find_library(
  LIBUBE_LIBRARIES ube
  HINTS ${PC_LIBUBE_LIBDIR} ${PC_LIBUBE_LIBRARY_DIRS}
)
list(APPEND LIBUBE_LIBRARIES ${CMAKE_DL_LIBS})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  Libube
  REQUIRED_VARS LIBUBE_LIBRARIES LIBUBE_INCLUDE_DIR
  VERSION_VAR LIBUBE_VERSION_STRING
)

mark_as_advanced(LIBUBE_INCLUDE_DIR LIBUBE_LIBRARIES)
