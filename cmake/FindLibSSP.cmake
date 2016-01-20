#
# Copyright 2016 by Idiap Research Institute, http://www.idiap.ch
#
# See the file COPYING for the licence associated with this software.
#
# Author(s):
#   Phil Garner, January 2016
#

find_package(PkgConfig)
pkg_check_modules(PC_LIBSSP QUIET libssp)

set(LIBSSP_DEFINITIONS ${PC_LIBSSP_CFLAGS_OTHER})

find_path(
  LIBSSP_INCLUDE_DIR ssp.h
  HINTS ${PC_LIBSSP_INCLUDEDIR} ${PC_LIBSSP_INCLUDE_DIRS}
  PATH_SUFFIXES libssp
)

find_library(
  LIBSSP_LIBRARY ssp
  HINTS ${PC_LIBSSP_LIBDIR} ${PC_LIBSSP_LIBRARY_DIRS}
)

set(LIBSSP_LIBRARIES ${LIBSSP_LIBRARY})
set(LIBSSP_INCLUDE_DIRS ${LIBSSP_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  LibSSP DEFAULT_MSG
  LIBSSP_LIBRARY LIBSSP_INCLUDE_DIR
)

mark_as_advanced(LIBSSP_INCLUDE_DIR LIBSSP_LIBRARY)
