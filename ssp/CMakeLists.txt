#
# Copyright 2015 by Idiap Research Institute, http://www.idiap.ch
#
# See the file COPYING for the licence associated with this software.
#
# Author(s):
#   Phil Garner, December 2013
#

set(HEADERS
  ssp.h
  cochlea.h
  warp.h
  filter.h
  )

add_library(ssp-shared SHARED
  ssp.cpp
  ar.cpp
  pitch.cpp
  arcodec.cpp
  cochlea.cpp
  filter.cpp
  )
set_target_properties(ssp-shared
  PROPERTIES OUTPUT_NAME "ssp"
  PUBLIC_HEADER "${HEADERS}"
  )
target_link_libraries(ssp-shared ${TARGET_LIBS})
list(APPEND LIBSSP_TARGETS ssp-shared)

# HTK file module
add_library(htk-lib MODULE htkfile.cpp)
target_link_libraries(htk-lib)
set_target_properties(htk-lib
  PROPERTIES OUTPUT_NAME "htk"
  )
list(APPEND MODULE_TARGETS htk-lib)

# Normal install path for the library
install(
  TARGETS ${LIBSSP_TARGETS}
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  PUBLIC_HEADER DESTINATION include/ssp
  )

# Slightly different location for the modules
install(
  TARGETS ${MODULE_TARGETS}
  LIBRARY DESTINATION lib/lube
  )
