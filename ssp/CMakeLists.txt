#
# Copyright 2015 by Idiap Research Institute, http://www.idiap.ch
#
# See the file COPYING for the licence associated with this software.
#
# Author(s):
#   Phil Garner, December 2013
#

add_library(shared-lib SHARED
  ssp.cpp
  ar.cpp
  pitch.cpp
  arcodec.cpp
  cochlea.cpp
  )
set_target_properties(shared-lib
  PROPERTIES OUTPUT_NAME "ssp"
  PUBLIC_HEADER "ssp.h"
  )
target_link_libraries(shared-lib ${TARGET_LIBS})
list(APPEND LIBSSP_TARGETS shared-lib)

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
  PUBLIC_HEADER DESTINATION include
  )

# Slightly different location for the modules
install(
  TARGETS ${MODULE_TARGETS}
  LIBRARY DESTINATION lib/lube
  )
