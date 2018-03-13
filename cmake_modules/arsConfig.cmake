# - Try to find Library ars
# Once done, this will define
#
#  ars_FOUND - system has ars module
#  ars_INCLUDE_DIRS - the ars include directories
#  ars_LIBRARY_DIRS - the ars library directories
#  ars_LIBRARIES - link these to use ars


# Uses  directory to search mrf_segmentation directory!
set(ars_PREFIX_DIR /usr/local)
message(STATUS "Searching ars in directory ${ars_PREFIX_DIR}." )

# Searches include directory /usr/local/include/ars
find_path(ars_INCLUDE_DIR ars ${ars_PREFIX_DIR}/include)
message(STATUS "    ars_INCLUDE_DIR ${ars_INCLUDE_DIR}." )
set(ars_INCLUDE_DIRS ${ars_INCLUDE_DIR})
  
# Searches library librimagraph.a in /usr/local/lib
find_path(ars_LIBRARY_DIR librimagraph.a ${ars_PREFIX_DIR}/lib)
message(STATUS "    ars_LIBRARY_DIR ${ars_LIBRARY_DIR}." )
set(ars_LIBRARY_DIRS ${ars_PREFIX_DIR}/lib)

# Sets the names of library components (actually A name and A component)
find_library(ars_LIBRARY ars ${ars_LIBRARY_DIRS})
message(STATUS "    ars_LIBRARY ${ars_LIBRARY}." )
set(ars_LIBRARIES ${ars_LIBRARY})

if(("${ars_INCLUDE_DIR}" STREQUAL "ars_INCLUDE_DIR-NOTFOUND") OR
   ("${ars_LIBRARY_DIRS}" STREQUAL "ars_LIBRARY_DIRS-NOTFOUND") OR
   ("${ars_LIBRARY}" STREQUAL "ars_LIBRARY-NOTFOUND")
  )
  message(STATUS "Library ars NOT found")
  unset(ars_FOUND)
  unset(ars_INCLUDE_DIR)
  unset(ars_LIBRARY_DIR)
  unset(ars_LIBRARY)
  unset(ars_LIBRARIES)
endif()

mark_as_advanced(ars_INCLUDE_DIRS ars_LIBRARY_DIRS ars_LIBRARIES)

