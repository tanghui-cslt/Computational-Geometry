#
# Try to find MOSEK
# Once done this will define
#
# MOSEK_FOUND           - system has MOSEK
# MOSEK_INCLUDE_DIRS    - the MOSEK include directories
# MOSEK_LIBRARIES       - Link these to use MOSEK
#

FIND_PATH(MOSEK_INCLUDE_DIR mosek.h
  PATHS 
  /usr/local/mosek/7/tools/platform/osx64x86/h/
  /usr/local/mosek/8/tools/platform/osx64x86/h/
    )

SET(SEARCH_PATHS "${MOSEK_INCLUDE_DIR}" "${MOSEK_INCLUDE_DIR}/../bin" "${MOSEK_INCLUDE_DIR}/lib")

set(MOSEK_LIBRARIES)
FIND_LIBRARY(MOSEK_LIBRARIES  NAMES mosek64 PATHS ${SEARCH_PATHS} NO_DEFAULT_PATH DPATH_SUFFIXES a lib dylib)

if(MOSEK_LIBRARIES AND MOSEK_INCLUDE_DIR)
message(STATUS "Found mosek: ${MOSEK_LIBRARIES}")
set(MOSEK_FOUND TRUE)
endif(MOSEK_LIBRARIES AND MOSEK_INCLUDE_DIR)

IF (MOSEK_FOUND)
   message(STATUS "Found MOSEK: ${MOSEK_INCLUDE_DIR}")
   SET(MOSEK_INCLUDE_DIRS ${MOSEK_INCLUDE_DIR} )
ELSE (MOSEK_FOUND)
    #add_definitions(-DIGL_NO_MOSEK)
    #message(WARNING "could NOT find MOSEK")
ENDIF (MOSEK_FOUND)
