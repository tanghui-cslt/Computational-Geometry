cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)
project(test-open)
include_directories (
  OpenMesh-8.0
  ${CMAKE_CURRENT_SOURCE_DIR}
)
if (WIN32 AND  BUILD_APPS )
    # prepare bundle generation cmake file and add a build target for it
    configure_file ("${CMAKE_SOURCE_DIR}/cmake/fixbundle.cmake.win.in"
       "${CMAKE_BINARY_DIR}/fixbundle.win.cmake" @ONLY IMMEDIATE)

    if ( NOT "${CMAKE_GENERATOR}" MATCHES "MinGW Makefiles" )
      # let bundle generation depend on all targets
      add_dependencies (fixbundle QtViewer DecimaterGui)
    endif()
 
 endif()
  
add_subdirectory(OpenMesh-8.0 output)
acg_append_files (headers "*.hh" .)
acg_append_files (sources "*.cc" .)
#acg_add_executable (${targetName} ${headers} ${sources})
#acg_add_executable (${targetName} ${headers} ${sources})
# target_link_libraries (${targetName}
  # OpenMeshCore
  # OpenMeshTools
# )

set(TRIANGLE_DIR "triangle")
  add_subdirectory("${TRIANGLE_DIR}" "triangle")
  list(APPEND LIBIGL_TRIANGLE_EXTRA_LIBRARIES "triangle")
  # list(APPEND LIBIGL_INCLUDE_DIRS ${TRIANGLE_DIR})
  # list(APPEND LIBIGL_TRIANGLE_EXTRA_LIBRARIES "triangle")
  # list(APPEND LIBIGL_EXTRA_LIBRARIES ${LIBIGL_TRIANGLE_EXTRA_LIBRARIES})

  # if(LIBIGL_USE_STATIC_LIBRARY)
    # CompileIGL_Module("triangle" "")
    # target_include_directories(igl_triangle PRIVATE ${TRIANGLE_DIR})
  # endif()
message(STATUS "${LIBIGL_TRIANGLE_EXTRA_LIBRARIES}")