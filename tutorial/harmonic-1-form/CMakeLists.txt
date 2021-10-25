cmake_minimum_required(VERSION 2.8.12)

project(Harmonic-1-form)
message(STATUS "CMAKE_C_COMPILER: ${CMAKE_C_COMPILER}")
include_directories("./include")
AUX_SOURCE_DIRECTORY(src DIR_SRCS)

set(CMAKE_CXX_STANDARD 11)
add_executable( 2-Harmonic-1-form ${DIR_SRCS})
target_include_directories(2-Harmonic-1-form PRIVATE ${LIBIGL_INCLUDE_DIRS})
target_compile_definitions(2-Harmonic-1-form PRIVATE ${LIBIGL_DEFINITIONS})
target_link_libraries(2-Harmonic-1-form ${LIBIGL_LIBRARIES} ${LIBIGL_VIEWER_EXTRA_LIBRARIES} ${LIBIGL_OPENGL_EXTRA_LIBTATIES} ${LIBIGL_OPENGL_GLFW_EXTRA_LIBARARIES} )