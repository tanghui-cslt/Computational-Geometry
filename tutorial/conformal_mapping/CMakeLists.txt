cmake_minimum_required(VERSION 2.8.12)

project(Conformal_Mapping)

set(CMAKE_CXX_STANDARD 11)

add_executable(1-Conformal_Mapping main.cpp)

target_include_directories(1-Conformal_Mapping PRIVATE ${LIBIGL_INCLUDE_DIRS})

target_compile_definitions(1-Conformal_Mapping PRIVATE ${LIBIGL_DEFINITIONS})

target_link_libraries(1-Conformal_Mapping ${LIBIGL_LIBRARIES} ${LIBIGL_VIEWER_EXTRA_LIBRARIES} ${LIBIGL_OPENGL_EXTRA_LIBTATIES} ${LIBIGL_OPENGL_GLFW_EXTRA_LIBARARIES} )
