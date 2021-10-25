#pragma once
#ifndef _CONVERT_TO_LIBIGL_H_
#define _CONVERT_TO_LIBIGL_H_

#include <igl/viewer/Viewer.h>
#include "openmesh_traits.h"
typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits> MyMesh;

void convert_libigl(MyMesh &t, igl::viewer::Viewer &viewer);
void init(MyMesh &t);
bool pre_draw(igl::viewer::Viewer &viewer);
bool init_new_window(igl::viewer::Viewer &viewer);
#endif // !_CONVERT_TO_LIBIGL_H_

