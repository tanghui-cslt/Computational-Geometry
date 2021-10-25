#pragma once
#ifndef _LIBIGL_SHOW_LINE_H_
#define _LIBIGL_SHOW_LINE_H_

#include <igl/viewer/Viewer.h>
#include "openmesh_traits.h"
#include <queue>

typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits> MyMesh;
typedef std::queue<OpenMesh::FaceHandle>  q_faceHandle;

int libigl_lines(MyMesh *mesh,igl::viewer::Viewer &viewer);
int libigl_show_lines(MyMesh *mesh,igl::viewer::Viewer &viewer);
int libigl_temp_show_lines(MyMesh *mesh, igl::viewer::Viewer &viewer);
int libigl_touch_fire(MyMesh *mesh, igl::viewer::Viewer &viewer, q_faceHandle &Q_Fh, Eigen::MatrixXd &Mesh_Color);

#endif // ! _LIBIGL_SHOW_LINE_H_

