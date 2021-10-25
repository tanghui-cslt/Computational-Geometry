#pragma once
#ifndef _CALC_UV_H_
#define _CALC_UV_H_

#include "closed_loop.h"
#include <igl/viewer/Viewer.h>
#include <string>
void init(MyMesh *mesh, int id, igl::viewer::Viewer &viewer, std::string uv_sym);
void build_uv(MyMesh *mesh, int id, igl::viewer::Viewer &viewer);
#endif