#pragma once
#ifndef _CLOSED_LOOP_H_
#define _CLOSED_LOOP_H_


#include "openmesh_traits.h"
#include <vector>
#include <iostream>

typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits> MyMesh;

bool find_loop(MyMesh *mesh,int &id);
void spanning_tree(MyMesh *mesh);
void set_mesh_true(MyMesh *mesh);
void set_loop_flag(MyMesh *mesh, int &id);		//标记当前路径
void save_false_flag(MyMesh *mesh);
std::vector< std::vector<EdgeHandle> > get_loop();

#endif