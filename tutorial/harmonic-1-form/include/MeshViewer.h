#pragma once
#ifndef _Mesh_VIEWER_
#define _Mesh_VIEWER_

#include <igl/viewer/Viewer.h>
#include <igl/viewer/Viewer.h>
#include <nanogui/formhelper.h>
#include <nanogui/screen.h>
#include <iostream>
#include <queue>

#include "LibiglShowLines.h"
#include "FindClosedLoop.h"
#include "HarmonicForm.h"


template <typename TMesh>  void convert_libigl(TMesh &t, igl::viewer::Viewer &viewer);
template <typename TMesh>  void init(TMesh &t);
bool pre_draw(igl::viewer::Viewer &viewer);
bool init_new_window(igl::viewer::Viewer &viewer);

using namespace std;
typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits> MyMesh;
Eigen::MatrixXd V;
Eigen::MatrixXi F;
HarmonicForm<MyMesh> h_form;
Eigen::MatrixXd Mesh_Color;

MyMesh *mesh = NULL;

bool touch_flag = false;
bool find_loop_flag = false;
bool puring_flag = false;
bool show_loop_flag = false;
bool calc_uv_flag = false;
Eigen::MatrixXd face_color;
bool open_close_flag = false;	//等待循环完成
int id = 0;						//寻找过程中的路径的变化
int loop_id = 0;				//输入显示第loop_id个路径
queue<OpenMesh::FaceHandle> Q_Fh;


template <typename TMesh>  void convert_libigl(TMesh &t, igl::viewer::Viewer &viewer)
{

	size_t v_row = t.n_vertices();
	size_t f_row = t.n_faces();

	V.resize(v_row, 3);
	F.resize(f_row, 3);
	Mesh_Color.resize(f_row, 3);
	face_color.resize(f_row, 3);
	int count = 0;

	for (auto it = t.vertices_begin(); it != t.vertices_end(); it++)
	{

		auto point = t.point(*it);
		V(count, 0) = point.data()[0];
		V(count, 1) = point.data()[1];
		V(count, 2) = point.data()[2];
		count++;
	}

	count = 0;
	for (auto it = t.faces_begin(); it != t.faces_end(); it++)
	{
		int col = 0;

		for (auto i = t.fv_begin(*it); i != t.fv_end(*it); i++)
		{
			auto vertex = *i;
			F(count, col) = vertex.idx();
			col++;
		}
		face_color.row(count) << 1, 1, 1;
		count++;
	}
	viewer.data.set_mesh(V, F);

}
template <typename TMesh>  void init(TMesh &t)
{
	OpenMesh::FaceHandle Init_Fh = *(t.faces_begin());
	t.data(Init_Fh).set_labeled(true);

	Q_Fh.push(Init_Fh);
	mesh = &t;
	Mesh_Color.row(Init_Fh.idx()) << 0, 1, 0;
}


bool init_new_window(igl::viewer::Viewer &viewer)
{

	viewer.ngui->addWindow(Eigen::Vector2i(300, 10), "Harmanic 1 form");
	viewer.ngui->addButton("calc homology group's basis ", [&]() {
		libigl_touch_fire<MyMesh>(mesh, viewer, Q_Fh, Mesh_Color);
		int value = libigl_show_lines<MyMesh>(mesh, viewer);
		while (value != 0)
		{
			//cout << "puring ok\n";
			value = libigl_show_lines<MyMesh>(mesh, viewer);
		}


		puring_flag = false;
		spanning_tree<MyMesh>(mesh);
		open_close_flag = find_loop<MyMesh>(mesh, id);		//找第id个loop，并判断是否找完

		find_loop_flag = true;

		//puring_flag = true; 
	});


	viewer.ngui->addButton("calc harmonic 1 form", [&]()
	{
		calc_uv_flag = true;
	}
	);

	viewer.ngui->addVariable<int >("show generator:", [&](int val) {
		loop_id = val; // set
		set_loop_flag<MyMesh>(mesh, loop_id);		//标记当前路径
		show_loop_flag = true;
		viewer.ngui->refresh();
	}, [&]() {
		return loop_id; // get
	});

	viewer.ngui->addButton("show u", [&]()
	{
		viewer.data.clear();
		viewer.data.set_mesh(V, F);
		viewer.data.set_colors(face_color);
		h_form.init(mesh, loop_id, viewer, "u");


	});
	viewer.ngui->addButton("show v", [&]()
	{
		viewer.data.clear();
		viewer.data.set_mesh(V, F);
		viewer.data.set_colors(face_color);
		h_form.init(mesh, loop_id, viewer, "v");


	});

	viewer.ngui->addButton("show uv", [&]()
	{
		viewer.data.clear();
		viewer.data.set_mesh(V, F);
		viewer.data.set_colors(face_color);
		h_form.init(mesh, loop_id, viewer, "uv");


	});

	viewer.screen->performLayout();

	return false;
}


bool pre_draw(igl::viewer::Viewer &viewer)
{
	if (viewer.core.is_animating)
	{

		if (puring_flag)		//提纯
		{

			if (libigl_show_lines<MyMesh>(mesh, viewer) == 0)
			{
				//cout << "puring ok\n";
				puring_flag = false;
				spanning_tree<MyMesh>(mesh);
				open_close_flag = find_loop<MyMesh>(mesh, id);		//找第id个loop，并判断是否找完

				find_loop_flag = true;
			}
		}
		if (find_loop_flag)		//寻找路径
		{
			do {
				if (libigl_lines<MyMesh>(mesh, viewer) == 0)
				{
					save_false_flag<MyMesh>(mesh);					//保存找到的loop

					if (open_close_flag)
					{
						find_loop_flag = false;
						cout << "find all loops!!\n";
					}
					else {
						open_close_flag = find_loop<MyMesh>(mesh, id);
					}

				}

			} while (find_loop_flag);
		}
		if (calc_uv_flag)
		{

			h_form.build_uv(mesh, viewer, 0);

			calc_uv_flag = false;
		}
		if (show_loop_flag)		//显示路径
		{
			if (libigl_temp_show_lines<MyMesh>(mesh, viewer) == 0)
			{
				show_loop_flag = false;
			}
		}
	}
	viewer.ngui->refresh();
	return false;
}

#endif // !_CONVERT_TO_LIBIGL_H_

