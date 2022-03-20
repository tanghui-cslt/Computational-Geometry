#pragma once
#ifndef _LIBIGL_SHOW_LINE_H_
#define _LIBIGL_SHOW_LINE_H_

#include <igl/viewer/Viewer.h>
#include "openmesh_traits.h"
#include <queue>


typedef std::queue<OpenMesh::FaceHandle>  q_faceHandle;

template <typename TMesh> int libigl_lines(TMesh *mesh,igl::viewer::Viewer &viewer);
template <typename TMesh> int libigl_show_lines(TMesh *mesh,igl::viewer::Viewer &viewer);
template <typename TMesh> int libigl_temp_show_lines(TMesh *mesh, igl::viewer::Viewer &viewer);
template <typename TMesh> int libigl_touch_fire(TMesh *mesh, igl::viewer::Viewer &viewer, q_faceHandle &Q_Fh, Eigen::MatrixXd &Mesh_Color);




template <typename TMesh>
int libigl_touch_fire(TMesh *mesh, igl::viewer::Viewer &viewer, q_faceHandle &Q_Fh, Eigen::MatrixXd &Mesh_Color)
{
	int touch_flag = 0;
	while (!Q_Fh.empty())
	{
		OpenMesh::FaceHandle Init_Fh = Q_Fh.front();
		Q_Fh.pop();
		//cout << "orig = " << Init_Fh.idx() << endl;
		for (auto it = mesh->fh_begin(Init_Fh); it != mesh->fh_end(Init_Fh); it++)
		{
			//auto he = it.handle();
			HalfedgeHandle hhf = mesh->opposite_halfedge_handle(*it);
			FaceHandle fh = mesh->face_handle(hhf);
			EdgeHandle eh = mesh->edge_handle(*it);

			if (mesh->data(eh).is_labeled()) continue;
			else
			{
				if (!(mesh->data(fh).is_labeled()))
				{
					mesh->data(fh).set_labeled(true);
					Q_Fh.push(fh);
					mesh->data(eh).set_labeled(true);
					Mesh_Color.row(fh.idx()) << 1, 0, 0;
				}
			}
		}
		//viewer.data.set_colors(Mesh_Color);
	}
	return touch_flag;
}

template <typename TMesh> bool isbranch(TMesh *mesh, EdgeHandle eh, int * num) {
	if (mesh->data(eh).is_labeled() == 1)	return false;
	VertexHandle v1 = mesh->to_vertex_handle(mesh->halfedge_handle(eh, 0));
	VertexHandle v2 = mesh->from_vertex_handle(mesh->halfedge_handle(eh, 0));
	int a = v1.idx();
	int b = v2.idx();

	return (num[a] == 1 || num[b] == 1);
}

template <typename TMesh>
int libigl_show_lines(TMesh *mesh, igl::viewer::Viewer &viewer)
{
	int branch_flag = 0;
	Eigen::MatrixXd P1, P2;
	Eigen::MatrixXd V = viewer.data.V;
	Eigen::MatrixXi F = viewer.data.F;
	do
	{
		branch_flag = 0;
		

		int  num_label_edge = 0;
		int * num = NULL;
		num = new int[mesh->n_vertices()];
		memset(num, 0, sizeof(int)*mesh->n_vertices());
	
		for (auto ie = mesh->edges_begin(); ie != mesh->edges_end(); ie++)
		{

			if (mesh->data(*ie).is_labeled())
			{
				continue;
			}
			VertexHandle v1 = mesh->to_vertex_handle(mesh->halfedge_handle(*ie, 0));
			VertexHandle v2 = mesh->from_vertex_handle(mesh->halfedge_handle(*ie, 0));
			num[v1.idx()] ++;
			num[v2.idx()] ++;

			num_label_edge++;
		}

		P1.resize(num_label_edge, 3);
		P2.resize(num_label_edge, 3);
		int num_times = 0;

		for (auto ie = mesh->edges_begin(); ie != mesh->edges_end(); ie++)
		{
			if (mesh->data(*ie).is_labeled())
				continue;
			if (isbranch(mesh, *ie, num))
			{
				//Q_eh.push(*ie);
				mesh->data(*ie).set_labeled(true);
				branch_flag++;
			}

			VertexHandle v1 = mesh->to_vertex_handle(mesh->halfedge_handle(*ie, 0));
			VertexHandle v2 = mesh->from_vertex_handle(mesh->halfedge_handle(*ie, 0));

			P1.row(num_times) = V.row(v1.idx());
			P2.row(num_times) = V.row(v2.idx());
			num_times++;
		}
	} while (branch_flag != 0);
	viewer.data.clear();
	viewer.data.set_mesh(V, F);
	viewer.data.add_edges(P1, P2, Eigen::RowVector3d(1, 0, 0));
	return branch_flag;
}

template <typename TMesh>
int libigl_lines(TMesh *mesh, igl::viewer::Viewer &viewer)
{
	int branch_flag = 0;

	Eigen::MatrixXd V = viewer.data.V;
	Eigen::MatrixXi F = viewer.data.F;

	int  num_label_edge = 0;
	int * num = NULL;
	num = new int[mesh->n_vertices()];
	memset(num, 0, sizeof(int)*mesh->n_vertices());
	for (auto ie = mesh->edges_begin(); ie != mesh->edges_end(); ie++)
	{

		if (mesh->data(*ie).is_labeled())
		{
			continue;
		}
		VertexHandle v1 = mesh->to_vertex_handle(mesh->halfedge_handle(*ie, 0));
		VertexHandle v2 = mesh->from_vertex_handle(mesh->halfedge_handle(*ie, 0));
		num[v1.idx()] ++;
		num[v2.idx()] ++;

		num_label_edge++;
	}

	
	for (auto ie = mesh->edges_begin(); ie != mesh->edges_end(); ie++)
	{
		if (mesh->data(*ie).is_labeled())
			continue;
		if (isbranch(mesh, *ie,num))
		{
			//			Q_eh.push(*ie);
			mesh->data(*ie).set_labeled(true);
			branch_flag++;
		}

	}

	return branch_flag;
}

template <typename TMesh>
int libigl_temp_show_lines(TMesh *mesh, igl::viewer::Viewer &viewer)
{
	int branch_flag = 0;

	Eigen::MatrixXd V = viewer.data.V;
	Eigen::MatrixXi F = viewer.data.F;

	int  num_label_edge = 0;
	
	for (auto ie = mesh->edges_begin(); ie != mesh->edges_end(); ie++)
	{
		if (mesh->data(*ie).is_labeled())
		{
			continue;
		}
		num_label_edge++;
	}

	Eigen::MatrixXd P1(num_label_edge, 3);
	Eigen::MatrixXd P2(num_label_edge, 3);
	//cout << "loop_len=" << num_label_edge << endl;
	int num_times = 0;
	for (auto ie = mesh->edges_begin(); ie != mesh->edges_end(); ie++)
	{
		if (mesh->data(*ie).is_labeled())
			continue;

		VertexHandle v1 = mesh->to_vertex_handle(mesh->halfedge_handle(*ie, 0));
		VertexHandle v2 = mesh->from_vertex_handle(mesh->halfedge_handle(*ie, 0));

		//cout << "P1" << endl;
		P1.row(num_times) = V.row(v1.idx());
		P2.row(num_times) = V.row(v2.idx());
		num_times++;
	}

	viewer.data.clear();
	viewer.core.show_texture = false;
	viewer.data.set_mesh(V, F);

	viewer.data.add_edges(P1, P2, Eigen::RowVector3d(1, 0, 0));
	return branch_flag;
}
#endif // ! _LIBIGL_SHOW_LINE_H_

