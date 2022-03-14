#ifndef _HARMONIC_FORM_
#define _HARMONIC_FORM_

#include <igl/viewer/Viewer.h>
#include <igl/cotmatrix.h>
#include <igl/PNG/readPNG.h>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/Cholesky>
#include <queue>
#include <vector>

#include <iostream>
#include <igl/per_face_normals.h>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include "FindClosedLoop.h"

using namespace std;
using namespace Eigen;


/*! \brief HarmonicForm class
*
*	class for calculating harmonic form
*
*/
template<class TMesh> class HarmonicForm 
{
private:
	
	vector<vector<EdgeHandle> > vector_all_loop;
	vector<EdgeHandle> vector_loop;
	MatrixXd V;
	MatrixXi F;
	
	int v_size;													//number of v
	int e_size;													//number of e
	int f_size;
	int loop_size = 0;											//loop's length 

	bool curr_label = false;									//label of each edge
	SparseMatrix<double> d0;									//d0
	SparseMatrix<double, Eigen::RowMajor> star1;				//*1
	vector<VertexHandle> vector_labled_boundary;				//当前选中的loop
	SparseMatrix<double> L;										//L = cotangent
	VectorXd characteristic_w;									//characteristic 1 form
	MatrixXd harmonic_w;										//harmonic 1 form

	vector< vector<double> >all_edge_v;							//dual harmonic 1 form
	vector< vector<double> > all_ver_u;
	vector< vector<double> > all_ver_v;
	
	VectorXd area;												//area
	MatrixXd N_faces;											//normal vector

	
public:
	bool _is_boundry_vertex(VertexHandle vh, TMesh *mesh, vector<EdgeHandle> loop);
	bool _is_boundry_edge(EdgeHandle eh);
	void _set_clear(TMesh *mesh);
	void _init_all_data(TMesh *mesh, igl::viewer::Viewer &viewer, int id);
	SparseMatrix<double>  _calc_laplacian_matrix(TMesh *mesh);
	void _find_boundary(TMesh *mesh,  int id);
	VectorXd _constrained_laplacian_matrix(TMesh *mesh, int id);
	VectorXd _solve_0_form(SparseMatrix<double, Eigen::RowMajor> charac_matrix, VectorXd b);
	VectorXd _solve_alpha(TMesh *mesh);
	void _calc_uv_on_vertex(TMesh *mesh, int id, string uv_sym);
	void _calc_closed_1_form(TMesh *mesh,  VectorXd xy);
	void _calc_harmonic_1_form(TMesh *mesh, int id);
	void _calc_u_on_vertex(TMesh *mesh, int id);
	VectorXd _calc_area(TMesh *mesh);
	MatrixXd _calc_vector_2_form(TMesh *mesh, Eigen::MatrixXd N, Eigen::VectorXd area, int id);
	void _calc_v_on_vertex( TMesh *mesh, int id);
	void _calc_dual_harmonic_1_form(TMesh *mesh, Eigen::MatrixXd star_w);
	void _save_uv(TMesh *mesh,  string uv);
	void build_uv(TMesh *mesh, igl::viewer::Viewer &viewer, int id);
	void _set_id_uv(int id, string uv_sym, MatrixXd &uv);
	void init(TMesh *mesh, int id, igl::viewer::Viewer &viewer, string uv_sym);
	
	
};



template<class TMesh>
bool HarmonicForm<TMesh>::_is_boundry_vertex(VertexHandle vh, TMesh *mesh, vector<EdgeHandle> loop)
{
	bool flag = false;
	size_t size = loop.size();
	for (size_t i = 0; i < size; i++)
	{
		auto half_0 = mesh->halfedge_handle(loop[i], 0);
		auto v0 = mesh->from_vertex_handle(half_0);
		auto v1 = mesh->to_vertex_handle(half_0);

		if (v1 == vh ||
			v0 == vh)
		{
			flag = true;
			return flag;
		}
	}
	return flag;

}

template<class TMesh>
bool HarmonicForm<TMesh>::_is_boundry_edge(EdgeHandle eh)
{
	for (auto ei = vector_loop.begin(); ei != vector_loop.end(); ei++)
		if (*ei == eh)		return true;
	return false;
}

template<class TMesh>
void HarmonicForm<TMesh>::_set_clear(TMesh *mesh)
{
	for (auto ve = vector_labled_boundary.begin(); ve != vector_labled_boundary.end(); ve++)
	{
		mesh->data(*ve).set_u(-1);
		mesh->data(*ve).set_v(-1);
		mesh->data(*ve).set_boundary(-1);
		mesh->data(*ve).set_boundary_clear();
	}
}

template<class TMesh>
void HarmonicForm<TMesh>::_init_all_data(TMesh *mesh, igl::viewer::Viewer &viewer, int id)
{

	v_size = (int)mesh->n_vertices();
	e_size = (int)mesh->n_edges();
	f_size = (int)mesh->n_faces();


	vector_all_loop = get_loop<TMesh>();

	loop_size = (int)vector_all_loop.size();

	d0.resize(e_size, v_size);
	star1.resize(e_size, e_size);


	characteristic_w.resize(e_size);
	harmonic_w.resize(loop_size, e_size);



	N_faces.resize(f_size, 3);



	V = viewer.data.V;
	F = viewer.data.F;
	igl::cotmatrix(V, F, L);
	igl::per_face_normals(viewer.data.V, viewer.data.F, N_faces);
}

template<class TMesh>
SparseMatrix<double>  HarmonicForm<TMesh>::_calc_laplacian_matrix(TMesh *mesh)
{
	SparseMatrix<double> lapacian_matrix(v_size, v_size);

	//build d0, V->E
	for (auto eh = mesh->edges_begin(); eh != mesh->edges_end(); eh++)
	{
		auto half_handle = mesh->halfedge_handle(*eh, 0);
		auto v0 = mesh->to_vertex_handle(half_handle);
		auto v1 = mesh->from_vertex_handle(half_handle);
		int e_id = eh->idx();
		int v0_id = v0.idx();
		int v1_id = v1.idx();
		d0.insert(e_id, v0_id) = 1;
		d0.insert(e_id, v1_id) = -1;
	}
	// build star, E->E*
	for (auto eh = mesh->edges_begin(); eh != mesh->edges_end(); eh++)
	{
		auto half_handle = mesh->halfedge_handle(*eh, 0);
		auto v0 = mesh->to_vertex_handle(half_handle);
		auto v1 = mesh->from_vertex_handle(half_handle);
		int e_id = eh->idx();
		int v0_id = v0.idx();
		int v1_id = v1.idx();
		star1.insert(e_id, e_id) = (L.coeff(v0_id, v1_id)) / 2;
	}
	Eigen::SparseMatrix<double> d0_transpose = d0.transpose();
	//lapacian = d0^T * star1 * d0
	lapacian_matrix = d0_transpose * star1*d0;
	return lapacian_matrix;
}

//for finding the boundaries of closed surface.
template<class TMesh>
void HarmonicForm<TMesh>::_find_boundary(TMesh *mesh, int id)
{

	id = id % loop_size;

	vector_loop = vector_all_loop[id];

	int len = (int)vector_loop.size();


	_set_clear(mesh);
	vector<VertexHandle> vector_boundary;

	//select one edge
	//get halfedge which is No. 0 of the edge

	EdgeHandle vh = vector_loop[0];							
	auto half_0 = mesh->halfedge_handle(vh, 0);				
	auto temp_half_0 = half_0;								
	//get halfedge which is No. 1 of the edge
	auto half_1 = mesh->halfedge_handle(vh, 1);				
	auto temp_half_1 = half_1;								
															
	int forward_id = 0;
	int backward_id = 0;
	for (int i = 0; i < len; i++)
	{
		forward_id = (forward_id + 1) % len;
		backward_id = (backward_id - 1 + len) % len;
		EdgeHandle pre_eh = vector_loop[backward_id];
		EdgeHandle next_eh = vector_loop[forward_id];

		auto ori_vh = mesh->to_vertex_handle(temp_half_0);
		int ori_id = ori_vh.idx();

		temp_half_0 = mesh->next_halfedge_handle(temp_half_0);		//这个半边的下一个半边

		// find the 1-ring of ori_id in ccw rotation 
		while (mesh->edge_handle(temp_half_0) != pre_eh &&			//判断是否找到边界
			mesh->edge_handle(temp_half_0) != next_eh)
		{
			auto v0 = mesh->to_vertex_handle(temp_half_0);				//保存边界


			if (!_is_boundry_vertex(v0, mesh, vector_loop))
			{
				mesh->data(v0).set_boundary(0, ori_id);
				mesh->data(v0).set_boundary(0);
				vector_boundary.push_back(v0);
			}

			temp_half_0 = mesh->opposite_halfedge_handle(temp_half_0);
			temp_half_0 = mesh->next_halfedge_handle(temp_half_0);
		}

		//1 号半边

		ori_vh = mesh->to_vertex_handle(temp_half_1);
		ori_id = ori_vh.idx();

		temp_half_1 = mesh->next_halfedge_handle(temp_half_1);		//这个半边的下一个半边
		VertexHandle v1 = mesh->to_vertex_handle(temp_half_1);
		// find the 1-ring of ori_id in cw rotation 
		while (mesh->edge_handle(temp_half_1) != pre_eh &&			//判断是否找到边界
			mesh->edge_handle(temp_half_1) != next_eh)
		{
			v1 = mesh->to_vertex_handle(temp_half_1);				//保存边界


			if (!_is_boundry_vertex(v1, mesh, vector_loop))

			{
				mesh->data(v1).set_boundary(1, ori_id);

				mesh->data(v1).set_boundary(1);
				vector_boundary.push_back(v1);
			}

			temp_half_1 = mesh->opposite_halfedge_handle(temp_half_1);
			temp_half_1 = mesh->next_halfedge_handle(temp_half_1);
		}
	}

	vector_labled_boundary.assign(vector_boundary.begin(), vector_boundary.end());
}

//the left boundary vertices are 0, the right boundary vertices are 1;
template<class TMesh>
VectorXd HarmonicForm<TMesh>::_constrained_laplacian_matrix(TMesh *mesh, int id)
{
	id = id % loop_size;

	vector_loop = vector_all_loop[id];

	int len = (int)vector_loop.size();
	VectorXd b(v_size - len);
	SparseMatrix<double, Eigen::RowMajor> charac_matrix(v_size - len, v_size - len);
	int count = 0;	//new id
	SparseMatrix<double> closed_matrix = _calc_laplacian_matrix(mesh);
	for (auto iv = mesh->vertices_begin(); iv != mesh->vertices_end(); ++iv)
	{
		//new id
		if (_is_boundry_vertex(*iv, mesh, vector_loop))
			mesh->data(*iv).set_new_id(0);

		else
		{
			mesh->data(*iv).set_new_id(count);
			count++;
		}
	}
	//calc characteristic 1 form, and b 
	for (auto iv = mesh->vertices_begin(); iv != mesh->vertices_end(); ++iv)
	{
		if (_is_boundry_vertex(*iv, mesh, vector_loop))	continue;
		int new_id = mesh->data(*iv).get_new_id();
		int id = iv->idx();
		double temp_a = 0;
		for (auto ve = mesh->vv_begin(*iv); ve != mesh->vv_end(*iv); ++ve)//l-ring 之和为0
		{
			auto ano_vh = ve;

			int ano_vid = ano_vh->idx();
			temp_a = closed_matrix.coeff(iv->idx(), ano_vid);
			int new_ano_id = mesh->data(*ano_vh).get_new_id();
			if (_is_boundry_vertex(*ano_vh, mesh, vector_loop))	//1-ring有边界点，根据之前的值，确定边界的情况
			{
				int temp_value = mesh->data(*iv).get_boundary();
				if (temp_value == 1)
					b(new_id) += (-1)*temp_a;
			}

			else
				charac_matrix.insert(new_id, new_ano_id) = temp_a;
		}
		charac_matrix.insert(new_id, new_id) = closed_matrix.coeff(iv->idx(), iv->idx());

	}
	VectorXd xyz = _solve_0_form(charac_matrix, b);
	return xyz;
}

template<class TMesh>
VectorXd HarmonicForm<TMesh>::_solve_0_form(SparseMatrix<double, Eigen::RowMajor> charac_matrix, VectorXd b)				// 求解特征1形式
{
	VectorXd xyz;
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>MatricesCholesky(charac_matrix);
	xyz = MatricesCholesky.solve(b);
	return xyz;
}

template<class TMesh>
VectorXd HarmonicForm<TMesh>::_solve_alpha(TMesh *mesh)
{
	VectorXd alpha(v_size);
	//1 
	Eigen::SparseMatrix<double> d0_transpose = d0.transpose();
	Eigen::VectorXd temp_closed_b = d0_transpose * star1*characteristic_w;
	Eigen::SparseMatrix<double> la = d0_transpose * star1*d0;

	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>MatricesCholesky(la);
	alpha = MatricesCholesky.solve(temp_closed_b);
	return alpha;
}

template<class TMesh>
void HarmonicForm<TMesh>::_calc_uv_on_vertex(TMesh *mesh, int id, string uv_sym)
{
	queue<VertexHandle> Q_v;

	auto vi = mesh->vertices_begin();


	mesh->data(*vi).set_labeled(!curr_label);
	Q_v.push(*vi);
	if (uv_sym == "u")
		mesh->data(*vi).set_u(0);
	else if (uv_sym == "v")
		mesh->data(*vi).set_v(0);


	int v_id = vi->idx();
	//uv.row(v_id) << 0, 0;

	int a = 0;
	vector<EdgeHandle> vec_eh;


	while (!Q_v.empty())
	{
		auto v_h = Q_v.front();
		double value_id;
		if (uv_sym == "u")				value_id = mesh->data(v_h).get_u();
		else if (uv_sym == "v")			value_id = mesh->data(v_h).get_v();

		Q_v.pop();

		bool ver_flag = _is_boundry_vertex(v_h, mesh, vector_loop);
		if (ver_flag)			continue;
		int a = 0;
		for (auto ve = mesh->ve_begin(v_h); ve != mesh->ve_end(v_h); ve++)
		{
			int edge_id = ve->idx();
			auto half_e = mesh->halfedge_handle(*ve, 0);
			auto v1 = mesh->to_vertex_handle(half_e);//1
			auto v2 = mesh->from_vertex_handle(half_e);//-1
			auto another_v = v1;
			double add_or_subtract = 1;					//target 1
			if (another_v == v_h)
			{
				another_v = v2;
				add_or_subtract = -1;					// source -1

			}
			int another_id = another_v.idx();
			double temp_value = 0;
			if (uv_sym == "u")				temp_value = harmonic_w(id, edge_id);// mesh->data(*ve).get_w();
			else if (uv_sym == "v")			temp_value = mesh->data(*ve).get_v();

			double temp_edge_value = 0;

			temp_edge_value = value_id + add_or_subtract * (temp_value);
			if (mesh->data(another_v).is_labeled() == !curr_label)
				continue;
			if (uv_sym == "u")
			{
				mesh->data(another_v).set_u(temp_edge_value);
				//uv.row(another_id) << temp_edge_value, 0;
			}
			else if (uv_sym == "v")
			{
				mesh->data(another_v).set_v(temp_edge_value);
				//double u = mesh->data(another_v).get_u();
				//uv.row(another_id) << u, temp_edge_value;				
			}
			mesh->data(another_v).set_labeled(!curr_label);
			Q_v.push(another_v);
		}
	}
	curr_label = !curr_label;   //reverise it when calcing another edge
}

template<class TMesh>
void HarmonicForm<TMesh>::_calc_closed_1_form(TMesh *mesh, VectorXd xyz)
{
	VectorXd f_v(v_size);

	for (auto iv = mesh->vertices_begin(); iv != mesh->vertices_end(); ++iv)
	{
		int v_id = iv->idx();
		bool flag = _is_boundry_vertex(*iv, mesh, vector_loop);
		double value = xyz(mesh->data(*iv).get_new_id());
		f_v(v_id) = flag == 0 ? value : 0;

	}

	for (auto ei = mesh->edges_begin(); ei != mesh->edges_end(); ei++)
	{
		int edge_id = ei->idx();
		auto half_e = mesh->halfedge_handle(*ei, 0);
		auto v_to = mesh->to_vertex_handle(half_e);//1
		auto v_from = mesh->from_vertex_handle(half_e);//-1
		int id_to = v_to.idx();
		int id_from = v_from.idx();
		int boundary_to = mesh->data(v_to).get_boundary(id_from);
		int boundary_from = mesh->data(v_from).get_boundary(id_to);
		double value_to = f_v(id_to);
		double value_from = f_v(id_from);

		// 0 左边界 1 右边界 
		if (_is_boundry_edge(*ei))
			characteristic_w(edge_id) = 0;
		else if (boundary_from == 0 && _is_boundry_vertex(v_to, mesh, vector_loop))
			characteristic_w(edge_id) = -value_from;
		else if (boundary_from == 1 && _is_boundry_vertex(v_to, mesh, vector_loop))
			characteristic_w(edge_id) = 1 - value_from;
		else if (boundary_to == 0 && _is_boundry_vertex(v_from, mesh, vector_loop))
			characteristic_w(edge_id) = value_to;
		else if (boundary_to == 1 && _is_boundry_vertex(v_from, mesh, vector_loop))
			characteristic_w(edge_id) = value_to - 1;
		else
			characteristic_w(edge_id) = value_to - value_from;
	}

}

template<class TMesh>
void HarmonicForm<TMesh>::_calc_harmonic_1_form(TMesh *mesh, int id)
{
	VectorXd alpha = _solve_alpha(mesh);

	Eigen::VectorXd temp_w = d0 * alpha;
	// w = da + (d*)b + r
	harmonic_w.row(id) = characteristic_w - temp_w;

	_calc_uv_on_vertex(mesh, id, "u");

}

template<class TMesh>
void HarmonicForm<TMesh>::_calc_u_on_vertex(TMesh *mesh, int id)
{
	_find_boundary(mesh, id);
	VectorXd xyz = _constrained_laplacian_matrix(mesh, id);
	_calc_closed_1_form(mesh, xyz);
	_calc_harmonic_1_form(mesh, id);
}

template<class TMesh>
VectorXd HarmonicForm<TMesh>::_calc_area(TMesh *mesh)
{
	Eigen::VectorXd area(f_size);
	for (auto fi = mesh->faces_begin(); fi != mesh->faces_end(); ++fi)
	{
		double a[3];
		int i = 0;
		for (auto eh = mesh->fe_begin(*fi); eh != mesh->fe_end(*fi); ++eh)
		{
			a[i] = mesh->calc_edge_length(*eh);
			i++;
		}
		double data1 = a[0] + a[1] + a[2];
		double data2 = a[0] + a[1] - a[2];
		double data3 = a[0] - a[1] + a[2];
		double data4 = (-1)*a[0] + a[1] + a[2];
		double single_area = 0.25*sqrt(data1*data2*data3*data4);
		area(fi->idx()) = single_area;
	}
	return area;
}

template<class TMesh>
MatrixXd HarmonicForm<TMesh>::_calc_vector_2_form(TMesh *mesh, Eigen::MatrixXd N, Eigen::VectorXd area, int id)
{
	MatrixXd star_w(f_size, 3), wedge(f_size, 3);
	//star_w.resize(f_size, 3);

	for (auto fi = mesh->faces_begin(); fi != mesh->faces_end(); ++fi)
	{
		OpenMesh::Vec3d v_point[3];
		int faceId = fi->idx();
		Eigen::Vector3d faceN = N.row(faceId);


		double single_area = area(faceId);
		wedge.row(faceId) << 0, 0, 0;
		double temp = 0;
		for (auto hi = mesh->fh_begin(*fi); hi != mesh->fh_end(*fi); ++hi)
		{
			auto eh = mesh->edge_handle(*hi);
			int e_id = eh.idx();
			auto halfEdge1 = mesh->halfedge_handle(eh, 0);
			//double temp = all_edge_u[id][eId];
			temp = harmonic_w(id, e_id);
			if (halfEdge1 != *hi)
			{
				temp = (-1) * temp;
			}
			auto nextHandle = mesh->next_halfedge_handle(*hi);
			VertexHandle nextVer = mesh->to_vertex_handle(nextHandle);
			auto point = mesh->point(nextVer);
			Eigen::Vector3d  v_point(point[0], point[1], point[2]);

			wedge.row(faceId) += (-1)*temp / (single_area)*v_point.cross(faceN);

		}
		Eigen::Vector3d temp_wedge = wedge.row(faceId);
		star_w.row(faceId) = faceN.cross(temp_wedge);

	}

	return star_w;
}

template<class TMesh>
void HarmonicForm<TMesh>::_calc_dual_harmonic_1_form(TMesh *mesh, Eigen::MatrixXd star_w)
{
	vector<double > v;
	for (auto ei = mesh->edges_begin(); ei != mesh->edges_end(); ++ei)
	{
		auto halfEdge0 = mesh->halfedge_handle(*ei, 0);
		auto halfEdge1 = mesh->halfedge_handle(*ei, 1);

		auto faceHandle0 = mesh->face_handle(halfEdge0);
		auto faceHandle1 = mesh->face_handle(halfEdge1);


		int idFace0 = faceHandle0.idx();
		int idFace1 = faceHandle1.idx();

		auto v_to = mesh->to_vertex_handle(halfEdge0);
		auto v_from = mesh->from_vertex_handle(halfEdge0);
		auto point = mesh->point(v_to) - mesh->point(v_from);

		Eigen::Vector3d edge_point(point[0], point[1], point[2]);

		double temp = star_w.row(idFace0).dot(edge_point) + star_w.row(idFace1).dot(edge_point);
		temp = temp / 4.0;
		v.push_back(temp);
		mesh->data(*ei).set_v(temp);

	}
	all_edge_v.push_back(v);			//所有的v

}

template<class TMesh>
void HarmonicForm<TMesh>::_calc_v_on_vertex(TMesh *mesh, int id)
{

	area = _calc_area(mesh);

	MatrixXd star_w = _calc_vector_2_form(mesh, N_faces, area, id);

	_calc_dual_harmonic_1_form(mesh, star_w);

	_calc_uv_on_vertex(mesh, id, "v");

}

template<class TMesh>
void HarmonicForm<TMesh>::_save_uv(TMesh *mesh, string uv_sym)
{
	vector<double> ver_u;
	vector<double> ver_v;
	if (uv_sym == "u")
	{
		for (auto vi = mesh->vertices_begin(); vi != mesh->vertices_end(); ++vi)
		{
			ver_u.push_back(mesh->data(*vi).get_u());
		}
		all_ver_u.push_back(ver_u);
	}
	else if (uv_sym == "v")
	{
		for (auto vi = mesh->vertices_begin(); vi != mesh->vertices_end(); ++vi)
		{
			ver_v.push_back(mesh->data(*vi).get_v());
		}
		all_ver_v.push_back(ver_v);
	}
}


template<class TMesh>
void HarmonicForm<TMesh>::build_uv(TMesh *mesh, igl::viewer::Viewer &viewer, int id)
{
	_init_all_data(mesh, viewer, id);
	for (int id = 0; id < loop_size; id++)
	{
		_calc_u_on_vertex(mesh, id);
		_save_uv(mesh, "u");
	}

	for (int id = 0; id < loop_size; id++)
	{
		_calc_v_on_vertex(mesh, id);

		_save_uv(mesh, "v");
	}

}

template<class TMesh>
void HarmonicForm<TMesh>::_set_id_uv(int id, string uv_sym, MatrixXd & uv)
{
	if (uv_sym == "uv")
		for (int i = 0; i < v_size; i++)
			uv.row(i) << all_ver_u[id][i], all_ver_v[id][i];

	else if (uv_sym == "u")
		for (int i = 0; i < v_size; i++)
			uv.row(i) << all_ver_u[id][i], 0.1;

	else if (uv_sym == "v")
		for (int i = 0; i < v_size; i++)
			uv.row(i) << 0.1, all_ver_v[id][i];
}

template<class TMesh>
void HarmonicForm<TMesh>::init(TMesh *mesh, int id, igl::viewer::Viewer &viewer, string uv_sym)
{

	Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> R, G, B, A;
	bool flag = igl::png::readPNG(TUTORIAL_SHARED_PATH "/checker_512.png", R, G, B, A);
	MatrixXd uv(v_size, 2);
	_set_id_uv(id, uv_sym, uv);
	viewer.data.set_uv(uv);
	//Draw checkerboard texture
	viewer.core.show_texture = true;
	// Disable wireframe
	viewer.core.show_lines = false;
	viewer.data.set_texture(R, G, B);

}



#endif // !_HARMONIC_FORM_
