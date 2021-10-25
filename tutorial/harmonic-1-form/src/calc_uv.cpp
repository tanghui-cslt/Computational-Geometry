#include <igl/viewer/Viewer.h>
#include <igl/cotmatrix.h>
#include <igl/PNG/readPNG.h>
#include <Eigen/Sparse>
#include <Eigen/SparseCore>
#include <Eigen/Cholesky>
#include <queue>
#include <vector>

#include <iostream>
#include <igl/per_face_normals.h>
#include "calc_uv.h"

using namespace std;

static Eigen::SparseMatrix<double, Eigen::RowMajor> charac_matrix;		//计算特征1形式的矩阵
static vector<vector<EdgeHandle> > vector_all_loop;						//所有的loop的vector
static vector<EdgeHandle> vector_loop;									//当前选中的loop的vector
static int v_size;													//点的个数
static int e_size;													//边的个数
static int f_size;													//面的个数
static Eigen::VectorXd xyz;												//计算特征1形式每个点的值
static double min_u = 0, max_u = 0;										//点的最大最小值
static int loop_size = 0;												//当前loop的长度
static Eigen::VectorXd b;												//特征1形式的右端项
static bool curr_label = false;											//当前所有的边所处的状态
static Eigen::SparseMatrix<double> d0;									//d0矩阵
static Eigen::SparseMatrix<double, Eigen::RowMajor> star1;				//*1 矩阵
static vector<VertexHandle> vector_labled_boundary;						//当前选中的loop
static Eigen::SparseMatrix<double> L;									//L = cotangent
static Eigen::VectorXd characteristic_w;								//所有边的特征1形式characteristic_w
static Eigen::VectorXd f_v;												//每个点特征1形式的值
static Eigen::VectorXd harmonic_w;										//每个点的调和1形式
static Eigen::VectorXd alpha;											//alpha矩阵的值
static Eigen::SparseMatrix<double> closed_matrix;						//闭矩阵的左端项
static Eigen::VectorXd closed_b;										//闭矩阵的右端项
static Eigen::MatrixXd uv;												// u每个点的u值，v每个点的v值								
static vector< vector<double> >all_edge_u;
static vector< vector<double> >all_edge_v;	//所有的u值，所有的v值
static vector< vector<double> > all_ver_u;
static vector< vector<double> > all_ver_v;
static int loop_num;
static int now_loop_id;
static Eigen::VectorXd area;									//每个面的面积
static Eigen::MatrixXd N_faces;									//法向量
static Eigen::MatrixXd star_w;									//*w 存放2g个结果
static Eigen::MatrixXd wedge;									// w 存放2g个结果

inline bool is_vertex_exist(VertexHandle vh, MyMesh *mesh, vector<EdgeHandle> loop = vector_loop)
{
	bool flag = false;
	size_t size = loop.size();
	for (size_t i = 0; i < size; i++)
	{
		auto half_0 = mesh->halfedge_handle(loop[i], 0);				//0号半边
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

inline bool is_edge_exist(EdgeHandle eh)
{
	for (auto ei = vector_loop.begin(); ei != vector_loop.end(); ei++)
		if (*ei == eh)		return true;
	return false;
}
void set_false_vector(MyMesh *mesh)
{
	for (auto ve = vector_labled_boundary.begin(); ve != vector_labled_boundary.end(); ve++)
	{
		mesh->data(*ve).set_u(-1);
		mesh->data(*ve).set_boundary(-1);
		mesh->data(*ve).set_boundary_clear();
	}
}

void uv_init_data(MyMesh *mesh, int id)
{

	v_size = (int)mesh->n_vertices();
	e_size = (int)mesh->n_edges();
	f_size = (int)mesh->n_faces();
	min_u = 0;
	max_u = 0;

	//初始化矩阵
	vector_all_loop = get_loop();

	loop_size = (int)vector_all_loop.size();

	id = id % loop_size;

	vector_loop = vector_all_loop[id];

	int len = (int)vector_loop.size();

	//cout << "loop_size = " << loop_size << " v " << v_size << " e " << e_size << " f " << f_size << endl;
	Eigen::SparseMatrix<double>temp_hor(v_size - len, v_size - len);

	charac_matrix = temp_hor;
	d0.resize(e_size, v_size);
	star1.resize(e_size, e_size);
	b = Eigen::VectorXd::Zero(v_size - len);
	f_v.resize(v_size);
	characteristic_w.resize(e_size);
	harmonic_w.resize(e_size);
	alpha.resize(v_size);
	//closed_matrix.resize(v_size, v_size);
	closed_b.resize(v_size);
	uv.resize(v_size, 2);
	star_w.resize(f_size, 3);
	wedge.resize(f_size, 3);
	N_faces.resize(f_size, 3);
}

static void uv_init_matrix_d_star(MyMesh *mesh)
{
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
	closed_matrix = d0_transpose*star1*d0;
}

// 找到loop的左右的两个边界
void uv_splice_loop(MyMesh *mesh, igl::viewer::Viewer &viewer)
{
	set_false_vector(mesh);
	vector<VertexHandle> vector_boundary;
	int len = (int)vector_loop.size();

	//	id = id % len;

	EdgeHandle vh = vector_loop[0];							//取出一条边
	auto half_0 = mesh->halfedge_handle(vh, 0);				//0号半边
	auto temp_half_0 = half_0;								//可移动半边

	auto half_1 = mesh->halfedge_handle(vh, 1);				//1号半边
	auto temp_half_1 = half_1;								//可移动的1号半边
															//cout << "v id = ";
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
		//cout << "id = " << ori_id << " ";

		temp_half_0 = mesh->next_halfedge_handle(temp_half_0);		//这个半边的下一个半边

		while (mesh->edge_handle(temp_half_0) != pre_eh &&			//判断是否找到边界
			mesh->edge_handle(temp_half_0) != next_eh)
		{
			auto v0 = mesh->to_vertex_handle(temp_half_0);				//保存边界
			/*if (v0.idx() == 1700)
			{
				cout << " num 0 halfedge of 1700 ori=" << v0.idx() << " ori = " << ori_id << endl;
			}*/

			if (//mesh->data(v0).get_boundary() == 0 ||
				is_vertex_exist(v0, mesh, vector_loop))
			{
				//cout << "!!!!!!!!!!!!\n";
			}
			else
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

		while (mesh->edge_handle(temp_half_1) != pre_eh &&								//判断是否找到边界
			mesh->edge_handle(temp_half_1) != next_eh)
		{
			v1 = mesh->to_vertex_handle(temp_half_1);				//保存边界


			if (//mesh->data(v1).get_boundary() == 1 ||
				is_vertex_exist(v1, mesh, vector_loop))
			{
			}
			else
			{
				//cout << "------boundary--\n";
				mesh->data(v1).set_boundary(1, ori_id);
				//cout << " boundary : " << mesh->data(v1).get_boundary() << " ";
				mesh->data(v1).set_boundary(1);
				vector_boundary.push_back(v1);
			}

			temp_half_1 = mesh->opposite_halfedge_handle(temp_half_1);
			temp_half_1 = mesh->next_halfedge_handle(temp_half_1);
		}
	}

	vector_labled_boundary.assign(vector_boundary.begin(), vector_boundary.end());
}


void uv_build_matrix_first_field(MyMesh *mesh, igl::viewer::Viewer &viewer)
{
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;
	V = viewer.data.V;
	F = viewer.data.F;

	igl::cotmatrix(V, F, L);
	int count = 0;	//new id
	uv_init_matrix_d_star(mesh);
	for (auto iv = mesh->vertices_begin(); iv != mesh->vertices_end(); ++iv)
	{
		//new id
		if (is_vertex_exist(*iv, mesh, vector_loop)) {
			continue;
		}
		mesh->data(*iv).set_new_id(count);
		count++;
	}
	// 计算左端矩阵特征1形式  & 右端矩阵b
	for (auto iv = mesh->vertices_begin(); iv != mesh->vertices_end(); ++iv)
	{
		if (is_vertex_exist(*iv, mesh, vector_loop))	continue;
		int new_id = mesh->data(*iv).get_new_id();
		int id = iv->idx();
		/* -------
				test_a data
		*/
		int test_a = 0;
		int test_1 = 0;
		int test_0 = 0;
		for (auto ve = mesh->vv_begin(*iv); ve != mesh->vv_end(*iv); ++ve)//l-ring 之和为0
		{
			auto ano_vh = ve;

			int ano_vid = ano_vh->idx();
			double temp_a = closed_matrix.coeff(iv->idx(), ano_vid);
			int new_ano_id = mesh->data(*ano_vh).get_new_id();
			if (is_vertex_exist(*ano_vh, mesh, vector_loop))	//1-ring有边界点，根据之前的值，确定边界的情况
			{
				test_a++;

				int temp_value = mesh->data(*iv).get_boundary(ano_vid);
				if (temp_value == 1)
				{
					b(new_id) += (-1)*temp_a;
					test_1++;
					//	cout << " b = " << b(new_id)<<" id="<<ano_vid;
				}
				else {// 测试有一部分为0， 一部分为1的情况
					test_0++;
				}
				/*if (test_1 * test_0 != 0)
				{
					cout << "id = " << id << "  1-ring the number of zero =" << test_a << " 1=" << test_1 << " 0 =" << test_0 << endl;
				}*/
			}
			else
			{
				charac_matrix.insert(new_id, new_ano_id) = temp_a;
			}
		}
		charac_matrix.insert(new_id, new_id) = closed_matrix.coeff(iv->idx(), iv->idx());

	}
}
static void solve_charasteristic()				// 求解特征1形式
{
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>MatricesCholesky(charac_matrix);
	xyz = MatricesCholesky.solve(b);
}
static void solve_alpha(MyMesh *mesh)
{
	//1 
	Eigen::SparseMatrix<double> d0_transpose = d0.transpose();
	Eigen::VectorXd temp_closed_b = d0_transpose*star1*characteristic_w;
	Eigen::SparseMatrix<double> la = d0_transpose*star1*d0;

	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>MatricesCholesky(la);
	alpha = MatricesCholesky.solve(temp_closed_b);
}

//id = 0 计算每个点的u值
//id = 1 计算每个点的v值
static void traveal_mesh_uv(MyMesh *mesh, int id = 0)
{
	queue<VertexHandle> Q_v;

	auto vi = mesh->vertices_begin();

	if (id != 0 && id != 1)
	{
		cout << "id = " << id << " ,输入id有错,无法计算每个点的uv值\n";
		return;
	}
	mesh->data(*vi).set_labeled(!curr_label);
	Q_v.push(*vi);
	if (id == 0)
	{
		mesh->data(*vi).set_u(0);
		//ver_u.push_back(0);
	}
	else
	{
		mesh->data(*vi).set_v(0);
		//ver_v.push_back(0);
	}


	int v_id = vi->idx();
	uv.row(v_id) << 0, 0;
	// -----
	//cout << "zero id = ";
	int a = 0;
	vector<EdgeHandle> vec_eh;

	// -----
	while (!Q_v.empty())
	{
		auto v_h = Q_v.front();
		double value_id;
		if (id == 0)		value_id = mesh->data(v_h).get_u();
		else			value_id = mesh->data(v_h).get_v();

		Q_v.pop();
		//if (fabs(u) < 1e-7) { cout << v_h.idx() << " "; }
		bool ver_flag = is_vertex_exist(v_h, mesh, vector_loop);
		if (ver_flag )
		{
			continue;
		}
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
			if (id == 0)		temp_value = mesh->data(*ve).get_w();
			else				temp_value = mesh->data(*ve).get_v();

			double temp_edge_value = 0;

			temp_edge_value = value_id + add_or_subtract*(temp_value);
			if (mesh->data(another_v).is_labeled() == !curr_label)
			{
				continue;
			}
			if (id == 0)
			{
				//double u = mesh->data(another_v).get_u();
				//int boundary = mesh->data(v_h).get_boundary();
				//bool boundary_flag = is_vertex_exist(another_v, mesh, vector_loop);
				//
				//if (boundary == 1 && boundary_flag)
				//{
				//	//cout << " u=" << u << " ";
				//	
				//}
				//
				//else
				//{
				//	mesh->data(another_v).set_u(temp_edge_value);
				//}
				mesh->data(another_v).set_u(temp_edge_value);
				uv.row(another_id) << temp_edge_value, 0;
				//ver_u.push_back(temp_edge_value);
			}
			else
			{
				
				//if(fabs(u) <0.00000001 )

				mesh->data(another_v).set_v(temp_edge_value);
				double u = mesh->data(another_v).get_u();
				uv.row(another_id) << u, temp_edge_value;
				//ver_v.push_back(temp_edge_value);
			}
			mesh->data(another_v).set_labeled(!curr_label);

			Q_v.push(another_v);
		}
	}


	curr_label = !curr_label;   //当前遍历与之前的不同
}

static void calc_w(MyMesh *mesh)
{
	solve_charasteristic();
	for (auto iv = mesh->vertices_begin(); iv != mesh->vertices_end(); ++iv)
	{
		int v_id = iv->idx();
		if (!(is_vertex_exist(*iv, mesh, vector_loop)))
		{
			double value = xyz(mesh->data(*iv).get_new_id());
			f_v(v_id) = value;
		}
		else
		{
			f_v(v_id) = 0;
		}

	}

	//int a = 0;
	//int b = 0, c = 0, d = 0, 
	int len = 0;
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
		if (is_edge_exist(*ei))
			characteristic_w(edge_id) = 0, len++;
		else if (boundary_from == 0 && is_vertex_exist(v_to, mesh))
			characteristic_w(edge_id) = -value_from;
		else if (boundary_from == 1 && is_vertex_exist(v_to, mesh))
			characteristic_w(edge_id) = 1 - value_from;
		else if (boundary_to == 0 && is_vertex_exist(v_from, mesh))
			characteristic_w(edge_id) = value_to;
		else if (boundary_to == 1 && is_vertex_exist(v_from, mesh))
			characteristic_w(edge_id) = value_to - 1;
		else
			characteristic_w(edge_id) = value_to - value_from;
	}

	for (auto fi = mesh->faces_begin(); fi != mesh->faces_end(); fi++)
	{
		double sum = 0;
		for (auto fh = mesh->fh_begin(*fi); fh != mesh->fh_end(*fi); fh++)
		{
			auto eh = mesh->edge_handle(*fh);
			auto half_e = mesh->halfedge_handle(eh, 0);
			int  id = eh.idx();
			double temp = characteristic_w(id);
			if (half_e != *fh)
			{
				temp = (-1)*temp;
			}
			sum += temp;
		}
		if (fabs(sum) > 1e-4)
		{
			getchar();
		}

	}
}
void calc_u(MyMesh *mesh)
{
	solve_alpha(mesh);
	Eigen::VectorXd temp_w = d0*alpha;
	vector<double > u;
	for (auto eh = mesh->edges_begin(); eh != mesh->edges_end(); eh++)
	{
		int e_id = eh->idx();
		double temp = characteristic_w(e_id) - temp_w(e_id);
		u.push_back(temp);
		harmonic_w(e_id) = temp;
		mesh->data(*eh).set_w(temp);
	}
	all_edge_u.push_back(u);				//所有的u的值
	//static int a = 0;
	//for (auto fi = mesh->faces_begin(); fi != mesh->faces_end(); fi++)
	//{
	//	a++;
	//	double sum = 0;
	//	for (auto fh = mesh->fh_begin(*fi); fh != mesh->fh_end(*fi); fh++)
	//	{
	//		auto eh = mesh->edge_handle(*fh);
	//		auto half_e = mesh->halfedge_handle(eh, 0);
	//		int  id = eh.idx();
	//		double temp = mesh->data(eh).get_w();
	//		if (half_e != *fh)
	//		{
	//			temp = (-1)*temp;
	//		}
	//		sum += temp;
	//	}
	//	//cout << sum << " ";
	//	if (fabs(sum) > 1e-4)
	//	{
	//		getchar();
	//	}

	//}
	traveal_mesh_uv(mesh);


	/*min_u = 0;
	max_u = 0;
	for (auto iv = mesh->vertices_begin(); iv != mesh->vertices_end(); ++iv)
	{
		int v_id = iv->idx();
		double value = mesh->data(*iv).get_u();
		if (value > max_u)	max_u = value;
		if (value < min_u)	min_u = value;

	}*/
	//cout << "min " << min_u << " max = " << max_u << "\n";

}

Eigen::VectorXd calc_area(MyMesh *mesh)
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

Eigen::MatrixXd calc_wedge(MyMesh *mesh, Eigen::MatrixXd N, Eigen::VectorXd area)
{
	//Eigen::MatrixXd star_w(f_size, 3);
	//Eigen::MatrixXd wedge(f_size, 3);
	//auto test_w = characteristic_w;			//w为闭1形式
	auto test_w = harmonic_w;				//w为调和一形式
	for (auto fi = mesh->faces_begin(); fi != mesh->faces_end(); ++fi)
	{
		OpenMesh::Vec3d v_point[3];
		int faceId = fi->idx();
		Eigen::Vector3d faceN = N.row(faceId);

		//cout << "faceN = " << faceN << endl;
		double single_area = area(faceId);
		wedge.row(faceId) << 0, 0, 0;
		for (auto hi = mesh->fh_begin(*fi); hi != mesh->fh_end(*fi); ++hi)
		{
			auto eh = mesh->edge_handle(*hi);
			int eId = eh.idx();
			auto halfEdge1 = mesh->halfedge_handle(eh, 0);
			double temp = test_w(eId);
			if (halfEdge1 != *hi)
			{
				temp = (-1) * temp;
			}
			auto nextHandle = mesh->next_halfedge_handle(*hi);
			VertexHandle nextVer = mesh->to_vertex_handle(nextHandle);
			auto point = mesh->point(nextVer);
			Eigen::Vector3d  v_point(point[0], point[1], point[2]);

			wedge.row(faceId) += (-1)*temp / (2.0*single_area)*v_point.cross(faceN);
		}
		Eigen::Vector3d temp_wedge = wedge.row(faceId);
		star_w.row(faceId) = faceN.cross(temp_wedge);
	}

	return star_w;
}
void calc_edge_v(MyMesh *mesh, Eigen::MatrixXd star_w)
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

		Eigen::Vector3d edge_point(point[0], point[1], point[3]);

		double temp = star_w.row(idFace0).dot(edge_point) + star_w.row(idFace1).dot(edge_point);
		v.push_back(temp / 2.0);
		mesh->data(*ei).set_v(temp / 2.0);

	}
	all_edge_v.push_back(v);			//所有的v
}
void optimized_v(MyMesh *mesh, Eigen::MatrixXd N_faces, Eigen::MatrixXd star_w)
{
	cout << loop_size << " ";

	Eigen::MatrixXd A(loop_size, loop_size);
	for (size_t i = 0; i < loop_size; i++)
	{
		for (size_t j = i; j < loop_size; j++)
		{
			double s = 0;
			for (auto fi = mesh->faces_begin(); fi != mesh->faces_end(); ++fi)
			{

			}
		}
	}
}
void calc_v(igl::viewer::Viewer & viewer, MyMesh *mesh)
{
	//Eigen::MatrixXd N_faces;
	//Eigen::MatrixXd star_w(f_size, 3);
	igl::per_face_normals(viewer.data.V, viewer.data.F, N_faces);

	area = calc_area(mesh);

	star_w = calc_wedge(mesh, N_faces, area);

	calc_edge_v(mesh, star_w);

	//optimized_v(mesh, N_faces, star_w);

	traveal_mesh_uv(mesh, 1);
	double min_v = 0;
	double max_v = 0;
	for (auto iv = mesh->vertices_begin(); iv != mesh->vertices_end(); ++iv)
	{
		int v_id = iv->idx();
		double value = mesh->data(*iv).get_v();
		if (value > max_v)	max_v = value;
		if (value < min_v)	min_v = value;

	}
	//cout << "minv " << min_v << " maxv = " << max_v << "\n";
}
Eigen::VectorXd normalized_vec(Eigen::VectorXd temp_vec)
{
	double min = 0, max = 1;
	for (int i = 0; i < temp_vec.size(); i++)
	{
		double temp = temp_vec(i);
		if (min > temp)	min = temp;
		if (max < temp) max = temp;
	}
	double len = max - min;
	Eigen::VectorXd new_vec(temp_vec.size());
	for (int i = 0; i < temp_vec.size(); i++)
	{
		new_vec(i) = (temp_vec(i) - min) / len;
	}

	return new_vec;
}
void save_uv(MyMesh *mesh, int id)
{
	int i = 0;
	vector<double> ver_u;
	vector<double> ver_v;
	for (auto vi = mesh->vertices_begin(); vi != mesh->vertices_end(); ++vi)
	{
		ver_u.push_back(mesh->data(*vi).get_u());
		ver_v.push_back(mesh->data(*vi).get_v());
		//ver_v.push_back(0);
		i++;
	}
	all_ver_u.push_back(ver_u);
	all_ver_v.push_back(ver_v);
}

void build_uv(MyMesh *mesh, int id, igl::viewer::Viewer &viewer)
{
	uv_init_data(mesh, id);
	uv_splice_loop(mesh, viewer);
	uv_build_matrix_first_field(mesh, viewer);

	loop_num = id;
	now_loop_id = id;
	calc_w(mesh);
	calc_u(mesh);
	calc_v(viewer, mesh);
	save_uv(mesh, id);
}

void set_id_uv(int id,string uv_sym)
{	
	if (uv_sym == "uv")
	{
		for (int i = 0; i < v_size; i++)
		{
			uv.row(i) << all_ver_u[id][i], all_ver_v[id][i];
		}

	}
	else if (uv_sym == "u")
	{
		for (int i = 0; i < v_size; i++)
		{
			uv.row(i) << all_ver_u[id][i], 0.1;
		}
	}
}
void init(MyMesh *mesh, int id, igl::viewer::Viewer &viewer, string uv_sym)
{

	//optimized_v(mesh, N_faces, star_w);
	Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> R, G, B, A;
	bool flag = igl::png::readPNG(TUTORIAL_SHARED_PATH "/checker_512.png", R, G, B, A);
	/*Eigen::VectorXd temp_u = uv.col(0);
	Eigen::VectorXd temp_v = normalized_vec(temp_u);
	uv.col(0) = temp_v;

	temp_u = uv.col(1);
	temp_v = normalized_vec(temp_u);*/
	//uv.col(1) = temp_v;
	/*for (int i = 0; i < temp_v.size(); i++)
	{
		double temp_a = temp_v[i];
		if (temp_a > 1 || temp_a < 0)	cout << " value = " << temp_a;
	}*/
	set_id_uv(id, uv_sym);
	viewer.data.set_uv(uv);
	//Draw checkerboard texture
	viewer.core.show_texture = true;
	// Disable wireframe
	viewer.core.show_lines = false;
	viewer.data.set_texture(R, G, B);
}



