#pragma once
#ifndef _CLOSED_LOOP_H_
#define _CLOSED_LOOP_H_


#include <vector>
#include <iostream>
#include "openmesh_traits.h"
using namespace std;
vector<EdgeHandle> vector_tree_eh;				//最小生成树的路径
vector<EdgeHandle> vector_all_eh;				//所有的路径 或者 所有的路径 - 最小生成树的路径
vector<EdgeHandle> vector_false_eh;				//上一次标记为false的路径
vector< vector<EdgeHandle> > vector_all_loop;		//所有loop的路径

template <typename TMesh>  bool find_loop(TMesh *mesh,int &id);
template <typename TMesh>  void spanning_tree(TMesh *mesh);
template <typename TMesh>  void set_mesh_true(TMesh *mesh);
template <typename TMesh>  void set_loop_flag(TMesh *mesh, int &id);		//标记当前路径
template <typename TMesh>  void save_false_flag(TMesh *mesh);
template <typename TMesh>  std::vector< std::vector<EdgeHandle> > get_loop();


template <typename TMesh>
void set_mesh_true(TMesh *mesh)				//上一次被标记为true， 标记为false；
{
	for (auto eh = mesh->edges_begin(); eh != mesh->edges_end(); ++eh)
	{
		mesh->data(*eh).set_labeled(true);
	}
}

template <typename TMesh>  void set_loop_true(TMesh *mesh)						//整个mesh标记为true
{
	for (auto vc_it = vector_false_eh.cbegin(); vc_it != vector_false_eh.cend(); vc_it++)
	{
		mesh->data(*vc_it).set_labeled(true);
	}
}
template <typename TMesh> bool is_closed(TMesh *mesh, EdgeHandle eh, int edge_id) //判断是否形成了一个环
{
	bool closed_flag = false;
	int times = 0, flag1 = 0, flag2 = 0;
	auto v1 = mesh->from_vertex_handle(mesh->halfedge_handle(eh, edge_id));
	auto v2 = mesh->to_vertex_handle(mesh->halfedge_handle(eh, edge_id));

	for (auto vc_it = vector_tree_eh.cbegin(); vc_it != vector_tree_eh.cend(); vc_it++)
	{
		auto temp_v1 = mesh->from_vertex_handle(mesh->halfedge_handle(*vc_it, edge_id));
		auto temp_v2 = mesh->to_vertex_handle(mesh->halfedge_handle(*vc_it, edge_id));
		if ((temp_v1 == v1 || temp_v2 == v1) && flag1 == 0)
		{
			times++;
			flag1 = 1;
		}
		else if ((temp_v1 == v2 || temp_v2 == v2) && flag2 == 0)
		{
			times++;
			flag2 = 1;
		}
	}

	//assert(times < 3);

	if (times == 2)
	{
		closed_flag = true;
	}
	return closed_flag;
}

template <typename TMesh> void find_path(TMesh *mesh, EdgeHandle eh, int last_id, int edge_id)//递归找最小生成树的路径
{

	mesh->data(eh).set_labeled(true);				//将mesh上的该边去掉

	if (is_closed(mesh, eh, edge_id))
		return;

	std::vector<EdgeHandle>::iterator iter;

	iter = find(vector_all_eh.begin(), vector_all_eh.end(), eh);
	if (iter != vector_all_eh.end()) {				//在vector中删除该点
		vector_all_eh.erase(iter);
	}


	vector_tree_eh.push_back(eh);					//压到vector中
	auto v1 = mesh->to_vertex_handle(mesh->halfedge_handle(eh, edge_id));
	auto v2 = mesh->from_vertex_handle(mesh->halfedge_handle(eh, edge_id));
	//cout << "v1 = " << v1.idx() << "v2 = " << v2.idx() << endl;
	if (v1.idx() == last_id)
		v1 = v2;

	for (auto ve = mesh->ve_begin(v1); ve != mesh->ve_end(v1); ++ve)
	{
		EdgeHandle now_eh = *ve;
		if (!mesh->data(now_eh).is_labeled())
		{
			find_path(mesh, now_eh, v1.idx(), edge_id);

		}
	}
}
template <typename TMesh> void spanning_tree(TMesh *mesh) // 生成最小生成树
{
	int i = 0;
	int edge_id = 0;
	EdgeHandle eh;
	for (auto ie = mesh->edges_begin(); ie != mesh->edges_end(); ++ie)
	{

		if (!mesh->data(*ie).is_labeled())
		{
			vector_all_eh.push_back(*ie);
		}
	}

	eh = vector_all_eh.front();
	auto last_id = mesh->to_vertex_handle(mesh->halfedge_handle(eh, 0)).idx();
	find_path(mesh, eh, last_id, edge_id);
	//min(1, 2);
}



template <typename TMesh> bool find_loop(TMesh *mesh, int &id)		//找到所有的环路。从中可以提取出loop。
{
	//assert(id < vector_all_eh.size());

	set_loop_true(mesh);

	for (auto ve_it = vector_tree_eh.begin(); ve_it != vector_tree_eh.end(); ve_it++)
	{
		mesh->data(*ve_it).set_labeled(false);
	}
	vector_false_eh.assign(vector_tree_eh.begin(), vector_tree_eh.end());	//复制整个数组

	EdgeHandle eh = vector_all_eh[id];
	mesh->data(eh).set_labeled(false);		//显示出来
	vector_false_eh.push_back(eh);


	id++;

	if (id == vector_all_eh.size())	return true;
	return false;
}

template <typename TMesh> void set_loop_flag(TMesh *mesh, int &id)		//标记当前路径
{
	set_loop_true(mesh);
	int size = (int)vector_all_eh.size();
	//assert(size > 0);
	id = id % vector_all_eh.size();


	for (auto v_it = vector_all_loop[id].begin(); v_it != vector_all_loop[id].end(); v_it++)
	{
		mesh->data(*v_it).set_labeled(false);
	}

	vector_false_eh.assign(vector_all_loop[id].begin(), vector_all_loop[id].end());	//复制整个数组

}

template <typename TMesh> void save_false_flag(TMesh *mesh)
{
	vector<EdgeHandle> single_loop;
	for (auto v_it = vector_false_eh.begin(); v_it != vector_false_eh.end(); v_it++)
	{
		if (!mesh->data(*v_it).is_labeled())
		{
			single_loop.push_back(*v_it);
		}
	}
	int len = (int)single_loop.size();

	// vector中存放的边是连续的
	vector<EdgeHandle> loop;
	loop.push_back(single_loop[len - 1]);
	single_loop.pop_back();

	for (int i = 0; i < len - 1; i++)
	{
		EdgeHandle ve = loop[i];
		auto half_e = mesh->halfedge_handle(ve, 0);
		auto v1 = mesh->to_vertex_handle(half_e);
		auto v2 = mesh->from_vertex_handle(half_e);

		for (auto iv = single_loop.begin(); iv != single_loop.end(); iv++)
		{

			EdgeHandle ie = *iv;
			auto temp_half_e = mesh->halfedge_handle(ie, 0);
			auto temp_v1 = mesh->to_vertex_handle(temp_half_e);
			auto temp_v2 = mesh->from_vertex_handle(temp_half_e);
			if (v1 == temp_v2 || v1 == temp_v1
				|| v2 == temp_v1 || v2 == temp_v2)
			{
				loop.push_back(ie);

				single_loop.erase(iv);
				break;
			}
		}

	}
	EdgeHandle ve = loop[len - 1];
	auto half_e = mesh->halfedge_handle(ve, 0);
	auto v1 = mesh->to_vertex_handle(half_e);
	auto v2 = mesh->from_vertex_handle(half_e);
	//std::cout << v1.idx() << " " << v2.idx() ;
	vector_all_loop.push_back(loop);
}

template <typename TMesh>  vector< vector<EdgeHandle> > get_loop()
{
	return vector_all_loop;
}

#endif