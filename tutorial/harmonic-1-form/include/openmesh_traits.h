#pragma once
#ifndef _OPENMESH_TRAITS_H_
#define _OPENMESH_TRAITS_H_

#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/Traits.hh>
#include <vector>
//#include <iostream>
using namespace OpenMesh;
struct MyTraits :public OpenMesh::DefaultTraits
{
	FaceTraits
	{
	private:
		bool _is_labeled;
	public:
		FaceT() :_is_labeled(false) {};
		const bool &is_labeled() const { return _is_labeled; }
		void set_labeled(const bool & _l) { _is_labeled = _l; }
	};
	EdgeTraits
	{
	private:
		bool _is_labeled;
		double  w;
		double  v;
	public:
		EdgeT() :_is_labeled(false),w(0),v(0) {};
		const bool &is_labeled() const { return _is_labeled; }
		void set_labeled(const bool & _l) { _is_labeled = _l; }

		void set_w(const double &_w) { w = _w; }
		const double get_w() { return w; }

		void set_v(const double &_v) { v = _v; }
		const double get_v() { return v; }

	};
	VertexTraits
	{
	private:
		double u;
		double v;
		bool _is_labeled;
		int new_id;
		int boundary;
		std::vector <int> vec_boundary;
		std::vector <int> vec_b_id;//boundary对应的边界点
	public:
		VertexT() : u(-1), v(0), _is_labeled(false),new_id(0),boundary(-1) {};

		const bool &is_labeled() { return _is_labeled; }
		void set_labeled(const bool & _l) { _is_labeled = _l; }

		void set_u(const double &_u) { u = _u; }
		const double get_u() { return u; }

		void set_v(const double &_v) { v = _v; }
		const double get_v() { return v; }

		void set_new_id(const int &id) { new_id = id; }
		const int get_new_id() { return new_id; }

		void set_boundary(const int &value) { boundary = value; }
		const int get_boundary() { return boundary; }

		void set_boundary(const int &value, const int &id) 
		{
			vec_boundary.push_back(value); vec_b_id.push_back(id); boundary = value;
		}
		const int get_boundary(const int & id) 
		{
			for (int i = 0; i < vec_b_id.size();i++)
			{
				if (vec_b_id[i] == id)
				{
					return vec_boundary[i];
				}
			}
			return -1; 
		}
		std::vector<int> get_boundary_id()
		{
			return vec_boundary;
		}
		const void set_boundary_clear() { vec_boundary.clear(); vec_b_id.clear(); }
	};
};
#endif // !_OPENMESH_TRAITS_H_
