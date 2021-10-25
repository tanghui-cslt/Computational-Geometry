#include <iostream>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <igl/viewer/Viewer.h>
#include <nanogui/formhelper.h>
#include <nanogui/screen.h>
#include "convert_to_libigl.h"
#include "openmesh_traits.h"

using namespace std;
using namespace igl;

int main()
{
	MyMesh mesh;
	OpenMesh::IO::Options opt;
	
	//const char* filename = "torus_cut_graph.obj";
	//const char* filename = TUTORIAL_SHARED_PATH "/3holes.off";
	const char* filename = TUTORIAL_SHARED_PATH "/2-torus.obj";
	if (!OpenMesh::IO::read_mesh(mesh, filename, opt))
	{
		cerr << "Error: Cannot read mesh from " << filename << endl;
		
	}
	//if (!OpenMesh::IO::read_mesh(mesh, filename, opt))
	//{
	//	cerr << "Error: Cannot read mesh from " << filename << endl;

	//}
	igl::viewer::Viewer viewer;
	convert_libigl(mesh,viewer);
	
	init(mesh);
	
	viewer.callback_init = &init_new_window;
	viewer.callback_pre_draw = &pre_draw;
	viewer.core.is_animating = true;
	viewer.core.animation_max_fps = 100.;

	//Disable wireframe
	viewer.core.show_lines = false;

	

	viewer.launch();
	

	return 0;
}