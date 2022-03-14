#include <iostream>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <igl/viewer/Viewer.h>
#include <nanogui/formhelper.h>
#include <nanogui/screen.h>
#include "MeshViewer.h"
#include "openmesh_traits.h"
#include "HarmonicForm.h"

using namespace std;
using namespace igl;
using namespace OpenMesh;
typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits> MyMesh;

int main()
{
	MyMesh mesh;

	OpenMesh::IO::Options opt;
	
	//const char* filename = "torus_cut_graph.obj";
	//const char* filename = TUTORIAL_SHARED_PATH "/2-torus.obj";
	const char* filename = TUTORIAL_SHARED_PATH "/torus_cut_graph.obj";
	if (!OpenMesh::IO::read_mesh(mesh, filename, opt))
	{
		cerr << "Error: Cannot read mesh from " << filename << endl;
		
	}

	igl::viewer::Viewer viewer;
	convert_libigl<MyMesh>(mesh,viewer);
	
	init<MyMesh>(mesh);
	

	viewer.callback_init = &init_new_window;
	viewer.callback_pre_draw = &pre_draw;
	viewer.core.is_animating = true;
	viewer.core.animation_max_fps = 100.;

	//Disable wireframe
	viewer.core.show_lines = false;

	

	viewer.launch();
	

	return 0;
}