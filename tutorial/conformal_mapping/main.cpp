#include <igl/readOFF.h>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/cotmatrix.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/viewer/Viewer.h>
#include <igl/per_vertex_normals.h>
#include <igl/read_triangle_mesh.h>
#include <iostream>
#include <nanogui/formhelper.h> 
#include <nanogui/screen.h>		

#include <time.h>

using namespace std;
using namespace Eigen;

Eigen::MatrixXd V;
Eigen::MatrixXi F;

Eigen::MatrixXd Laplacian_matrix;						
Eigen::MatrixXd N_vertices;								
SparseMatrix<double> L;									//cot Matrix

bool stop_flag = true;
double energy_init = 0;								
double energy_times_i = 0;							
double energy_times_i1 = 0;							

bool show_gauss = false;
void calc_laplacian();								

double  step_size = 0.1;							
double threshold_value = 1e-5;						
bool boolVariable = true;
static int times = 0;								
									
clock_t start, finish;								
double totaltime;

bool GUI(igl::viewer::Viewer &viewer)
{

	viewer.ngui->setLabelFontSize(20);
	viewer.ngui->addWindow(Eigen::Vector2i(300, 10), "Control Panel");

	// ... or using a custom callback

	viewer.ngui->addButton("gauss mapping", [&]() {
		viewer.data.clear();
		viewer.data.set_mesh(N_vertices, F);
		viewer.core.align_camera_center(N_vertices, F);
	});

	viewer.ngui->addButton("calc Laplacian equation", [&]() {
		viewer.core.is_animating = true;
	});

	viewer.ngui->addButton("PAUSE", [&]() {
		viewer.core.is_animating = false;
	});

	
	viewer.ngui->setLabelFontSize(20);
	viewer.ngui->addVariable("step size :", step_size);
	viewer.ngui->addVariable("threshold value:", threshold_value);
	//viewer.ngui->addVariable("iteration times:", times);
	viewer.ngui->addVariable("total energy:", energy_times_i1);

	finish = clock();
	totaltime = (double)(finish - start) / CLOCKS_PER_SEC;
	viewer.ngui->addVariable("run time (unit of time : s):", totaltime);
	// Generate menu
	viewer.screen->performLayout();

	return false;
}

bool key_down(igl::viewer::Viewer &viewer, unsigned char key, int mods)
{
	switch (key)
	{
	case ' ':
		viewer.core.is_animating = !viewer.core.is_animating;
		return true;
	}
	return false;
}

bool pre_draw(igl::viewer::Viewer & viewer)
{

	if (viewer.core.is_animating)
	{

		if (stop_flag)
		{
			times++;
			calc_laplacian();
			//	GUI(viewer);
			viewer.data.clear();
			viewer.data.set_mesh(N_vertices, F);
			viewer.core.align_camera_center(N_vertices, F);

			double diff_e = energy_times_i - energy_times_i1;
			if (fabs(diff_e)< 1e-5)
			{
				cout << "Stop" << endl;
				stop_flag = false;


			}			
			finish = clock();
			totaltime = (double)(finish - start) / CLOCKS_PER_SEC;
		}
	}
	//GUI(viewer);
	finish = clock();
	totaltime = (double)(finish - start) / CLOCKS_PER_SEC;
	//viewer.ngui->addVariable("run time (unit of time : s):", totaltime);
	viewer.ngui->refresh();
	return false;
}

int main(int argc, char *argv[])
{
	// Load a mesh in OFF format
	//std::string filename = TUTORIAL_SHARED_PATH "/.off";
	//std::string filename =TUTORIAL_SHARED_PATH "/truck.obj";
	//std::string filename = "cube.off";

	start = clock();
	//igl::readOFF(filename, V, F);
	
	igl::readOBJ(TUTORIAL_SHARED_PATH "/cow.obj", V, F);

	// Plot the mesh
	igl::viewer::Viewer viewer;
	
	igl::per_vertex_normals(V, F, N_vertices);
	igl::cotmatrix(V, F, L);

	Laplacian_matrix.resize(V.rows(), 3);


	viewer.callback_pre_draw = &pre_draw;
	viewer.callback_key_down = &key_down;

	viewer.callback_init = &GUI;

	viewer.data.clear();
	viewer.data.set_mesh(V, F);
	viewer.core.align_camera_center(V, F);
	viewer.launch();

}



void calc_laplacian()
{

	int row_len = L.rows();
	int col_len = L.cols();
	energy_times_i = energy_times_i1;
	energy_times_i1 = 0;
	Vector3d  centroid_vector(0, 0, 0);

	
	for (int i = 0; i < row_len; i++)
	{
		double sum_energy = 0;
		Laplacian_matrix(i, 0) = 0;
		Laplacian_matrix(i, 1) = 0;
		Laplacian_matrix(i, 2) = 0;

		for (int j = 0; j < col_len; j++)
		{

			double cot = L.coeff(i, j);
			if (fabs(cot) < 1e-4)
			{
			}
			else if (i != j)
			{
				cot += L.coeff(j, i);
				cot /= 2;
				Vector3d temp_N_vertices = N_vertices.row(i) - N_vertices.row(j);

				//calc Laplacian operator. 
				
				Laplacian_matrix(i, 0) += cot * temp_N_vertices(0);
				Laplacian_matrix(i, 1) += cot * temp_N_vertices(1);
				Laplacian_matrix(i, 2) += cot * temp_N_vertices(2);

				
				sum_energy += cot*(temp_N_vertices).dot(temp_N_vertices);
			}
		}
	
		energy_times_i1 += sum_energy;
		

		Vector3d temp_Laplacian = Laplacian_matrix.row(i);
		Vector3d temp_N_vertices = N_vertices.row(i);

		
		double mod = temp_Laplacian.dot(temp_N_vertices);
	
		Vector3d normal_vector = mod*temp_N_vertices;
		Vector3d tangent_vector;

		tangent_vector = temp_Laplacian - normal_vector;

		Vector3d temp_v = N_vertices.row(i);

		temp_v(0) = N_vertices(i, 0) - tangent_vector(0)*step_size;
		temp_v(1) = N_vertices(i, 1) - tangent_vector(1)*step_size;
		temp_v(2) = N_vertices(i, 2) - tangent_vector(2)*step_size;



		for (int i = 0; i < row_len; i++)
		{
			centroid_vector(0) += N_vertices(i, 0);
			centroid_vector(1) += N_vertices(i, 1);
			centroid_vector(2) += N_vertices(i, 2);
		}
		centroid_vector /= row_len;
		temp_v -= centroid_vector;

		temp_v = temp_v / temp_v.norm();

		N_vertices.row(i) = temp_v;

	}


}
