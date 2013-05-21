// Simple cell statistics demonstration code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

#include "voro++.hh"
using namespace voro;

// This function returns a random floating point number between 0 and 1
double rnd() {return double(rand())/RAND_MAX;}

int main() {
	double x,y,z;
	voronoicell v;

	// Initialize the Voronoi cell to be a cube of side length 2, centered
	// on the origin
	v.init(-1,1,-1,1,-1,1);

	// Remove one edge of the cell with a single plane cut
	v.plane(1,1,0,2);

	// Output the Voronoi cell to a file in gnuplot format
	v.draw_gnuplot(0,0,0,"simple_cell.gnu");

	// Output vertex-based statistics
	printf("Total vertices      : %d\n",v.p);
	printf("Vertex positions    : ");v.output_vertices();puts("");
	printf("Vertex orders       : ");v.output_vertex_orders();puts("");
	printf("Max rad. sq. vertex : %g\n\n",0.25*v.max_radius_squared());

	// Output edge-based statistics
	printf("Total edges         : %d\n",v.number_of_edges());
	printf("Total edge distance : %g\n",v.total_edge_distance());
	printf("Face perimeters     : ");v.output_face_perimeters();puts("\n");

	// Output face-based statistics
	printf("Total faces         : %d\n",v.number_of_faces());
	printf("Surface area        : %g\n",v.surface_area());
	printf("Face freq. table    : ");v.output_face_freq_table();puts("");
	printf("Face orders         : ");v.output_face_orders();puts("");
	printf("Face areas          : ");v.output_face_areas();puts("");
	printf("Face normals        : ");v.output_normals();puts("");
	printf("Face vertices       : ");v.output_face_vertices();puts("\n");

	// Output volume-based statistics
	v.centroid(x,y,z);
	printf("Volume              : %g\n"
	       "Centroid vector     : (%g,%g,%g)\n",v.volume(),x,y,z);

}
