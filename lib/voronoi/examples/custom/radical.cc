// Radical Voronoi tessellation example code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

#include "voro++.hh"
using namespace voro;

// Set up constants for the container geometry
const double x_min=-3,x_max=3;
const double y_min=-3,y_max=3;
const double z_min=0,z_max=6;

// Set up the number of blocks that the container is divided
// into.
const int n_x=3,n_y=3,n_z=3;

int main() {

	// Create a container with the geometry given above, and make it
	// non-periodic in each of the three coordinates. Allocate space for
	// eight particles within each computational block. Import
	// the monodisperse test packing and output the Voronoi
	// tessellation in gnuplot and POV-Ray formats.
	container con(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,
			false,false,false,8);
	con.import("pack_six_cube");
	con.draw_cells_gnuplot("pack_six_cube.gnu");
	con.draw_cells_pov("pack_six_cube_v.pov");
	con.draw_particles_pov("pack_six_cube_p.pov");

	// Create a container for polydisperse particles using the same
	// geometry as above. Import the polydisperse test packing and
	// output the Voronoi radical tessellation in gnuplot and POV-Ray
	// formats.
	container_poly conp(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,
			false,false,false,8);
	conp.import("pack_six_cube_poly");
	conp.draw_cells_gnuplot("pack_six_cube_poly.gnu");
	conp.draw_cells_pov("pack_six_cube_poly_v.pov");
	conp.draw_particles_pov("pack_six_cube_poly_p.pov");
}
