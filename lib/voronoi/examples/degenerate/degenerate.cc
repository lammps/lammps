// Degenerate Voronoi cell example code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

#include "voro++.hh"
using namespace voro;

const double pi=3.1415926535897932384626433832795;

// The number of planes to be cut around each coordinate axis
const int n=32;
const double step=2*pi/n;

// The angle (in radians) of the cutting planes from horizontal
const double theta=pi/4-0.25;

int main() {
	double x,y,z,phi;
	voronoicell v;

	// Initialize the Voronoi cell to be a cube of side length 2, centered
	// on the origin
	v.init(-1,1,-1,1,-1,1);

	// Plane cutting
	for(phi=0;phi<2*pi-0.5*step;phi+=step) {
		x=cos(theta);y=cos(phi)*sin(theta);z=sin(phi)*sin(theta);
		v.plane(x,y,z,1);
		v.plane(-x,y,z,1);
		v.plane(y,x,z,1);
		v.plane(y,-x,z,1);
		v.plane(y,z,x,1);
		v.plane(y,z,-x,1);
	}

	// Check that the relation table is correct, and that there are no
	// duplicate edges
	v.check_relations();
	v.check_duplicates();

	// Output the Voronoi cell to a file in Gnuplot format
	v.draw_gnuplot(0,0,0,"degenerate.gnu");

	// Output the Voronoi cell to a file in POV-Ray format
	v.draw_pov(0,0,0,"degenerate_v.pov");
}
