// Cell cutting region example code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

#include "voro++.hh"
using namespace voro;

const double pi=3.1415926535897932384626433832795;

// This constant sets the tolerance in the bisection search algorithm
const double tolwidth=1e-7;

// This constant determines the density of points to test
const double theta_step=pi/200;

int main() {
	double x,y,z,r,rmin,rmax;
	double theta,phi,phi_step;
	voronoicell v;
	FILE *fp=safe_fopen("cell_cut_region.gnu","w");

	// Initialize the Voronoi cell to be an octahedron and make a single
	// plane cut to add some variation
	v.init_octahedron(1);
	v.plane(0.4,0.3,1,0.1);

	// Output the cell in gnuplot format
	v.draw_gnuplot(0,0,0,"cell.gnu");

	// Now test over direction vectors from the center of the sphere. For
	// each vector, carry out a search to find the maximum distance along
	// that vector that a plane will intersect with cell, and save it to
	// the output file.
	for(theta=theta_step*0.5;theta<pi;theta+=theta_step) {
		phi_step=2*pi/(int(2*pi*sin(theta)/theta_step)+1);
		for(phi=phi_step*0.5;phi<2*pi;phi+=phi_step) {

			// Calculate a direction to look along
			x=sin(theta)*cos(phi);
			y=sin(theta)*sin(phi);
			z=cos(theta);

			// Now carry out a bisection search. Here, we initialize
			// a minimum and a maximum guess for the distance
			// along this vector. Keep multiplying rmax by two until
			// the plane no longer makes a cut.
			rmin=0;rmax=1;
			while (v.plane_intersects(x,y,z,rmax)) rmax*=2;

			// We now know that the distance is somewhere between
			// rmin and rmax. Test the point halfway in-between
			// these two. If it intersects, then move rmin to this
			// point; otherwise, move rmax there. At each stage the
			// bounding interval is divided by two. Exit when the
			// width of the interval is smaller than the tolerance.
			while (rmax-rmin>tolwidth) {
				r=(rmax+rmin)*0.5;
				if (v.plane_intersects(x,y,z,r)) rmin=r;
				else rmax=r;
			}

			// Output this point to file
			r=(rmax+rmin)*0.5;
			x*=r;y*=r;z*=r;
			fprintf(fp,"%g %g %g\n",x,y,z);
		}
	}

	fclose(fp);
}
