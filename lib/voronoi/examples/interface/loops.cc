// Example code demonstrating the loop classes
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

#include "voro++.hh"
using namespace voro;

// Constants determining the configuration of the tori
const double dis=1.25,mjrad=2.5,mirad=0.95,trad=mjrad+mirad;

// Set the number of particles that are going to be randomly introduced
const int particles=100000;

// This function returns a random double between 0 and 1
double rnd() {return double(rand())/RAND_MAX;}

int main() {
	int i;
	double x,y,z,r;
	voronoicell c;

	// Create a container as a non-periodic 10 by 10 by 10 box
	container con(-5,5,-5,5,-5,5,26,26,26,false,false,false,8);
	particle_order po;

	// Randomly add particles into the container
	for(i=0;i<particles;i++) {
		x=10*rnd()-5;
		y=10*rnd()-5;
		z=10*rnd()-5;

		// If the particle lies within the first torus, store it in the
		// ordering class when adding to the container
		r=sqrt((x-dis)*(x-dis)+y*y);
		if((r-mjrad)*(r-mjrad)+z*z<mirad) con.put(po,i,x,y,z);
		else con.put(i,x,y,z);
	}

	// Compute Voronoi cells for the first torus. Here, the points
	// previously stored in the ordering class are looped over.
	FILE *f1=safe_fopen("loops1_m.pov","w");
	FILE *f2=safe_fopen("loops1_v.pov","w");
	c_loop_order clo(con,po);
	if(clo.start()) do if(con.compute_cell(c,clo)) {

		// Get the position of the current particle under consideration
		clo.pos(x,y,z);

		// Save a POV-Ray mesh to one file and a cylinder/sphere
		// representation to the other file
		c.draw_pov_mesh(x,y,z,f1);
		c.draw_pov(x,y,z,f2);
	} while (clo.inc());
	fclose(f1);
	fclose(f2);

	// Compute Voronoi cells for the second torus. Here, the subset loop is
	// used to search over the blocks overlapping the torus, and then each
	// particle is individually tested.
	f1=safe_fopen("loops2_m.pov","w");
	f2=safe_fopen("loops2_v.pov","w");
	c_loop_subset cls(con);
	cls.setup_box(-dis-trad,-dis+trad,-mirad,mirad,-trad,trad,false);
	if(cls.start()) do {

		// Get the position of the current particle under consideration
		cls.pos(x,y,z);

		// Test whether this point is within the torus, and if so,
		// compute and save the Voronoi cell
		r=sqrt((x+dis)*(x+dis)+z*z);
		if((r-mjrad)*(r-mjrad)+y*y<mirad&&con.compute_cell(c,cls)) {
			c.draw_pov_mesh(x,y,z,f1);
			c.draw_pov(x,y,z,f2);
		}
	} while (cls.inc());
	fclose(f1);
	fclose(f2);
}
