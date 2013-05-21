// Timing test example code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

#include <ctime>
using namespace std;

#include "voro++.hh"
using namespace voro;

// Set up constants for the container geometry
const double x_min=-1,x_max=1;
const double y_min=-1,y_max=1;
const double z_min=-1,z_max=1;

// Set up the number of blocks that the container is divided into. If the
// preprocessor variable NNN hasn't been passed to the code, then initialize it
// to a good value. Otherwise, use the value that has been passed.
#ifndef NNN
#define NNN 26
#endif
const int n_x=NNN,n_y=NNN,n_z=NNN;

// Set the number of particles that are going to be randomly introduced
const int particles=100000;

// This function returns a random double between 0 and 1
double rnd() {return double(rand())/RAND_MAX;}

int main() {
	clock_t start,end;
	int i;double x,y,z;

	// Create a container with the geometry given above, and make it
	// periodic in each of the three coordinates. Allocate space for eight
	// particles within each computational block.
	container con(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,
			true,true,true,8);

	//Randomly add particles into the container
	for(i=0;i<particles;i++) {
		x=x_min+rnd()*(x_max-x_min);
		y=y_min+rnd()*(y_max-y_min);
		z=z_min+rnd()*(z_max-z_min);
		con.put(i,x,y,z);
	}

	// Store the initial clock time
	start=clock();

	// Carry out a dummy computation of all cells in the entire container
	con.compute_all_cells();

	// Calculate the elapsed time and print it
	end=clock();
	double runtime=double(end-start)/CLOCKS_PER_SEC;
	printf("%g\n",runtime);
}
