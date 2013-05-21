// Direct C++ interface example code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

#include <vector>
using namespace std;

#include "voro++.hh"
using namespace voro;

void draw_polygon(FILE *fp,vector<int> &f_vert,vector<double> &v,int j);

int main() {
	unsigned int i,j;
	int id,nx,ny,nz;
	double x,y,z;
	voronoicell_neighbor c;
	vector<int> neigh,f_vert;
	vector<double> v;

	// Create a pre-container class to import the input file and guess the
	// best computational grid size to use.
	pre_container pcon(-3,3,-3,3,0,6,false,false,false);
	pcon.import("pack_six_cube");
	pcon.guess_optimal(nx,ny,nz);

	// Set up the container class and import the particles from the
	// pre-container
	container con(-3,3,-3,3,0,6,nx,ny,nz,false,false,false,8);
	pcon.setup(con);

	// Open the output files
	FILE *fp4=safe_fopen("polygons4_v.pov","w"),
	     *fp5=safe_fopen("polygons5_v.pov","w"),
	     *fp6=safe_fopen("polygons6_v.pov","w");

	// Loop over all particles in the container and compute each Voronoi
	// cell
	c_loop_all cl(con);
	if(cl.start()) do if(con.compute_cell(c,cl)) {
		cl.pos(x,y,z);id=cl.pid();

		// Gather information about the computed Voronoi cell
		c.neighbors(neigh);
		c.face_vertices(f_vert);
		c.vertices(x,y,z,v);

		// Loop over all faces of the Voronoi cell
		for(i=0,j=0;i<neigh.size();i++) {

			// Draw all quadrilaterals, pentagons, and hexagons.
			// Skip if the neighbor information is smaller than
			// this particle's ID, to avoid double counting. This
			// also removes faces that touch the walls, since the
			// neighbor information is set to negative numbers for
			// these cases.
			if(neigh[i]>id) {
				switch(f_vert[j]) {
					case 4: draw_polygon(fp4,f_vert,v,j);
						break;
					case 5: draw_polygon(fp5,f_vert,v,j);
						break;
					case 6: draw_polygon(fp6,f_vert,v,j);
				}
			}

			// Skip to the next entry in the face vertex list
			j+=f_vert[j]+1;
		}
	} while (cl.inc());

	// Close the output files
	fclose(fp4);
	fclose(fp5);
	fclose(fp6);

	// Draw the outline of the domain
	con.draw_domain_pov("polygons_d.pov");
}

void draw_polygon(FILE *fp,vector<int> &f_vert,vector<double> &v,int j) {
	static char s[6][128];
	int k,l,n=f_vert[j];

	// Create POV-Ray vector strings for each of the vertices
	for(k=0;k<n;k++) {
		l=3*f_vert[j+k+1];
		sprintf(s[k],"<%g,%g,%g>",v[l],v[l+1],v[l+2]);
	}

	// Draw the interior of the polygon
	fputs("union{\n",fp);
	for(k=2;k<n;k++) fprintf(fp,"\ttriangle{%s,%s,%s}\n",s[0],s[k-1],s[k]);
	fputs("\ttexture{t1}\n}\n",fp);

	// Draw the outline of the polygon
	fputs("union{\n",fp);
	for(k=0;k<n;k++) {
		l=(k+1)%n;
		fprintf(fp,"\tcylinder{%s,%s,r}\n\tsphere{%s,r}\n",
			s[k],s[l],s[l]);
	}
	fputs("\ttexture{t2}\n}\n",fp);
}

