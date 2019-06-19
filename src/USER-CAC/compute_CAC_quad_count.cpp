/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include <cstring>
#include <cmath>
#include "compute_CAC_quad_count.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "pair.h"
#include "memory.h"
#include "error.h"

#define QUADSCALE  0 //adhoc scaling used to offset the fact that quadrature points
#define QUADSCALE2  0 //adhoc scaling used to offset the fact that quadrature points
#define MAXNEIGH   1000 //used to divide weights so they aren't such large numbers
//for elements are more expensive due to the tendency to have more nodal interpolation
//for their neighbors; optimized for a case with a block of atoms inside a block of FE

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeCACQuadCount::ComputeCACQuadCount(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  quad_count(NULL)
{
  if (narg != 3) error->all(FLERR,"Illegal compute cac/quad/count command");

  peratom_flag = 1;
  size_peratom_cols = 0;

  nmax = 0;
}

/* ---------------------------------------------------------------------- */

ComputeCACQuadCount::~ComputeCACQuadCount()
{
  memory->destroy(quad_count);
}

/* ---------------------------------------------------------------------- */

void ComputeCACQuadCount::init()
{
  int count = 0;
	if(!atom->CAC_flag){
		error->all(FLERR,"compute cac/quad/count requires a CAC atom style");
	}
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"cac/quad/count") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute cac/quad/count");
}

/* ---------------------------------------------------------------------- */

void ComputeCACQuadCount::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow ke array if necessary

  if (atom->nmax > nmax) {
    memory->destroy(quad_count);
    nmax = atom->nmax;
    memory->create(quad_count,nmax,"compute_CAC_quad_count: quad_count");
    vector_atom = quad_count;
  }
  int nlocal = atom->nlocal;
  // compute quadrature counts for each CAC element in the group
  for (int i = 0; i < nlocal; i++) {
	  quad_count[i] = 1;
  }

  int *mask = atom->mask;
  int *element_type = atom->element_type;
  int **neighbor_weights = atom->neighbor_weights;
  double ****nodal_positions = atom->nodal_positions;
  int *poly_count = atom->poly_count;
  int **node_types = atom->node_types;
  int **element_scale = atom->element_scale;
  //int surface_counts[3];
  double interior_scales[3];
  int n1, n2, n3;
  //enable passing quadrature rank that is not 2
  int  quad = atom->quadrature_node_count;
	if(nlocal!=atom->weight_count&&atom->neigh_weight_flag!=0)
	error->one(FLERR, "weight_counts don't match nlocal");
  for (int i = 0; i < nlocal; i++) {
	  if(atom->neigh_weight_flag==0){
			//if(1){
	  current_element_scale = element_scale[i];
	  current_nodal_positions = nodal_positions[i];
	  current_element_type = element_type[i];
	  current_poly_count = poly_count[i];
	  if (current_element_type == 0) quad_count[i]=1;
      if (current_element_type != 0) {
        compute_surface_depths(interior_scales[0], interior_scales[1], interior_scales[2],
				  n1, n2, n3, 1);

        quad_count[i]= quad*quad*quad + 2 * n1*quad*quad + 2 * n2*quad*quad +
				  +2 * n3*quad*quad + 4 * n1*n2*quad + 4 * n3*n2*quad + 4 * n1*n3*quad
				  + 8 * n1*n2*n3;
        quad_count[i] *= current_poly_count;
	  }
	}
	else{
    //quad_count[i]=neighbor_weights[i][0]+QUADSCALE*neighbor_weights[i][1]+QUADSCALE2*neighbor_weights[i][2];
		if(atom->outer_neigh_flag)
		quad_count[i]=neighbor_weights[i][2];
		//quad_count[i]=neighbor_weights[i][0]+neighbor_weights[i][1];
		else
		quad_count[i]=neighbor_weights[i][0]+neighbor_weights[i][1]+QUADSCALE*neighbor_weights[i][1]+QUADSCALE2*neighbor_weights[i][2];
	}
  }
	atom->neigh_weight_flag=0;
}



/* ----------------------------------------------------------------------
surface quadrature location routine used to predict counts
------------------------------------------------------------------------- */

void ComputeCACQuadCount::compute_surface_depths(double &scalex, double &scaley, double &scalez,
	int &countx, int &county, int &countz, int flag) {
	int poly = 0;
	double unit_cell_mapped[3];
	
	double rcut=neighbor->cutneighmax - neighbor->skin;
	//flag determines the current element type and corresponding procedure to calculate parameters for 
	//surface penetration depth to be used when computing force density with influences from neighboring
	//elements

	unit_cell_mapped[0] = 2 / double(current_element_scale[0]);
	unit_cell_mapped[1] = 2 / double(current_element_scale[1]);
	unit_cell_mapped[2] = 2 / double(current_element_scale[2]);
	double ds_x = (current_nodal_positions[0][poly][0] - current_nodal_positions[1][poly][0])*
		(current_nodal_positions[0][poly][0] - current_nodal_positions[1][poly][0]);
	double ds_y = (current_nodal_positions[0][poly][1] - current_nodal_positions[1][poly][1])*
		(current_nodal_positions[0][poly][1] - current_nodal_positions[1][poly][1]);
	double ds_z = (current_nodal_positions[0][poly][2] - current_nodal_positions[1][poly][2])*
		(current_nodal_positions[0][poly][2] - current_nodal_positions[1][poly][2]);
	double ds_surf = 2 * rcut / sqrt(ds_x + ds_y + ds_z);
	ds_surf = unit_cell_mapped[0] * (int)(ds_surf / unit_cell_mapped[0]) + unit_cell_mapped[0];

	double dt_x = (current_nodal_positions[0][poly][0] - current_nodal_positions[3][poly][0])*
		(current_nodal_positions[0][poly][0] - current_nodal_positions[3][poly][0]);
	double dt_y = (current_nodal_positions[0][poly][1] - current_nodal_positions[3][poly][1])*
		(current_nodal_positions[0][poly][1] - current_nodal_positions[3][poly][1]);
	double dt_z = (current_nodal_positions[0][poly][2] - current_nodal_positions[3][poly][2])*
		(current_nodal_positions[0][poly][2] - current_nodal_positions[3][poly][2]);

	double dt_surf = 2 * rcut / sqrt(dt_x + dt_y + dt_z);
	dt_surf = unit_cell_mapped[1] * (int)(dt_surf / unit_cell_mapped[1]) + unit_cell_mapped[1];

	double dw_x = (current_nodal_positions[0][poly][0] - current_nodal_positions[4][poly][0])*
		(current_nodal_positions[0][poly][0] - current_nodal_positions[4][poly][0]);
	double dw_y = (current_nodal_positions[0][poly][1] - current_nodal_positions[4][poly][1])*
		(current_nodal_positions[0][poly][1] - current_nodal_positions[3][poly][1]);
	double dw_z = (current_nodal_positions[0][poly][2] - current_nodal_positions[4][poly][2])*
		(current_nodal_positions[0][poly][2] - current_nodal_positions[4][poly][2]);

	double dw_surf = 2 * rcut / sqrt(dw_x + dw_y + dw_z);
	dw_surf = unit_cell_mapped[2] * (int)(dw_surf / unit_cell_mapped[2]) + unit_cell_mapped[2];
	if (ds_surf > 1) {
		ds_surf = 1;
	}
	if (dt_surf > 1) {
		dt_surf = 1;
	}
	if (dw_surf > 1) {
		dw_surf = 1;
	}
    if (atom->one_layer_flag) {
		ds_surf = unit_cell_mapped[0];
		dt_surf = unit_cell_mapped[1];
		dw_surf = unit_cell_mapped[2];
	}

	scalex = 1 - ds_surf;
	scaley = 1 - dt_surf;
	scalez = 1 - dw_surf;

	countx = (int)(ds_surf / unit_cell_mapped[0]);
	county = (int)(dt_surf / unit_cell_mapped[1]);
	countz = (int)(dw_surf / unit_cell_mapped[2]);
  if(countx==0) countx=1;
	if(county==0) county=1;
	if(countz==0) countz=1;



}



/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeCACQuadCount::memory_usage()
{
  double bytes = nmax * sizeof(double);
  return bytes;
}
