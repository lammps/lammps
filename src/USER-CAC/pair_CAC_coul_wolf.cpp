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

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include "pair_CAC_coul_wolf.h"
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "update.h"
#include "neigh_list.h"
#include "integrate.h"
#include "respa.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "domain.h"
#include "asa_user.h"

#define MAXNEIGH1  50
#define MAXNEIGH2  10
//#include "math_extra.h"


using namespace LAMMPS_NS;
using namespace MathConst;
using namespace std;


/* ---------------------------------------------------------------------- */

PairCACCoulWolf::PairCACCoulWolf(LAMMPS *lmp) : PairCAC(lmp)
{
  restartinfo = 0;
  nmax = 0;
  outer_neighflag = 0;

  interior_scales = NULL;
  surface_counts = NULL;
  inner_neighbor_coords = NULL;

  inner_neighbor_types = NULL;
  inner_neighbor_charges = NULL;
  surface_counts_max[0] = 0;
  surface_counts_max[1] = 0;
  surface_counts_max[2] = 0;
  surface_counts_max_old[0] = 0;
  surface_counts_max_old[1] = 0;
  surface_counts_max_old[2] = 0;
}

/* ---------------------------------------------------------------------- */

PairCACCoulWolf::~PairCACCoulWolf() {
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);



    memory->destroy(mass_matrix);
    //memory->destroy(force_density_interior);
	memory->destroy(inner_neighbor_coords);

	memory->destroy(inner_neighbor_types);
	memory->destroy(inner_neighbor_charges);

  }
}

/* ---------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
 allocate all arrays
 ------------------------------------------------------------------------- */

void PairCACCoulWolf::allocate()
{
  allocated = 1;
  int n = atom->ntypes;
  max_nodes_per_element = atom->nodes_per_element;
  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");




  memory->create(mass_matrix, max_nodes_per_element, max_nodes_per_element,"pairCAC:mass_matrix");
  memory->create(mass_copy, max_nodes_per_element, max_nodes_per_element,"pairCAC:copy_mass_matrix");
  memory->create(force_column, max_nodes_per_element,3,"pairCAC:force_residue");
  memory->create(current_force_column, max_nodes_per_element,"pairCAC:current_force_residue");
  memory->create(current_nodal_forces, max_nodes_per_element,"pairCAC:current_nodal_force");
  memory->create(pivot, max_nodes_per_element+1,"pairCAC:pivots");
  memory->create(surf_set, 6, 2, "pairCAC:surf_set");
  memory->create(dof_set, 6, 4, "pairCAC:surf_set");
  memory->create(sort_surf_set, 6, 2, "pairCAC:surf_set");
  memory->create(sort_dof_set, 6, 4, "pairCAC:surf_set");
  quadrature_init(2);
}

/* ----------------------------------------------------------------------
global settings
------------------------------------------------------------------------- */
void PairCACCoulWolf::settings(int narg, char **arg) {
	if (narg <2 || narg>3) error->all(FLERR, "Illegal pair_style command");

	//cutmax = force->numeric(FLERR, arg[0]);

	//cut_global_s = force->numeric(FLERR,arg[1]);
	//neighrefresh = force->numeric(FLERR, arg[1]);
	//maxneigh_setting = force->numeric(FLERR, arg[2]);
	
	force->newton_pair = 0;
	
	alf = force->numeric(FLERR, arg[0]);
	cut_global_s = force->numeric(FLERR, arg[1]);
	if (narg == 3) {
		if (strcmp(arg[2], "one") == 0) atom->one_layer_flag=one_layer_flag = 1;
		else error->all(FLERR, "Unexpected argument in CAC/coul/wolf invocation");
	}
	cut_coul = cut_global_s;
	// reset cutoffs that have been explicitly set
	// initialize unit cell vectors

	
}

/* ----------------------------------------------------------------------
 set coeffs for one or more type pairs
 ------------------------------------------------------------------------- */

void PairCACCoulWolf::coeff(int narg, char **arg) {
	if (narg != 2) error->all(FLERR, "Incorrect args for pair coefficients");
	if (!allocated) allocate();

	int ilo, ihi, jlo, jhi;
	force->bounds(FLERR, arg[0], atom->ntypes, ilo, ihi);
	force->bounds(FLERR, arg[1], atom->ntypes, jlo, jhi);

	int count = 0;
	for (int i = ilo; i <= ihi; i++) {
		for (int j = MAX(jlo, i); j <= jhi; j++) {
			setflag[i][j] = 1;
			count++;
		}
	}

	if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
 init for one type pair i,j and corresponding j,i
 ------------------------------------------------------------------------- */

double PairCACCoulWolf::init_one(int i, int j) {


	
	return cut_global_s;
}

/* ---------------------------------------------------------------------- */


void PairCACCoulWolf::init_style()
{
	check_existence_flags();
  if (atom->tag_enable == 0)
    error->all(FLERR,"Pair style CAC_Buck requires atom IDs");
  maxneigh_quad_inner = MAXNEIGH2;
  maxneigh_quad_outer = MAXNEIGH1;
  if (!atom->q_flag)
	  error->all(FLERR, "Pair coul/wolf requires atom attribute q for charges");
  // need a full neighbor list

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->half = 0;
  //neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->CAC = 1;
  cut_coulsq = cut_coul*cut_coul;
  //surface selection array 
  surf_set[0][0] = 1;
  surf_set[0][1] = -1;
  surf_set[1][0] = 1;
  surf_set[1][1] = 1;
  surf_set[2][0] = 2;
  surf_set[2][1] = -1;
  surf_set[3][0] = 2;
  surf_set[3][1] = 1;
  surf_set[4][0] = 3;
  surf_set[4][1] = -1;
  surf_set[5][0] = 3;
  surf_set[5][1] = 1;

  //surface DOF array

  dof_set[0][0] = 0;
  dof_set[0][1] = 3;
  dof_set[0][2] = 4;
  dof_set[0][3] = 7;

  dof_set[1][0] = 1;
  dof_set[1][1] = 2;
  dof_set[1][2] = 5;
  dof_set[1][3] = 6;

  dof_set[2][0] = 0;
  dof_set[2][1] = 1;
  dof_set[2][2] = 4;
  dof_set[2][3] = 5;

  dof_set[3][0] = 2;
  dof_set[3][1] = 3;
  dof_set[3][2] = 6;
  dof_set[3][3] = 7;

  dof_set[4][0] = 0;
  dof_set[4][1] = 1;
  dof_set[4][2] = 2;
  dof_set[4][3] = 3;

  dof_set[5][0] = 4;
  dof_set[5][1] = 5;
  dof_set[5][2] = 6;
  dof_set[5][3] = 7;

  for (int si = 0; si < 6; si++) {
	  sort_dof_set[si][0] = dof_set[si][0];
	  sort_dof_set[si][1] = dof_set[si][1];
	  sort_dof_set[si][2] = dof_set[si][2];
	  sort_dof_set[si][3] = dof_set[si][3];
	  sort_surf_set[si][0] = surf_set[si][0];
	  sort_surf_set[si][1] = surf_set[si][1];
  }
  //minimization algorithm parameters
  //asacg_parm scgParm;
  //asa_parm sasaParm;

  memory->create(cgParm, 1, "pairCAC:cgParm");

  memory->create(asaParm, 1, "pairCAC:asaParm");
  memory->create(Objective, 1, "pairCAC:asaParm");
  // if you want to change parameter value, initialize strucs with default 
  asa_cg_default(cgParm);
  asa_default(asaParm);

  // if you want to change parameters, change them here: 
  cgParm->PrintParms = FALSE;
  cgParm->PrintLevel = 0;

  asaParm->PrintParms = FALSE;
  asaParm->PrintLevel = 0;
  asaParm->PrintFinal = 0;


}







//-----------------------------------------------------------------------


void PairCACCoulWolf::force_densities(int iii, double s, double t, double w, double coefficients,
	double &force_densityx, double &force_densityy, double &force_densityz) {

int internal;

double delx,dely,delz;

double r2inv;
double r6inv;
double shape_func;
double shape_func2;
double boxmap_matrix[3][3];
int neighborflag=0;
int outofbounds=0;
int timestep=update->ntimestep;
double unit_cell_mapped[3];
double scanning_unit_cell[3];
double box_positions[8][3];
double  fpair, forcecoul, factor_coul;
double prefactor;
double r;

double erfcc, erfcd, v_sh, dvdrr, e_self, e_shift, f_shift, qisq;
double *special_coul = force->special_coul;
double qqrd2e = force->qqrd2e;
int *type = atom->type;
double unit_cell[3];
double distancesq;
double current_position[3];
double scan_position[3];
double rcut;
int current_type = poly_counter;
int *element_type = atom->element_type;
double cbox_positions[3];

int flagm;
int neigh_count=0;
int neigh_index=0;
double cds[3];
double maxds=0;
double maxdt=0;
double maxdw=0;
int neighbor_cell_count[3];
int nodes_per_element;
int *nodes_count_list = atom->nodes_per_element_list;	
int neighbor_nodes_per_element;

//equivalent isoparametric cutoff range for a cube of rcut


unit_cell_mapped[0] = 2 / double(current_element_scale[0]);
unit_cell_mapped[1] = 2 / double(current_element_scale[1]);
unit_cell_mapped[2] = 2 / double(current_element_scale[2]);






unit_cell[0] = s;
unit_cell[1] = t;
unit_cell[2] = w;





//scan the surrounding unit cell locations in a cartesian grid
//of isoparametric space until the cutoff is exceeded
//for each grid scan


 scanning_unit_cell[0]=unit_cell[0];
 scanning_unit_cell[1]=unit_cell[1];
 scanning_unit_cell[2]=unit_cell[2];


int distanceflag=0;
    current_position[0]=0;
    current_position[1]=0;
    current_position[2]=0;

	if (!atomic_flag) {
		nodes_per_element = nodes_count_list[current_element_type];
		for (int kkk = 0; kkk < nodes_per_element; kkk++) {
			shape_func = shape_function(unit_cell[0], unit_cell[1], unit_cell[2], 2, kkk + 1);
			current_position[0] += current_nodal_positions[kkk][poly_counter][0] * shape_func;
			current_position[1] += current_nodal_positions[kkk][poly_counter][1] * shape_func;
			current_position[2] += current_nodal_positions[kkk][poly_counter][2] * shape_func;
		}
	}
	else {
		current_position[0] = s;
		current_position[1] = t;
		current_position[2] = w;
	}


	rcut = cut_global_s;
	int origin_type = type_array[poly_counter];
	

			int listtype;
			int listindex;
			int poly_index;
			int scan_type;
			int element_index;
			int *ilist, *jlist, *numneigh, **firstneigh;
			int neigh_max = inner_quad_lists_counts[iii][neigh_quad_counter];
			e_shift = erfc(alf*cut_coul) / cut_coul;
			f_shift = -(e_shift + 2.0*alf / MY_PIS * exp(-alf*alf*cut_coul*cut_coul)) /
				cut_coul;
			ilist = list->ilist;
			numneigh = list->numneigh;
			firstneigh = list->firstneigh;
			jlist = firstneigh[iii];
			double ****nodal_positions = atom->nodal_positions;
			int **node_types = atom->node_types;
			double **node_charges = atom->node_charges;
			double origin_element_charge= node_charges[iii][poly_counter];
			double neighbor_element_charge;
			qisq = origin_element_charge*origin_element_charge;
			e_self = -(e_shift / 2.0 + alf / MY_PIS) * qisq*qqrd2e;
			
			if(neigh_max>local_inner_max){
      memory->grow(inner_neighbor_coords, neigh_max, 3, "Pair_CAC_coul_wolf:neighbor_coords");

			memory->grow(inner_neighbor_types, neigh_max, "Pair_CAC_coul_wolf:neighbor_types");
			memory->grow(inner_neighbor_charges, neigh_max, "Pair_CAC_coul_wolf:neighbor_charges");
	     local_inner_max=neigh_max;
	     }      
			for (int l = 0; l < neigh_max; l++) {
			scanning_unit_cell[0] = inner_quad_lists_ucell[iii][neigh_quad_counter][l][0];
		    scanning_unit_cell[1] = inner_quad_lists_ucell[iii][neigh_quad_counter][l][1];
		    scanning_unit_cell[2] = inner_quad_lists_ucell[iii][neigh_quad_counter][l][2];
		     //listtype = quad_list_container[iii].inner_list2ucell[neigh_quad_counter].cell_indexes[l][0];
		     listindex = inner_quad_lists_index[iii][neigh_quad_counter][l][0];
		    poly_index = inner_quad_lists_index[iii][neigh_quad_counter][l][1];
		    element_index = listindex;
		    element_index &= NEIGHMASK;
		    inner_neighbor_types[l] = node_types[element_index][poly_index];
			inner_neighbor_charges[l] = node_charges[element_index][poly_index];
		    neigh_list_cord(inner_neighbor_coords[l][0], inner_neighbor_coords[l][1], inner_neighbor_coords[l][2],
			  element_index, poly_index, scanning_unit_cell[0], scanning_unit_cell[1], scanning_unit_cell[2]);

			}


			for (int l = 0; l < neigh_max; l++) {

	
				scan_type = inner_neighbor_types[l];
				scan_position[0] = inner_neighbor_coords[l][0];
				scan_position[1] = inner_neighbor_coords[l][1];
				scan_position[2] = inner_neighbor_coords[l][2];
				neighbor_element_charge = inner_neighbor_charges[l];
	
				delx = current_position[0] - scan_position[0];
				dely = current_position[1] - scan_position[1];
				delz = current_position[2] - scan_position[2];
				distancesq = delx*delx + dely*dely + delz*delz;
				if (distancesq < cut_coulsq) {
					r = sqrt(distancesq);
					prefactor = qqrd2e*origin_element_charge*neighbor_element_charge / r;
					erfcc = erfc(alf*r);
					erfcd = exp(-alf*alf*r*r);
					v_sh = (erfcc - e_shift*r) * prefactor;
					dvdrr = (erfcc / distancesq + 2.0*alf / MY_PIS * erfcd / r) + f_shift;
					forcecoul = dvdrr*distancesq*prefactor;
					//if (factor_coul < 1.0) forcecoul -= (1.0 - factor_coul)*prefactor;
					fpair = forcecoul / distancesq;

					force_densityx += delx*fpair;
					force_densityy += dely*fpair;
					force_densityz += delz*fpair;
				}
			}

		





//end of scanning loop


 //induce segfault to debug
 //segv=force_density[133][209];








}
