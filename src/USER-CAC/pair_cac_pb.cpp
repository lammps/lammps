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
#include "pair_cac_pb.h"
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

#define MAXNEIGHIN  100
#define MAXNEIGHOUT  200
#define EXPAND 10
//#include "math_extra.h"
using namespace LAMMPS_NS;
using namespace MathConst;


/* ---------------------------------------------------------------------- */

PairCACPb::PairCACPb(LAMMPS *lmp) : PairCAC(lmp)
{
  restartinfo = 0;
  nmax = 0;
  outer_neighflag = 0;
  inner_neighbor_coords = NULL;
  inner_neighbor_types = NULL;
  inner_neighbor_charges = NULL;
}

/* ---------------------------------------------------------------------- */

PairCACPb::~PairCACPb() {
  if (allocated) {
  memory->destroy(setflag);
  memory->destroy(cutsq);
	memory->destroy(cut);
	memory->destroy(a);
	memory->destroy(rho);
	memory->destroy(c);
	memory->destroy(rhoinv);
	memory->destroy(buck1);
	memory->destroy(buck2);
	memory->destroy(offset);
	memory->destroy(inner_neighbor_coords);
	memory->destroy(inner_neighbor_types);
	memory->destroy(inner_neighbor_charges);
  }
}

/* ---------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
 allocate all arrays
 ------------------------------------------------------------------------- */

void PairCACPb::allocate()
{
  allocated = 1;
  int n = atom->ntypes;
  max_nodes_per_element = atom->nodes_per_element;
  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(cut, n + 1, n + 1, "pair:cut_lj");
  memory->create(a, n + 1, n + 1, "pair:a");
  memory->create(rho, n + 1, n + 1, "pair:rho");
  memory->create(c, n + 1, n + 1, "pair:c");
  memory->create(rhoinv, n + 1, n + 1, "pair:rhoinv");
  memory->create(buck1, n + 1, n + 1, "pair:buck1");
  memory->create(buck2, n + 1, n + 1, "pair:buck2");
  memory->create(offset, n + 1, n + 1, "pair:offset");


  memory->create(mass_matrix,max_nodes_per_element, max_nodes_per_element,"pairCAC:mass_matrix");
  memory->create(mass_copy, max_nodes_per_element, max_nodes_per_element,"pairCAC:copy_mass_matrix");
  memory->create(force_column, max_nodes_per_element,3,"pairCAC:force_residue");
  memory->create(current_force_column, max_nodes_per_element,"pairCAC:current_force_residue");
  memory->create(current_nodal_forces, max_nodes_per_element,"pairCAC:current_nodal_force");
  memory->create(pivot, max_nodes_per_element+1,"pairCAC:pivots");
  quadrature_init(2);
}

/* ----------------------------------------------------------------------
global settings
------------------------------------------------------------------------- */
void PairCACPb::settings(int narg, char **arg) {
	if (narg <3 || narg>4) error->all(FLERR, "Illegal pair_style command");

	force->newton_pair = 0;
	cut_buck = force->numeric(FLERR, arg[0]);
	alf = force->numeric(FLERR, arg[1]);
	cut_coul = force->numeric(FLERR, arg[2]);
	 if (narg == 4) {
		if (strcmp(arg[3], "one") == 0) atom->one_layer_flag=one_layer_flag = 1;
		else error->all(FLERR, "Unexpected argument in CAC/Pb pair style invocation");

	}
	cut_global_s = cut_buck;
	if(cut_coul>cut_buck)
	 cut_global_s=cut_coul;

	if (allocated) {
		int i, j;
		for (i = 1; i <= atom->ntypes; i++)
			for (j = i; j <= atom->ntypes; j++)
				if (setflag[i][j]) cut[i][j] = cut_global_s;
	}
}

/* ----------------------------------------------------------------------
 set coeffs for one or more type pairs
 ------------------------------------------------------------------------- */

void PairCACPb::coeff(int narg, char **arg) {
	if (narg < 5 || narg > 6)
		error->all(FLERR, "Incorrect args for pair coefficients");
	if (!allocated) allocate();
	setflag[1][1] = 1;
	int ilo, ihi, jlo, jhi;
	force->bounds(FLERR, arg[0], atom->ntypes, ilo, ihi);
	force->bounds(FLERR, arg[1], atom->ntypes, jlo, jhi);

	double a_one = force->numeric(FLERR, arg[2]);
	double rho_one = force->numeric(FLERR, arg[3]);
	if (rho_one <= 0) error->all(FLERR, "Incorrect args for pair coefficients");
	double c_one = force->numeric(FLERR, arg[4]);

	double cut_one = cut_global_s;
	if (narg == 6) cut_one = force->numeric(FLERR, arg[5]);

	int count = 0;
	for (int i = ilo; i <= ihi; i++) {
		for (int j = MAX(jlo, i); j <= jhi; j++) {
			a[i][j] = a_one;
			rho[i][j] = rho_one;
			c[i][j] = c_one;
			cut[i][j] = cut_one;
			setflag[i][j] = 1;
			count++;
		}
	}

	if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
 init for one type pair i,j and corresponding j,i
 ------------------------------------------------------------------------- */

double PairCACPb::init_one(int i, int j) {

	//if (setflag[i][j] == 0) error->all(FLERR, "All pair coeffs are not set");

	rhoinv[i][j] = 1.0 / rho[i][j];
	buck1[i][j] = a[i][j] / rho[i][j];
	buck2[i][j] = 6.0*c[i][j];

	if (offset_flag && (cut[i][j] > 0.0)) {
		double rexp = exp(-cut[i][j] / rho[i][j]);
		offset[i][j] = a[i][j] * rexp - c[i][j] / pow(cut[i][j], 6.0);
	}
	else offset[i][j] = 0.0;

	a[j][i] = a[i][j];
	c[j][i] = c[i][j];
	rhoinv[j][i] = rhoinv[i][j];
	buck1[j][i] = buck1[i][j];
	buck2[j][i] = buck2[i][j];
	offset[j][i] = offset[i][j];

	// compute I,J contribution to long-range tail correction
	// count total # of atoms of type I and J via Allreduce

	if (tail_flag) {
		int *type = atom->type;
		int nlocal = atom->nlocal;

		double count[2], all[2];
		count[0] = count[1] = 0.0;
		for (int k = 0; k < nlocal; k++) {
			if (type[k] == i) count[0] += 1.0;
			if (type[k] == j) count[1] += 1.0;
		}
		MPI_Allreduce(count, all, 2, MPI_DOUBLE, MPI_SUM, world);

		double rho1 = rho[i][j];
		double rho2 = rho1*rho1;
		double rho3 = rho2*rho1;
		double rc = cut[i][j];
		double rc2 = rc*rc;
		double rc3 = rc2*rc;
		etail_ij = 2.0*MY_PI*all[0] * all[1] *
			(a[i][j] * exp(-rc / rho1)*rho1*(rc2 + 2.0*rho1*rc + 2.0*rho2) -
				c[i][j] / (3.0*rc3));
		ptail_ij = (-1 / 3.0)*2.0*MY_PI*all[0] * all[1] *
			(-a[i][j] * exp(-rc / rho1)*
			(rc3 + 3.0*rho1*rc2 + 6.0*rho2*rc + 6.0*rho3) + 2.0*c[i][j] / rc3);
	}
		

	return cut_global_s;
}

/* ---------------------------------------------------------------------- */


void PairCACPb::init_style()
{
	check_existence_flags();
  if (atom->tag_enable == 0)
    error->all(FLERR,"Pair style CAC_Buck requires atom IDs");
  if (!atom->q_flag)
	  error->all(FLERR, "Pair coul/wolf requires atom attribute q for charges");
  atom->max_neigh_inner_init = maxneigh_quad_inner = MAXNEIGHIN;
  atom->max_neigh_outer_init = maxneigh_quad_outer = MAXNEIGHOUT;
  // need a full neighbor list

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->half = 0;
  //neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->cac = 1;
}

//-----------------------------------------------------------------------

void PairCACPb::force_densities(int iii, double s, double t, double w, double coefficients,
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
double *special_lj = force->special_lj;
double  forcebuck, factor_lj, fpair;
double  fpair2, forcecoul, factor_coul;
double prefactor;
double r, rexp;
int *type = atom->type;
double unit_cell[3];
double distancesq;
double current_position[3];
double scan_position[3];
double rcut;
int current_type = poly_counter;
double erfcc, erfcd, v_sh, dvdrr, e_self, e_shift, f_shift, qisq;
double *special_coul = force->special_coul;
double qqrd2e = force->qqrd2e;
double cbox_positions[3];

int flagm;
int nodes_per_element;
int *nodes_count_list = atom->nodes_per_element_list;	
int neigh_count=0;
int neigh_index=0;
double cds[3];
double maxds=0;
double maxdt=0;
double maxdw=0;
int neighbor_cell_count[3];

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
			current_position[0] += current_nodal_positions[kkk][0] * shape_func;
			current_position[1] += current_nodal_positions[kkk][1] * shape_func;
			current_position[2] += current_nodal_positions[kkk][2] * shape_func;
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
	int neigh_max = inner_quad_lists_counts[pqi];
	e_shift = erfc(alf*cut_coul) / cut_coul;
	f_shift = -(e_shift + 2.0*alf / MY_PIS * exp(-alf*alf*cut_coul*cut_coul)) /
	cut_coul;
	double ****nodal_positions = atom->nodal_positions;
	int **node_types = atom->node_types;
	double **node_charges = atom->node_charges;
	double origin_element_charge = node_charges[iii][poly_counter];
	double neighbor_element_charge;
  int **inner_quad_indices = inner_quad_lists_index[pqi];
	qisq = origin_element_charge*origin_element_charge;

	if(neigh_max>local_inner_max){
	memory->grow(inner_neighbor_coords, neigh_max+EXPAND, 3, "Pair_CAC_pb:inner_neighbor_coords");
	memory->grow(inner_neighbor_types, neigh_max+EXPAND, "Pair_CAC_pb:inner_neighbor_types");
	memory->grow(inner_neighbor_charges, neigh_max+EXPAND, "Pair_CAC_pb:inner_neighbor_types");
	local_inner_max=neigh_max+EXPAND;
	}
	
	for (int l = 0; l < neigh_max; l++) {
	element_index = inner_quad_indices[l][0];
  poly_index = inner_quad_indices[l][1];
	inner_neighbor_types[l] = node_types[element_index][poly_index];
	inner_neighbor_charges[l] = node_charges[element_index][poly_index];
	}

      //interpolate virtual atom coordinates from shape functions corresponding to unit cells
      interpolation(iii);
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
				if (distancesq < cut_global_s*cut_global_s) {
					fpair = 0;
					fpair2 = 0;
					if (distancesq < cut_coul*cut_coul) {
						r = sqrt(distancesq);
						prefactor = qqrd2e*origin_element_charge*neighbor_element_charge / r;
						erfcc = erfc(alf*r);
						erfcd = exp(-alf*alf*r*r);
						v_sh = (erfcc - e_shift*r) * prefactor;
						dvdrr = (erfcc / distancesq + 2.0*alf / MY_PIS * erfcd / r) + f_shift;
						forcecoul = dvdrr*distancesq*prefactor;
						//if (factor_coul < 1.0) forcecoul -= (1.0 - factor_coul)*prefactor;
						fpair = forcecoul / distancesq;
					}
					if (distancesq < cut_buck*cut_buck) {
						if ((scan_type==1&&origin_type!=1)|| (scan_type != 1 && origin_type == 1)|| (scan_type != 1 && origin_type != 1)) {
							//buckingham force
							r2inv = 1.0 / distancesq;
							r6inv = r2inv*r2inv*r2inv;
							rexp = exp(-r*rhoinv[origin_type][scan_type]);
							forcebuck = buck1[origin_type][scan_type] * r*rexp - buck2[origin_type][scan_type] * r6inv;
							fpair2 = forcebuck*r2inv;
						}
					}
					//compute total pair force
					force_densityx += delx*(fpair + fpair2);
					force_densityy += dely*(fpair + fpair2);
					force_densityz += delz*(fpair + fpair2);
				}
				if(atom->CAC_virial){
		    virial_density[0] += 0.5*delx*delx*fpair;
		    virial_density[1] += 0.5*dely*dely*fpair;
		    virial_density[2] += 0.5*delz*delz*fpair;
		    virial_density[3] += 0.5*delx*dely*fpair;
		    virial_density[4] += 0.5*delx*delz*fpair;
		    virial_density[5] += 0.5*dely*delz*fpair;
		        }
        if (quad_eflag) 
				quadrature_energy += (a[origin_type][scan_type]*rexp - c[origin_type][scan_type]*r6inv -
                      offset[origin_type][scan_type])/2+v_sh/2;
			}


}
