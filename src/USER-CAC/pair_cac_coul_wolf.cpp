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

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "pair_cac_coul_wolf.h"
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

#define MAXNEIGHOUT  50
#define MAXNEIGHIN  10
#define EXPAND 10
using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

PairCACCoulWolf::PairCACCoulWolf(LAMMPS *lmp) : PairCAC(lmp)
{
  restartinfo = 0;
  nmax = 0;
  outer_neighflag = 0;
  flux_enable = 1;
}

/* ---------------------------------------------------------------------- */

PairCACCoulWolf::~PairCACCoulWolf() {
  if (allocated) {
  memory->destroy(setflag);
  memory->destroy(cutsq);
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
  quadrature_init(2);
}

/* ----------------------------------------------------------------------
global settings
------------------------------------------------------------------------- */
void PairCACCoulWolf::settings(int narg, char **arg) {
  if (narg <2 || narg>3) error->all(FLERR, "Illegal pair_style command");

  force->newton_pair = 0;

  alf = utils::numeric(FLERR, arg[0],false,lmp);
  cut_global_s = utils::numeric(FLERR, arg[1],false,lmp);
  if (narg == 3) {
    if (strcmp(arg[2], "one") == 0) atom->one_layer_flag=one_layer_flag = 1;
    else error->all(FLERR, "Unexpected argument in cac/coul/wolf invocation");
  }
  cut_coul = cut_global_s;
}

/* ----------------------------------------------------------------------
 set coeffs for one or more type pairs
 ------------------------------------------------------------------------- */

void PairCACCoulWolf::coeff(int narg, char **arg) {
  if (narg != 2) error->all(FLERR, "Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo, ihi, jlo, jhi;
  utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
  utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);

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
  PairCAC::init_style();
  atom->max_neigh_inner_init = maxneigh_quad_inner = MAXNEIGHIN;
  atom->max_neigh_outer_init = maxneigh_quad_outer = MAXNEIGHOUT;
  if (atom->tag_enable == 0)
    error->all(FLERR,"Pair style cac/coul/wolf requires atom IDs");
  atom->max_neigh_inner_init = maxneigh_quad_inner = MAXNEIGHIN;
  atom->max_neigh_outer_init = maxneigh_quad_outer = MAXNEIGHOUT;
  if (!atom->q_flag)
    error->all(FLERR, "Pair coul/wolf requires atom attribute q for charges");
  cut_coulsq = cut_coul*cut_coul;
  e_shift = erfc(alf*cut_coul) / cut_coul;
  f_shift = -(e_shift + 2.0*alf / MY_PIS * exp(-alf*alf*cut_coul*cut_coul)) /
  cut_coul;
}

//-----------------------------------------------------------------------


void PairCACCoulWolf::force_densities(int iii, double s, double t, double w, double coefficients,
  double &force_densityx, double &force_densityy, double &force_densityz) {

  double delx,dely,delz;

  double r2inv;
  double r6inv;
  double shape_func;
  double shape_func2;
  double fpair;
  double prefactor;
  double r;

  double e_self, qisq;
  double *special_coul = force->special_coul;
  double qqrd2e = force->qqrd2e;
  int *type = atom->type;
  double distancesq;
  double scan_position[3];
  double rcut;
  int current_type = poly_counter;
  int *element_type = atom->element_type;

  int nodes_per_element;
  int *nodes_count_list = atom->nodes_per_element_list;
  int neighbor_nodes_per_element;

  rcut = cut_global_s;
  int origin_type = type_array[poly_counter];

  int listtype;
  int listindex;
  int poly_index;
  int scan_type;
  int element_index;
  int *ilist, *numneigh, **firstneigh;
  int neigh_max = inner_quad_lists_counts[pqi];
  double ****nodal_positions = atom->nodal_positions;
  int **node_types = atom->node_types;
  double **node_charges = atom->node_charges;
  double origin_element_charge= node_charges[iii][poly_counter];
  double neighbor_element_charge;
  int **inner_quad_indices = inner_quad_lists_index[pqi];
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  qisq = origin_element_charge*origin_element_charge;
  e_self = -(e_shift / 2.0 + alf / MY_PIS) * qisq*qqrd2e;
  quadrature_energy += e_self;

  //allocate arrays that store neighbor information around just this quadrature point
  allocate_quad_memory();
  //set virtual neighbor types, etc.
  init_quad_arrays();
  //interpolate virtual atom coordinates from shape functions corresponding to unit cells
  interpolation(iii,s,t,w);
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
      factor_coul = special_coul[sbmask(inner_quad_indices[l][0])];
      fpair = pair_interaction_q(distancesq, origin_type, scan_type,
        origin_element_charge, neighbor_element_charge);

      force_densityx += delx*fpair;
      force_densityy += dely*fpair;
      force_densityz += delz*fpair;
      if(atom->CAC_virial){
      virial_density[0] += 0.5*delx*delx*fpair;
      virial_density[1] += 0.5*dely*dely*fpair;
      virial_density[2] += 0.5*delz*delz*fpair;
      virial_density[3] += 0.5*delx*dely*fpair;
      virial_density[4] += 0.5*delx*delz*fpair;
      virial_density[5] += 0.5*dely*delz*fpair;
      }
      if (quad_eflag)
        quadrature_energy += v_sh/2;
      //cac flux contribution due to current quadrature point and neighbor pair interactions
      if(quad_flux_flag){
        current_quad_flux(l,delx*fpair,dely*fpair,delz*fpair);
      }
    }
  }
  //end of force density loop
  
  //additional cac flux contributions due to neighbors interacting with neighbors
  //  in the vicinity of this quadrature point
  if (quad_flux_flag) {
    //compute_intersections();
    quad_neigh_flux();
  }
}

/* ---------------------------------------------------------------------- */

double PairCACCoulWolf::pair_interaction_q(double distancesq, int itype, int jtype
                                          , double qi, double qj)
{
  double fpair, forcecoul, prefactor, dvdrr;
  double r, erfcc, erfcd;
  double qqrd2e = force->qqrd2e;

  r = sqrt(distancesq);
  prefactor = qqrd2e*qi*qj / r;
  erfcc = erfc(alf*r);
  erfcd = exp(-alf*alf*r*r);
  v_sh = (erfcc - e_shift*r) * prefactor;
  dvdrr = (erfcc / distancesq + 2.0*alf / MY_PIS * erfcd / r) + f_shift;
  forcecoul = dvdrr*distancesq*prefactor;
  if (factor_coul < 1.0) forcecoul -= (1.0 - factor_coul)*prefactor;
  fpair = forcecoul / distancesq;
  if(quad_eflag)
    if (factor_coul < 1.0) v_sh -= (1.0-factor_coul)*prefactor;
  return fpair;
}
