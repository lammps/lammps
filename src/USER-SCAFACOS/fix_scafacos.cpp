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

/* ----------------------------------------------------------------------
   Contributing author: Rene Halver (JSC)
------------------------------------------------------------------------- */

#include <cstdio>
#include <cstring>
#include <cstdlib>
#include "fix_scafacos.h"
#include "atom.h"
#include "comm.h"
#include "update.h"
#include "neighbor.h"
#include "domain.h"
#include "force.h"
#include "pair.h"
#include "modify.h"
#include "compute.h"
#include "memory.h"
#include "error.h"

// ScaFaCoS library
#include <sstream>
#include <string>
#include "fcs.h"
#include "universe.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#define INVOKED_PERATOM 8

/* ---------------------------------------------------------------------- */

FixScafacos::FixScafacos(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  // form of fix string:
  // fix <#> <group> scafacos <method> [ tolerance <type> <value> ]


  if (narg < 4) error->all(FLERR,"Illegal fix scafacos command");

  // read method from fix input
  method = arg[3];

  int arg_index = 4;

  tolerance_set = false;

  while(arg_index < narg)
  {
    // check if tolerance option is set
    if (strcmp(arg[arg_index],"tolerance") == 0)
    {
      tolerance_set = true;
      ++arg_index; 
      if (strcmp(arg[arg_index],"energy") == 0)
        tolerance_type = FCS_TOLERANCE_TYPE_ENERGY;     
      else if (strcmp(arg[arg_index],"energy_rel") == 0)
        tolerance_type = FCS_TOLERANCE_TYPE_ENERGY_REL;     
      else if (strcmp(arg[arg_index],"field") == 0)
        tolerance_type = FCS_TOLERANCE_TYPE_FIELD;     
      else if (strcmp(arg[arg_index],"field_rel") == 0)
        tolerance_type = FCS_TOLERANCE_TYPE_FIELD_REL;     
      else if (strcmp(arg[arg_index],"potential") == 0)
        tolerance_type = FCS_TOLERANCE_TYPE_POTENTIAL;     
      else if (strcmp(arg[arg_index],"potential_rel") == 0)
        tolerance_type = FCS_TOLERANCE_TYPE_POTENTIAL_REL;     
      else
        error->all(FLERR,"Illegal fix scafacos command");
      ++arg_index;
      tolerance_value = atof(arg[arg_index]);
      ++arg_index; 
    }
    else
      error->all(FLERR,"Illegal fix scafacos command");
  }
}

/* ---------------------------------------------------------------------- */

FixScafacos::~FixScafacos()
{
  memory->destroy(pot);
  memory->destroy(field);
}

/* ---------------------------------------------------------------------- */

int FixScafacos::setmask()
{
  int mask = 0;
  //mask |= INITIAL_INTEGRATE;
  //mask |= FINAL_INTEGRATE;
  mask |= PRE_REVERSE;
  mask |= POST_FORCE;
  mask |= MIN_POST_FORCE;
  mask |= THERMO_ENERGY;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixScafacos::init()
{
  // error checks

  if (domain->dimension == 2)
    error->all(FLERR,"Fix scafacos requires 3d problem");

  rank = comm->me;

  int nlocal = 0;
  int proc_grid[3] = {comm->procgrid[0], comm->procgrid[1], comm->procgrid[2]};
  int myloc[3] = {comm->myloc[0], comm->myloc[1], comm->myloc[2]};

  // call ScaFaCoS init
  // TODO: check if universe->uworld is a good idea for computations
  result = fcs_init(&fcs,method.c_str(),universe->uworld);
  if (!check_result(result, rank)) return;

  setup_handle();

  nlocal = atom -> nlocal;
  // get particle data
  x = &atom->x[0][0];
  q = atom->q; 

  if (tolerance_set)
  {
    result = fcs_set_tolerance(fcs,tolerance_type,tolerance_value);
    if (!check_result(result, rank)) return;
  }

  // print out parameters within handle
  if (rank == 0) fcs_print_parameters(fcs);

  // call the tuning routine (needs to be redone, if critical system 
  // parameters should change)
  result = fcs_tune(fcs,nlocal,x,q);
  if (!check_result(result, rank)) return;

  // allocate arrays larger in order to avoid fast 
  // reallocation due to fluctuations
  local_array_size = (int)(1.25 * (double)nlocal);

  // allocate new result arrays
  memory->create(pot,local_array_size,"scafacos:potential");
  memory->create(field,local_array_size*3,"scafacos:field");
}

/* ---------------------------------------------------------------------- */

void FixScafacos::init_list(int id, NeighList *ptr)
{
}

/* ---------------------------------------------------------------------- */

void FixScafacos::setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixScafacos::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixScafacos::setup_pre_reverse(int eflag, int vflag)
{
  pre_reverse(eflag,vflag);
}

/* ----------------------------------------------------------------------
   integrate electronic degrees of freedom
------------------------------------------------------------------------- */

void FixScafacos::initial_integrate(int vflag) {}

/* ----------------------------------------------------------------------
   store eflag, so can use it in post_force to tally per-atom energies
------------------------------------------------------------------------- */

void FixScafacos::pre_reverse(int eflag, int vflag)
{
  int eflag_caller = eflag;
}

/* ---------------------------------------------------------------------- */

void FixScafacos::post_force(int vflag)
{

  int nlocal;

  nlocal = atom->nlocal;
  x = &atom->x[0][0];
  q = atom->q;

  // check if box has changed since the last call of fcs_tune, 
  // if it has, call fcs_tune
  if (box_has_changed())
  {
    setup_handle();

    // print out parameters within handle TODO: should be done in C style
    /*
    if (rank == 0)
    {
      std::cout << " updated ScaFaCoS handle: " << std::endl;
      fcs_print_parameters(fcs);
    }
    */

    // call the tuning routine (needs to be redone, if critical 
    // system parameters should change)
    result = fcs_tune(fcs,nlocal,x,q);
    if (!check_result(result, rank)) return;
  }
  
  // check if arrays for potentials and field are still large enough 
  // for the number of particles on process
  if (nlocal > local_array_size)
  {
    // allocate arrays larger in order to avoid fast 
    // reallocation due to fluctuations
    local_array_size = (int)(1.25 * (double)nlocal);

    // destroy old result arrays
    memory->destroy(pot);
    memory->destroy(field);

    // allocate result arrays
    memory->create(pot,local_array_size,"scafacos:potential");
    memory->create(field,3*local_array_size,"scafacos:field");
     
  }

  // set result vectors to zero 
  for ( int i = 0; i < nlocal; ++i)
  {
    pot[i] = 0.0;
    field[3*i] = field[3*i+1] = field[3*i+2] = 0.0;
  }
 
  // compute Coulomb
  fcs_run(fcs, nlocal, x, q, field, pot);
  if(!check_result(result,rank)) return;

  double E_coul_loc = 0.0;

  for (int i = 0; i < atom->nlocal; ++i)
  {
    //std::cout << atom->f[i][0] << " " << field[3*i] << " " << q[i] << std::endl;
    atom->f[i][0] += field[3*i] * q[i];
    atom->f[i][1] += field[3*i+1] * q[i];
    atom->f[i][2] += field[3*i+2] * q[i];
    E_coul_loc += 0.5 * q[i] * pot[i];
  } 

  force->pair->eng_coul += E_coul_loc;

}

/* ---------------------------------------------------------------------- */

void FixScafacos::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   integrate electronic degrees of freedom
------------------------------------------------------------------------- */

void FixScafacos::final_integrate() {}

/* ---------------------------------------------------------------------- */

void FixScafacos::reset_dt()
{
  //dtv = update->dt;
  //dtf = 0.5 * update->dt * force->ftm2v;
}

/* ----------------------------------------------------------------------
   DFTB energy from LATTE
------------------------------------------------------------------------- */

double FixScafacos::compute_scalar()
{
}

/* ----------------------------------------------------------------------
   memory usage of local arrays
------------------------------------------------------------------------- */

double FixScafacos::memory_usage()
{
  double bytes = 0.0;
  bytes += local_array_size * sizeof(double);
  bytes += local_array_size * 3 * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
    setup of ScaFaCoS handle with common parameters 
------------------------------------------------------------------------- */
void FixScafacos::setup_handle()
{
  // store periodicity
  periodicity[0] = domain->xperiodic;
  periodicity[1] = domain->yperiodic;
  periodicity[2] = domain->zperiodic;

  // store offset of the system
  offset[0] = domain->boundary[0][0];
  offset[1] = domain->boundary[1][0];
  offset[2] = domain->boundary[2][0];

  // calculate box vectors
  box_x[0] = domain->prd[0];
  box_x[1] = box_x[2] = 0.0;

  box_y[1] = domain->prd[1];
  box_y[0] = box_y[2] = 0.0;

  box_z[2] = domain->prd[2];
  box_z[1] = box_z[0] = 0.0;

  total_particles = atom->natoms;

  // TODO: for now disable short-range calculations within LAMMPS
  near_field_flag = 0;

  // enter all parameters required to ScaFaCoS handle
  result = fcs_set_box_a(fcs, box_x);
  if (!check_result(result, rank)) return;

  result = fcs_set_box_b(fcs, box_y);
  if (!check_result(result, rank)) return;

  result = fcs_set_box_c(fcs, box_z);
  if (!check_result(result, rank)) return;

  result = fcs_set_box_origin(fcs, offset);
  if (!check_result(result, rank)) return;

  result = fcs_set_periodicity(fcs, periodicity);
  if (!check_result(result, rank)) return; 

  result = fcs_set_near_field_flag(fcs, near_field_flag);
  if (!check_result(result, rank)) return;

  result = fcs_set_total_particles(fcs, atom->natoms);
  if (!check_result(result, rank)) return;
}

/* ----------------------------------------------------------------------
    check if box parameters changed, requiring a new call to fcs_tune
------------------------------------------------------------------------- */
bool FixScafacos::box_has_changed()
{
  bool changed = false;

  double n_periodicity[3];
  double n_offset[3];
  double n_box_x[3];
  double n_box_y[3];
  double n_box_z[3];

  int n_total_particles;

  // store periodicity
  n_periodicity[0] = domain->xperiodic;
  n_periodicity[1] = domain->yperiodic;
  n_periodicity[2] = domain->zperiodic;

  // store offset of the system
  n_offset[0] = domain->boundary[0][0];
  n_offset[1] = domain->boundary[1][0];
  n_offset[2] = domain->boundary[2][0];

  // calculate box vectors
  n_box_x[0] = domain->prd[0];
  n_box_x[1] = n_box_x[2] = 0.0;

  n_box_y[1] = domain->prd[1];
  n_box_y[0] = n_box_y[2] = 0.0;

  n_box_z[2] = domain->prd[2];
  n_box_z[1] = n_box_z[0] = 0.0;

  n_total_particles = atom->natoms;

  changed = changed ||
              ( n_periodicity[0] != periodicity[0] ) ||
              ( n_periodicity[1] != periodicity[1] ) ||
              ( n_periodicity[2] != periodicity[2] ) ||
              ( n_offset[0] != offset[0] ) ||
              ( n_offset[1] != offset[1] ) ||
              ( n_offset[2] != offset[2] ) ||
              ( n_box_x[0] != box_x[0] ) ||
              ( n_box_x[1] != box_x[1] ) ||
              ( n_box_x[2] != box_x[2] ) ||
              ( n_box_y[0] != box_y[0] ) ||
              ( n_box_y[1] != box_y[1] ) ||
              ( n_box_y[2] != box_y[2] ) ||
              ( n_box_z[0] != box_z[0] ) ||
              ( n_box_z[1] != box_z[1] ) ||
              ( n_box_z[2] != box_z[2] ) ||
              ( n_total_particles != total_particles );
  return changed;
}



/* ----------------------------------------------------------------------
    check of ScaFaCoS result
------------------------------------------------------------------------- */
bool FixScafacos::check_result(FCSResult result, int comm_rank) 
{
  if (result) 
  {
    printf("ScaFaCoS Error: Caught error on task %d.\n", comm_rank);
    std::string err_msg;
    std::stringstream ss;

    ss << fcs_result_get_function(result) << "\n" 
       << fcs_result_get_message(result) << "\n";
    err_msg = ss.str();

    error -> all(FLERR, err_msg.c_str());
    fcs_result_destroy(result);
  }
  return true;
}
 
