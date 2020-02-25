/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   Fix_num_diff was created by Charles Sievers (UC Davis)
------------------------------------------------------------------------- */

#include "fix_num_diff.h"
#include <algorithm>
#include <mpi.h>
#include <cstring>
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "comm.h"
#include "compute.h"
#include "domain.h"
#include "region.h"
#include "respa.h"
#include "pair.h"
#include "bond.h"
#include "angle.h"
#include "dihedral.h"
#include "improper.h"
#include "kspace.h"
#include "input.h"
#include "neighbor.h"
#include "timer.h"
#include "variable.h"
#include "memory.h"
#include "error.h"
#include "force.h"
#include "group.h"

using namespace LAMMPS_NS;
using namespace FixConst;

//enum{};

/* ---------------------------------------------------------------------- */

FixNumDiff::FixNumDiff(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), numdiff_forces(NULL), temp_f(NULL), temp_x(NULL), id_pe(NULL)
{
  respa_level_support = 1;
  ilevel_respa = 0;
  eflag = 1;
  energy = 0.0;
  size_peratom_cols = 3;

  if (narg < 5) error->all(FLERR,"Illegal fix numdiff command");

  nevery = force->inumeric(FLERR,arg[3]);
  del = force->numeric(FLERR,arg[4]);

  if (nevery <= 0)
    error->all(FLERR,"Illegal fix numdiff command");

  peratom_flag = 1;
  if (force->pair && force->pair->compute_flag) pair_compute_flag = 1;
  else pair_compute_flag = 0;
  if (force->kspace && force->kspace->compute_flag) kspace_compute_flag = 1;
  else kspace_compute_flag = 0;

  numdiff_forces = NULL;

  maxatom1 = 0;

  int n = strlen(id) + 6;
  id_pe = new char[n];
  strcpy(id_pe,id);
  strcat(id_pe,"_pe");

  char **newarg = new char*[10];
  newarg[0] = id_pe;
  newarg[1] = (char *) "all";
  newarg[2] = (char *) "pe";
  newarg[3] = (char *) "pair";
  newarg[4] = (char *) "bond";
  newarg[5] = (char *) "angle";
  newarg[6] = (char *) "dihedral";
  newarg[7] = (char *) "improper";
  newarg[8] = (char *) "kspace";
  newarg[9] = (char *) "fix";
  modify->add_compute(3,newarg);
  delete [] newarg;

  nmax = 0;
  memory->create(numdiff_forces,atom->natoms,3,"numdiff:numdiff_force");
  force_clear(numdiff_forces);
  array_atom = numdiff_forces;
}

/* ---------------------------------------------------------------------- */

FixNumDiff::~FixNumDiff()
{
  memory->destroy(numdiff_forces);
  memory->destroy(temp_f);
  memory->destroy(temp_x);

  modify->delete_compute(id_pe);
  delete [] id_pe;
}

/* ---------------------------------------------------------------------- */

int FixNumDiff::setmask()
{
  datamask_read = datamask_modify = 0;

  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixNumDiff::init()
{
  // check variables

  int icompute = modify->find_compute(id_pe);
  if (icompute < 0)
    error->all(FLERR,"Compute ID for fix numdiff does not exist");
  pe = modify->compute[icompute];

  if (strstr(update->integrate_style,"respa")) {
    ilevel_respa = ((Respa *) update->integrate)->nlevels-1;
    if (respa_level >= 0) ilevel_respa = MIN(respa_level,ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixNumDiff::post_force(int vflag)
{
  if (update->ntimestep % nevery) return;

  // energy and virial setup

  calculate_forces(vflag);

  // if (lmp->kokkos)
  //   atom->sync_modify(Host, (unsigned int) (F_MASK | MASK_MASK),
  //                     (unsigned int) F_MASK);
}

/* ---------------------------------------------------------------------- */

void FixNumDiff::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == ilevel_respa) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixNumDiff::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   create dynamical matrix
------------------------------------------------------------------------- */

void FixNumDiff::calculate_forces(int vflag)
{
  int local_idx; // local index
  bigint natoms = atom->natoms;
  double **f = atom->f;
  double **x = atom->x;
  double temp = 1.0 / 2 / del;
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  int flag = 0;
  int allflag = 0;

  if (atom->nmax > maxatom1) {
    memory->destroy(numdiff_forces);
    memory->destroy(temp_f);
    memory->destroy(temp_x);
    maxatom1 = atom->nmax;
    memory->create(numdiff_forces,maxatom1,3,"numdiff:numdiff_force");
    memory->create(temp_f,maxatom1,3,"numdiff:temp_f");
    memory->create(temp_x,maxatom1,3,"numdiff:temp_x");
    array_atom = numdiff_forces;
  }

  for (bigint i = 0; i < natoms; i++)
    for (int j = 0; j < 3; j++)
      temp_f[i][j] = f[i][j];

  for (bigint i = 0; i < natoms; i++)
    for (int j = 0; j < 3; j++)
      temp_x[i][j] = x[i][j];

  //initialize forces to all zeros
  force_clear(numdiff_forces);

  for (bigint i=1; i<=natoms; i++){
    local_idx = atom->map(i);
    if (mask[local_idx] & groupbit) flag = 1;
    else flag = 0;
    MPI_Allreduce(&flag, &allflag,1, MPI_INT, MPI_SUM,world);
    if (!allflag) continue;
    for (int alpha=0; alpha<3; alpha++){
      displace_atom(local_idx, alpha, 1);
      update_energy(vflag);
      if (local_idx >= 0 && local_idx < nlocal)
        numdiff_forces[local_idx][alpha] -= energy;

      displace_atom(local_idx,alpha,-2);
      update_energy(vflag);
      if (local_idx >= 0 && local_idx < nlocal) {
        numdiff_forces[local_idx][alpha] += energy;
        numdiff_forces[local_idx][alpha] *= temp;
      }
      displace_atom(local_idx,alpha,1);
    }
  }

  for (bigint i = 0; i < natoms; i++)
    for (int j = 0; j < 3; j++)
      f[i][j] = temp_f[i][j];

  for (bigint i = 0; i < natoms; i++)
    for (int j = 0; j < 3; j++)
      x[i][j] = temp_x[i][j];

}

/* ----------------------------------------------------------------------
  Displace atoms
   ---------------------------------------------------------------------- */

void FixNumDiff::displace_atom(int local_idx, int direction, int magnitude)
{
  if (local_idx < 0) return;

  double **x = atom->x;
  int *sametag = atom->sametag;
  int j = local_idx;
  x[local_idx][direction] += del*magnitude;

  while (sametag[j] >= 0){
    j = sametag[j];
    x[j][direction] += del*magnitude;
  }
}

/* ----------------------------------------------------------------------
   evaluate potential energy and forces
   may migrate atoms due to reneighboring
   return new energy, which should include nextra_global dof
   return negative gradient stored in atom->f
   return negative gradient for nextra_global dof in fextra
------------------------------------------------------------------------- */

void FixNumDiff::update_energy(int vflag)
{
  force_clear(atom->f);

  if (pair_compute_flag) {
    force->pair->compute(eflag,vflag);
    timer->stamp(Timer::PAIR);
  }
  if (atom->molecular) {
    if (force->bond) {
      force->bond->compute(eflag,vflag);
    } if (force->angle) {
      force->angle->compute(eflag,vflag);
    } if (force->dihedral) {
      force->dihedral->compute(eflag,vflag);
    } if (force->improper) {
      force->improper->compute(eflag,vflag);
    }
    timer->stamp(Timer::BOND);
  }
  if (kspace_compute_flag) {
    force->kspace->compute(eflag,vflag);
    timer->stamp(Timer::KSPACE);
  }

  energy = pe->compute_scalar();
}

/* ----------------------------------------------------------------------
   clear forces needed
------------------------------------------------------------------------- */

void FixNumDiff::force_clear(double **forces)
{

  size_t nbytes = sizeof(double) * atom->nlocal;
  if (force->newton) nbytes += sizeof(double) * atom->nghost;

  if (nbytes) {
    memset(&forces[0][0],0,3*nbytes);
  }
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

int FixNumDiff::pack_exchange(int i, double *buf)
{
  int n = 0;
  buf[n++] = numdiff_forces[i][0];
  buf[n++] = numdiff_forces[i][1];
  buf[n++] = numdiff_forces[i][2];
  return n;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
------------------------------------------------------------------------- */

int FixNumDiff::unpack_exchange(int nlocal, double *buf)
{
  int n = 0;
  numdiff_forces[nlocal][0] = buf[n++];
  numdiff_forces[nlocal][1] = buf[n++];
  numdiff_forces[nlocal][2] = buf[n++];
  return n;
}

/* ----------------------------------------------------------------------
   return force calculated by numerical difference
------------------------------------------------------------------------- */

double FixNumDiff::compute_array(int i, int j)
{
  return numdiff_forces[i][j];
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixNumDiff::memory_usage()
{
  bigint bytes = 0.0;
  bytes += 3 * atom->natoms * 3 * sizeof(double); // temp_f, temp_x, numdiff_f
  return bytes;
}
