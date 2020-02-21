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

#include "fix_num_diff.h"
#include <mpi.h>
#include <cstring>
#include <cmath>
#include "atom.h"
#include "atom_masks.h"
#include "update.h"
#include "modify.h"
#include "comm.h"
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
#include <algorithm>

//TODO: Add energy from potentials to single energy array and then reverse comm pack and unpack

using namespace LAMMPS_NS;
using namespace FixConst;

//enum{};

/* ---------------------------------------------------------------------- */

FixNumDiff::FixNumDiff(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), groupmap(NULL), local_forces(NULL), global_forces(NULL), temp_f(NULL)
{
  respa_level_support = 1;
  ilevel_respa = 0;
  eflag = 1;
  energy = 0;

  if (narg < 5) error->all(FLERR,"Illegal fix numdiff command");

  nevery = force->inumeric(FLERR,arg[3]);
  del = force->numeric(FLERR,arg[4]);

  if (nevery <= 0)
    error->all(FLERR,"Illegal fix ave/atom command");

  peratom_flag = 1;
  if (force->pair && force->pair->compute_flag) pair_compute_flag = 1;
  else pair_compute_flag = 0;
  if (force->kspace && force->kspace->compute_flag) kspace_compute_flag = 1;
  else kspace_compute_flag = 0;

  igroup = group->find(arg[1]);
  if (igroup == -1) error->all(FLERR,"Could not find numerical difference group ID");
  groupbit = group->bitmask[igroup];
  gcount = group->count(igroup);
  flen = gcount*3;
  memory->create(groupmap,atom->natoms,"total_group_map:total_group_map");
  local_forces = NULL;
  global_forces = NULL;
  array_flag = 1;
  size_array_cols = 3;
  size_array_rows = gcount;

  memory->create(local_forces,gcount,3,"fix_numdiff:local_forces");
  memory->create(global_forces,gcount,3,"fix_numdiff:global_forces");
  memory->grow(temp_f, atom->natoms, 3, "fix_numdiff:temporary_forces");

  for (bigint i = 0; i < gcount; i++)
    for (int j = 0; j < 3; j++) {
      global_forces[i][j] = 0;
      local_forces[i][j] = 0;
    }

  nmax = 0;
}

/* ---------------------------------------------------------------------- */

FixNumDiff::~FixNumDiff()
{
  memory->destroy(groupmap);
  memory->destroy(local_forces);
  memory->destroy(global_forces);
  memory->destroy(temp_f);
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

  if (gcount == atom->natoms)
    for (bigint i=0; i<atom->natoms; i++)
      groupmap[i] = i;
  else
    create_groupmap();

  if (strstr(update->integrate_style,"respa")) {
    ilevel_respa = ((Respa *) update->integrate)->nlevels-1;
    if (respa_level >= 0) ilevel_respa = MIN(respa_level,ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixNumDiff::min_setup(int vflag)
{
  post_force(vflag);
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
  bigint *gm = groupmap;
  double **f = atom->f;
  double temp = 1.0 / 2 / del;
  int nlocal = atom->nlocal;

  for (bigint i = 0; i < natoms; i++)
    for (int j = 0; j < 3; j++)
      temp_f[i][j] = f[i][j];

  //initialize forces to all zeros
  nd_force_clear(local_forces);
  nd_force_clear(global_forces);

  for (bigint i=1; i<=natoms; i++){
    local_idx = atom->map(i);
    if (gm[i-1] < 0)
      continue;
    for (int alpha=0; alpha<3; alpha++){
      displace_atom(local_idx, alpha, 1);
      update_force(vflag);
      if (local_idx >= 0 && local_idx < nlocal)
        local_forces[gm[i-1]][alpha] -= energy;

      displace_atom(local_idx,alpha,-2);
      update_force(vflag);
      if (local_idx >= 0 && local_idx < nlocal) {
        local_forces[gm[i-1]][alpha] += energy;
        local_forces[gm[i-1]][alpha] *= temp;
      }
      displace_atom(local_idx,alpha,1);
    }
    for (int k=0; k<gcount; k++)
      MPI_Allreduce(local_forces[k],global_forces[k],3,MPI_DOUBLE,MPI_SUM,world);
  }

  for (bigint i = 0; i < natoms; i++)
    for (int j = 0; j < 3; j++)
      f[i][j] = temp_f[i][j];

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

void FixNumDiff::update_force(int vflag)
{
  force_clear();

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
  // if (force->newton) {
  //   comm->reverse_comm();
  //   timer->stamp(Timer::COMM);
  // }
  compute_energy();


}

/* ----------------------------------------------------------------------
   clear force on own & ghost atoms
   clear other arrays as needed
------------------------------------------------------------------------- */

void FixNumDiff::force_clear()
{
  if (external_force_clear) return;

  // clear global force array
  // if either newton flag is set, also include ghosts

  size_t nbytes = sizeof(double) * atom->nlocal;
  if (force->newton) nbytes += sizeof(double) * atom->nghost;

  if (nbytes) {
    memset(&atom->f[0][0],0,3*nbytes);
  }
}

/* ----------------------------------------------------------------------
   clear forces needed
------------------------------------------------------------------------- */

void FixNumDiff::nd_force_clear(double **forces)
{

  size_t nbytes = sizeof(double) * flen;


  if (nbytes) {
      memset(&forces[0][0],0,nbytes);
  }
}

/* ---------------------------------------------------------------------- */
void FixNumDiff::compute_energy()
{
  double one = 0.0;
  if (pair_compute_flag)
    one += force->pair->eng_vdwl + force->pair->eng_coul;

  if (atom->molecular) {
    if (force->bond) one += force->bond->energy;
    if (force->angle) one += force->angle->energy;
    if (force->dihedral) one += force->dihedral->energy;
    if (force->improper) one += force->improper->energy;
  }

  MPI_Allreduce(&one,&energy,1,MPI_DOUBLE,MPI_SUM,world);

  if (kspace_compute_flag) energy += force->kspace->energy;

  if (pairflag && force->pair && force->pair->tail_flag) {
    double volume = domain->xprd * domain->yprd * domain->zprd;
    energy += force->pair->etail / volume;
  }
}

/* ---------------------------------------------------------------------- */

void FixNumDiff::create_groupmap()
{
  //Create a group map which maps atom order onto group
  // groupmap[global atom index-1] = output column/row

  int local_idx; // local index
  int gid = 0; //group index
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  bigint natoms = atom->natoms;
  int *recv = new int[comm->nprocs];
  int *displs = new int[comm->nprocs];
  bigint *temp_groupmap = new bigint[natoms];

  //find number of local atoms in the group (final_gid)
  for (bigint i=1; i<=natoms; i++){
    local_idx = atom->map(i);
    if ((local_idx >= 0) && (local_idx < nlocal) && mask[local_idx] & groupbit)
      gid += 1; // gid at the end of loop is final_Gid
  }
  //create an array of length final_gid
  bigint *sub_groupmap = new bigint[gid];

  gid = 0;
  //create a map between global atom id and group atom id for each proc
  for (bigint i=1; i<=natoms; i++){
    local_idx = atom->map(i);
    if ((local_idx >= 0) && (local_idx < nlocal) && mask[local_idx] & groupbit){
      sub_groupmap[gid] = i;
      gid += 1;
    }
  }

  //populate arrays for Allgatherv
  for (int i=0; i<comm->nprocs; i++){
    recv[i] = 0;
  }
  recv[comm->me] = gid;
  MPI_Allreduce(recv,displs,comm->nprocs,MPI_INT,MPI_SUM,world);
  for (int i=0; i<comm->nprocs; i++){
    recv[i]=displs[i];
    if (i>0) displs[i] = displs[i-1]+recv[i-1];
    else displs[i] = 0;
  }

  //combine subgroup maps into total temporary groupmap
  MPI_Allgatherv(sub_groupmap,gid,MPI_LMP_BIGINT,temp_groupmap,recv,displs,MPI_LMP_BIGINT,world);
  std::sort(temp_groupmap,temp_groupmap+gcount);

  //populate member groupmap based on temp groupmap
  bigint j = 0;
  for (bigint i=1; i<=natoms; i++){
    // flag groupmap contents that are in temp_groupmap
    if (j < gcount && i == temp_groupmap[j])
      groupmap[i-1] = j++;
    else
      groupmap[i-1] = -1;
  }

  //free that memory!
  delete[] recv;
  delete[] displs;
  delete[] sub_groupmap;
  delete[] temp_groupmap;
}

/* ----------------------------------------------------------------------
   return force calculated by numerical difference
------------------------------------------------------------------------- */

double FixNumDiff::compute_array(int i, int j)
{
  return global_forces[i][j];
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixNumDiff::memory_usage()
{
  bigint bytes = 0.0;
  bytes += atom->natoms * 3 * sizeof(double); // temp_f
  bytes += 2 * flen * sizeof(double); // local_f, global_f
  bytes += gcount * sizeof(bigint); // groupmap
  return bytes;
}
