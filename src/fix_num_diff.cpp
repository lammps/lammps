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
#include <cstdlib>
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

//TODO: Add energy from potentials to single energy array and then reverse comm pack and unpack

using namespace LAMMPS_NS;
using namespace FixConst;

//enum{};

/* ---------------------------------------------------------------------- */

FixNumDiff::FixNumDiff(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), local_forces(NULL), global_forces(NULL), groupmap(NULL), fp(NULL)
{
  respa_level_support = 1;
  ilevel_respa = 0;
  eflag = 1;

  if (narg < 5) error->all(FLERR,"Illegal fix numdiff command");

  nevery = force->inumeric(FLERR,arg[3]);
  del = force->numeric(FLERR,arg[4]);

  if (nevery <= 0)
    error->all(FLERR,"Illegal fix ave/atom command");

  peratom_flag = 1;

  igroup = group->find(arg[1]);
  if (igroup == -1) error->all(FLERR,"Could not find numerical difference group ID");
  groupbit = group->bitmask[igroup];
  gcount = group->count(igroup);
  flen = gcount*3;
  memory->create(groupmap,atom->natoms,"total_group_map:total_group_map");
  local_forces = NULL;
  global_forces = NULL;

  memory->grow(local_forces,gcount,3,"fix_numdiff:local_forces");
  memory->grow(global_forces,gcount,3,"fix_numdiff:global_forces");
  memory->grow(temp_f, atom->natoms, 3, "fix_numdiff:temporary_forces");
  memory->grow(e, atom->natoms, "fix_numdiff:per_atom_energy");
}

/* ---------------------------------------------------------------------- */

FixNumDiff::~FixNumDiff()
{
  if (fp && me == 0) fclose(fp);
  memory->destroy(groupmap);
  memory->destroy(local_forces);
  memory->destroy(global_forces);
  memory->destroy(temp_f);
  memory->destroy(e);
  fp = NULL;
}

/* ---------------------------------------------------------------------- */

int FixNumDiff::setmask()
{
  datamask_read = datamask_modify = 0;

  int mask = 0;
  mask |= POST_FORCE;
  mask |= THERMO_ENERGY;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixNumDiff::init()
{
  // check variables

  create_groupmap();

  if (strstr(update->integrate_style,"respa")) {
    ilevel_respa = ((Respa *) update->integrate)->nlevels-1;
    if (respa_level >= 0) ilevel_respa = MIN(respa_level,ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixNumDiff::setup(int vflag)
{
  ev_setup(eflag, vflag);

  //if all then skip communication groupmap population
  if (gcount == atom->natoms)
    for (bigint i=0; i<atom->natoms; i++)
      groupmap[i] = i;
  else
    create_groupmap();
  if (strstr(update->integrate_style,"verlet"))
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(ilevel_respa);
    post_force_respa(vflag,ilevel_respa,0);
    ((Respa *) update->integrate)->copy_f_flevel(ilevel_respa);
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

  size_t nbytes = sizeof(double) * natoms * 3;
  memcpy(temp_f, f, nbytes);

  //initialize forces to all zeros
  nd_force_clear(local_forces);
  nd_force_clear(global_forces);

  // if (comm->me == 0 && screen) {
  //   fprintf(screen,"Calculating Forces ...\n");
  //   fprintf(screen,"  Total # of atoms = " BIGINT_FORMAT "\n", natoms);
  //   fprintf(screen,"  Atoms in group = " BIGINT_FORMAT "\n", gcount);
  //   fprintf(screen,"  Total Force elements = " BIGINT_FORMAT "\n", (flen) );
  // }

  update->nsteps = 0;
  int prog = 0;
  for (bigint i=1; i<=natoms; i++){
    local_idx = atom->map(i);
    if (gm[i-1] < 0)
      continue;
    for (int alpha=0; alpha<3; alpha++){
      displace_atom(local_idx, alpha, 1);
      update_force(vflag);
      fprintf(screen, "%f\n", e[local_idx]);
      local_forces[gm[i-1]][alpha] -= e[local_idx];

      displace_atom(local_idx,alpha,-2);
      update_force(vflag);
      local_forces[gm[i-1]][alpha] -= -e[local_idx];
      local_forces[gm[i-1]][alpha] *= temp;
      displace_atom(local_idx,alpha,1);
    }
    //for (int k=0; k<3; k++)
    MPI_Reduce(local_forces,global_forces,flen,MPI_DOUBLE,MPI_SUM,0,world);
    if (comm->me == 0 && screen) {
      int p = 10 * gm[i-1] / gcount;
      if (p > prog) {
        prog = p;
        fprintf(screen," %d%%",p*10);
        fflush(screen);
      }
    }
  }

  memcpy(f, temp_f, nbytes);
  for (int i = 0; i < natoms; i++)
    for (int j = 0; j < 3; j++)
      fprintf(screen, "%f %f\n", f[i][j], temp_f[i][j]);

  if (comm->me == 0 && screen) fprintf(screen,"\n");

  if (screen && me ==0 ) fprintf(screen,"Finished Calculating Dynamical Matrix\n");
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
  int n_post_force = modify->n_post_force;

  if (pair_compute_flag) {
    double *eatom = force->pair->eatom;
    force->pair->compute(eflag,vflag);
    timer->stamp(Timer::PAIR);
  }
  if (atom->molecular) {
    if (force->bond) {
      double *eatom = force->bond->eatom;
      force->bond->compute(eflag,vflag);
    } if (force->angle) {
      double *eatom = force->angle->eatom;
      force->angle->compute(eflag,vflag);
    } if (force->dihedral) {
      double *eatom = force->dihedral->eatom;
      force->dihedral->compute(eflag,vflag);
    } if (force->improper) {
      double *eatom = force->improper->eatom;
      force->improper->compute(eflag,vflag);
    }
    timer->stamp(Timer::BOND);
  }
  if (kspace_compute_flag) {
    double *eatom = force->kspace->eatom;
    force->kspace->compute(eflag,vflag);
    timer->stamp(Timer::KSPACE);
  }
  if (force->newton) {
    comm->reverse_comm();
    timer->stamp(Timer::COMM);
  }

  // force modifications

  //if (n_post_force) modify->post_force(vflag);
  //timer->stamp(Timer::MODIFY);

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
  bytes += atom->natoms * sizeof(double); // e
  return bytes;
}
