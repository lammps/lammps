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
   Contributing author: Hendrik Heenen (Technical University of Munich)
                        (hendrik.heenen at mytum.com)
------------------------------------------------------------------------- */

#include "compute_temp_cs.h"
#include <mpi.h>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "atom.h"
#include "atom_vec.h"
#include "domain.h"
#include "update.h"
#include "force.h"
#include "group.h"
#include "modify.h"
#include "fix.h"
#include "fix_store.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeTempCS::ComputeTempCS(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg), vint(NULL), id_fix(NULL), fix(NULL)
{
  if (narg != 5) error->all(FLERR,"Illegal compute temp/cs command");

  if (atom->avec->bonds_allow == 0)
    error->all(FLERR,"Compute temp/cs used when bonds are not allowed");

  scalar_flag = vector_flag = 1;
  size_vector = 6;
  extscalar = 0;
  extvector = 1;
  tempflag = 1;
  tempbias = 1;
  extarray = 0;

  // find and define groupbits for core and shell groups

  cgroup = group->find(arg[3]);
  if (cgroup == -1)
    error->all(FLERR,"Cannot find specified group ID for core particles");
  groupbit_c = group->bitmask[cgroup];

  sgroup = group->find(arg[4]);
  if (sgroup == -1)
    error->all(FLERR,"Cannot find specified group ID for shell particles");
  groupbit_s = group->bitmask[sgroup];

  // create a new fix STORE style
  // id = compute-ID + COMPUTE_STORE, fix group = compute group

  int n = strlen(id) + strlen("_COMPUTE_STORE") + 1;
  id_fix = new char[n];
  strcpy(id_fix,id);
  strcat(id_fix,"_COMPUTE_STORE");

  char **newarg = new char*[6];
  newarg[0] = id_fix;
  newarg[1] = group->names[igroup];
  newarg[2] = (char *) "STORE";
  newarg[3] = (char *) "peratom";
  newarg[4] = (char *) "0";
  newarg[5] = (char *) "1";
  modify->add_fix(6,newarg);
  fix = (FixStore *) modify->fix[modify->nfix-1];
  delete [] newarg;

  // set fix store values = 0 for now
  // fill them in via setup() once Comm::borders() has been called
  // skip if resetting from restart file

  if (fix->restart_reset) {
    fix->restart_reset = 0;
    firstflag = 0;
  } else {
    double *partner = fix->vstore;
    int nlocal = atom->nlocal;
    for (int i = 0; i < nlocal; i++) partner[i] = ubuf(0).d;
    firstflag = 1;
  }

  // allocate memory

  vector = new double[6];
  maxatom = 0;
  vint = NULL;

  // set comm size needed by this Compute

  comm_reverse = 1;
}

/* ---------------------------------------------------------------------- */

ComputeTempCS::~ComputeTempCS()
{
  // check nfix in case all fixes have already been deleted

  if (modify->nfix) modify->delete_fix(id_fix);

  delete [] id_fix;
  delete [] vector;
  memory->destroy(vint);
}

/* ---------------------------------------------------------------------- */

void ComputeTempCS::init()
{
  if (comm->ghost_velocity == 0)
    error->all(FLERR,"Compute temp/cs requires ghost atoms store velocity");
}

/* ---------------------------------------------------------------------- */

void ComputeTempCS::setup()
{
  if (firstflag) {
    firstflag = 0;

    // insure # of core atoms = # of shell atoms

    int ncores = group->count(cgroup);
    nshells = group->count(sgroup);
    if (ncores != nshells)
      error->all(FLERR,"Number of core atoms != number of shell atoms");

    // for each C/S pair:
    // set partner IDs of both atoms if this atom stores bond between them
    // will set partner IDs for ghost atoms if needed by another proc
    // nall loop insures all ghost atom partner IDs are set before reverse comm

    int *num_bond = atom->num_bond;
    tagint **bond_atom = atom->bond_atom;
    tagint *tag = atom->tag;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    double *partner = fix->vstore;
    tagint partnerID;

    int nall = nlocal + atom->nghost;
    for (int i = nlocal; i < nall; i++) partner[i] = ubuf(0).d;

    int i,j,m,match;
    for (i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit_c || mask[i] & groupbit_s) {
        for (m = 0; m < num_bond[i]; m++) {
          partnerID = bond_atom[i][m];
          j = atom->map(partnerID);
          if (j == -1) error->one(FLERR,"Core/shell partner atom not found");
          match = 0;
          if (mask[i] & groupbit_c && mask[j] & groupbit_s) match = 1;
          if (mask[i] & groupbit_s && mask[j] & groupbit_c) match = 1;
          if (match) {
            partner[i] = ubuf(partnerID).d;
            partner[j] = ubuf(tag[i]).d;
          }
        }
      }
    }

    // reverse comm to acquire unknown partner IDs from ghost atoms
    // only needed if newton_bond = on

    if (force->newton_bond) comm->reverse_comm_compute(this);

    // check that all C/S partners were found

    int flag = 0;
    for (i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit_c || mask[i] & groupbit_s) {
        partnerID = (tagint) ubuf(partner[i]).i;
        if (partnerID == 0) flag = 1;
      }
    }

    int flagall;
    MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
    if (flagall) error->all(FLERR,"Core/shell partners were not all found");
  }

  // calculate DOF for temperature

  dof_compute();
}

/* ---------------------------------------------------------------------- */

void ComputeTempCS::dof_compute()
{
  adjust_dof_fix();
  int nper = domain->dimension;
  natoms_temp = group->count(igroup);
  dof = nper * natoms_temp;
  dof -= nper * nshells;
  dof -= extra_dof + fix_dof;
  if (dof > 0) tfactor = force->mvv2e / (dof * force->boltz);
  else tfactor = 0.0;
}

/* ---------------------------------------------------------------------- */

double ComputeTempCS::compute_scalar()
{
  double vthermal[3];

  invoked_scalar = update->ntimestep;

  vcm_pairs();

  // calculate thermal scalar in respect to atom velocities as center-of-mass
  // velocities of its according core/shell pairs

  double **v = atom->v;
  int *mask = atom->mask;
  int *type = atom->type;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;

  double t = 0.0;

  for (int i = 0; i < nlocal; i++){
    if (mask[i] & groupbit) {
      vthermal[0] = v[i][0] - vint[i][0];
      vthermal[1] = v[i][1] - vint[i][1];
      vthermal[2] = v[i][2] - vint[i][2];
      if (rmass)
        t += (vthermal[0]*vthermal[0] + vthermal[1]*vthermal[1] +
              vthermal[2]*vthermal[2]) * rmass[i];
      else
        t += (vthermal[0]*vthermal[0] + vthermal[1]*vthermal[1] +
              vthermal[2]*vthermal[2]) * mass[type[i]];
    }
  }

  MPI_Allreduce(&t,&scalar,1,MPI_DOUBLE,MPI_SUM,world);
  if (dynamic) dof_compute();
  if (dof < 0.0 && natoms_temp > 0.0)
    error->all(FLERR,"Temperature compute degrees of freedom < 0");
  scalar *= tfactor;
  return scalar;
}

/* ---------------------------------------------------------------------- */

void ComputeTempCS::compute_vector()
{
  invoked_vector = update->ntimestep;

  double **v = atom->v;
  int *mask = atom->mask;
  int *type = atom->type;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;
  double massone;

  double t[6];
  for (int i = 0; i < 6; i++) t[i] = 0.0;

  for (int i = 0; i < nlocal; i++){
    if (mask[i] & groupbit) {
      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];
      t[0] += massone * v[i][0]*v[i][0];
      t[1] += massone * v[i][1]*v[i][1];
      t[2] += massone * v[i][2]*v[i][2];
      t[3] += massone * v[i][0]*v[i][1];
      t[4] += massone * v[i][0]*v[i][2];
      t[5] += massone * v[i][1]*v[i][2];
    }
  }

  MPI_Allreduce(t,vector,6,MPI_DOUBLE,MPI_SUM,world);
  for (int i = 0; i < 6; i++) vector[i] *= force->mvv2e;
}

/* ---------------------------------------------------------------------- */

void ComputeTempCS::vcm_pairs()
{
  int i,j;
  double massone,masstwo;
  double vcm[3];

  // reallocate vint if necessary

  int nlocal = atom->nlocal;

  if (atom->nmax > maxatom) {
    memory->destroy(vint);
    maxatom = atom->nmax;
    memory->create(vint,maxatom,3,"temp/cs:vint");
  }

  // vcm = COM velocity of each CS pair
  // vint = internal velocity of each C/S atom, used as bias

  double **v = atom->v;
  int *mask = atom->mask;
  int *type = atom->type;
  double *mass = atom->mass;
  double *rmass = atom->rmass;

  double *partner = fix->vstore;
  tagint partnerID;

  for (i = 0; i < nlocal; i++) {
    if ((mask[i] & groupbit) &&
        (mask[i] & groupbit_c || mask[i] & groupbit_s)) {
      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];
      vcm[0] = v[i][0]*massone;
      vcm[1] = v[i][1]*massone;
      vcm[2] = v[i][2]*massone;

      partnerID = (tagint) ubuf(partner[i]).i;
      j = atom->map(partnerID);
      if (j == -1) error->one(FLERR,"Core/shell partner atom not found");

      if (rmass) masstwo = rmass[j];
      else masstwo = mass[type[j]];
      vcm[0] += v[j][0]*masstwo;
      vcm[1] += v[j][1]*masstwo;
      vcm[2] += v[j][2]*masstwo;
      vcm[0] /= (massone + masstwo);
      vcm[1] /= (massone + masstwo);
      vcm[2] /= (massone + masstwo);

      vint[i][0] = v[i][0] - vcm[0];
      vint[i][1] = v[i][1] - vcm[1];
      vint[i][2] = v[i][2] - vcm[2];
    } else vint[i][0] = vint[i][1] = vint[i][2] = 0.0;
  }
}

/* ----------------------------------------------------------------------
   remove velocity bias from atom I to leave thermal velocity
   thermal velocity in this case is COM velocity of C/S pair
------------------------------------------------------------------------- */

void ComputeTempCS::remove_bias(int i, double *v)
{
  v[0] -= vint[i][0];
  v[1] -= vint[i][1];
  v[2] -= vint[i][2];
}

/* ----------------------------------------------------------------------
   remove velocity bias from all atoms to leave thermal velocity
   thermal velocity in this case is COM velocity of C/S pair
------------------------------------------------------------------------- */

void ComputeTempCS::remove_bias_all()
{
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      v[i][0] -= vint[i][0];
      v[i][1] -= vint[i][1];
      v[i][2] -= vint[i][2];
    }
}

/* ----------------------------------------------------------------------
   reset thermal velocity of all atoms to be consistent with bias
   called from velocity command after it creates thermal velocities
   this resets each atom's velocity to COM velocity of C/S pair
------------------------------------------------------------------------- */

void ComputeTempCS::reapply_bias_all()
{
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  // recalculate current COM velocities

  vcm_pairs();

  // zero vint after using ti so that Velocity call to restore_bias_all()
  // will not further alter the velocities within a C/S pair

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      v[i][0] -= vint[i][0];
      v[i][1] -= vint[i][1];
      v[i][2] -= vint[i][2];
      vint[i][0] = 0.0;
      vint[i][1] = 0.0;
      vint[i][2] = 0.0;
    }
}

/* ----------------------------------------------------------------------
   add back in velocity bias to atom I removed by remove_bias()
   assume remove_bias() was previously called
------------------------------------------------------------------------- */

void ComputeTempCS::restore_bias(int i, double *v)
{
  v[0] += vint[i][0];
  v[1] += vint[i][1];
  v[2] += vint[i][2];
}

/* ----------------------------------------------------------------------
   add back in velocity bias to all atoms removed by remove_bias_all()
   assume remove_bias_all() was previously called
------------------------------------------------------------------------- */

void ComputeTempCS::restore_bias_all()
{
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      v[i][0] += vint[i][0];
      v[i][1] += vint[i][1];
      v[i][2] += vint[i][2];
    }
}

/* ---------------------------------------------------------------------- */

int ComputeTempCS::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  double *partner = fix->vstore;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) buf[m++] = partner[i];
  return m;
}

/* ---------------------------------------------------------------------- */

void ComputeTempCS::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;

  double *partner = fix->vstore;
  tagint partnerID;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    partnerID = (tagint) ubuf(buf[m++]).i;
    if (partnerID) partner[j] = ubuf(partnerID).d;
  }
}

/* ----------------------------------------------------------------------
   memory usage of local data
------------------------------------------------------------------------- */

double ComputeTempCS::memory_usage()
{
  double bytes = (bigint) maxatom * 3 * sizeof(double);
  return bytes;
}
