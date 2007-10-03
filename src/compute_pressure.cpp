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

#include "mpi.h"
#include "string.h"
#include "stdlib.h"
#include "compute_pressure.h"
#include "atom.h"
#include "domain.h"
#include "modify.h"
#include "fix.h"
#include "force.h"
#include "pair.h"
#include "bond.h"
#include "angle.h"
#include "dihedral.h"
#include "improper.h"
#include "kspace.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputePressure::ComputePressure(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 4) error->all("Illegal compute pressure command");

  if (igroup) error->all("Compute pressure must use group all");

  // store temperature ID used by pressure computation
  // insure it is valid for temperature computation

  npre = 1;
  id_pre = new char*[1];
  int n = strlen(arg[3]) + 1;
  id_pre[0] = new char[n];
  strcpy(id_pre[0],arg[3]);

  int icompute = modify->find_compute(id_pre[0]);

  if (icompute < 0) error->all("Could not find compute pressure temp ID");
  if (modify->compute[icompute]->tempflag == 0)
    error->all("Compute pressure temp ID does not compute temperature");

  scalar_flag = vector_flag = 1;
  size_vector = 6;
  extensive = 0;
  pressflag = 1;

  vector = new double[6];
  nvirial = 0;
  vptr = NULL;
}

/* ---------------------------------------------------------------------- */

ComputePressure::~ComputePressure()
{
  delete [] vector;
  delete [] vptr;
}

/* ---------------------------------------------------------------------- */

void ComputePressure::init()
{
  boltz = force->boltz;
  nktv2p = force->nktv2p;
  dimension = domain->dimension;

  // set temperature used by pressure

  int icompute = modify->find_compute(id_pre[0]);
  if (icompute < 0) error->all("Could not find compute pressure temp ID");
  temperature = modify->compute[icompute];

  // detect contributions to virial
  // vptr points to all virial[6] contributions

  delete [] vptr;
  nvirial = 0;
  vptr = NULL;

  if (force->pair) nvirial++;
  if (atom->molecular && force->bond) nvirial++;
  if (atom->molecular && force->angle) nvirial++;
  if (atom->molecular && force->dihedral) nvirial++;
  if (atom->molecular && force->improper) nvirial++;
  for (int i = 0; i < modify->nfix; i++)
    if (modify->fix[i]->virial_flag) nvirial++;

  if (nvirial) {
    vptr = new double*[nvirial];
    nvirial = 0;
    if (force->pair) vptr[nvirial++] = force->pair->virial;
    if (force->bond) vptr[nvirial++] = force->bond->virial;
    if (force->angle) vptr[nvirial++] = force->angle->virial;
    if (force->dihedral) vptr[nvirial++] = force->dihedral->virial;
    if (force->improper) vptr[nvirial++] = force->improper->virial;
    for (int i = 0; i < modify->nfix; i++)
      if (modify->fix[i]->virial_flag)
	vptr[nvirial++] = modify->fix[i]->virial;
  }

  // flag Kspace contribution separately, since not summed across procs

  kspaceflag = 0;
  if (force->kspace) {
    kspaceflag = 1;
    kspace_virial = force->kspace->virial;
  }
}

/* ----------------------------------------------------------------------
   compute total pressure, averaged over Pxx, Pyy, Pzz
   assume temperature has already been computed
------------------------------------------------------------------------- */

double ComputePressure::compute_scalar()
{
  if (dimension == 3) {
    inv_volume = 1.0 / (domain->xprd * domain->yprd * domain->zprd);
    virial_compute(3);
    scalar = (temperature->dof * boltz * temperature->scalar + 
	      virial[0] + virial[1] + virial[2]) / 3.0 * inv_volume * nktv2p;
  } else {
    inv_volume = 1.0 / (domain->xprd * domain->yprd);
    virial_compute(2);
    scalar = (temperature->dof * boltz * temperature->scalar + 
	      virial[0] + virial[1]) / 2.0 * inv_volume * nktv2p;
  }

  return scalar;
}

/* ----------------------------------------------------------------------
   compute pressure tensor
   assume KE tensor has already been computed
------------------------------------------------------------------------- */

void ComputePressure::compute_vector()
{
  if (dimension == 3) {
    inv_volume = 1.0 / (domain->xprd * domain->yprd * domain->zprd);
    virial_compute(6);
    double *ke_tensor = temperature->vector;
    for (int i = 0; i < 6; i++)
      vector[i] = (ke_tensor[i] + virial[i]) * inv_volume * nktv2p;
  } else {
    inv_volume = 1.0 / (domain->xprd * domain->yprd);
    virial_compute(4);
    double *ke_tensor = temperature->vector;
    vector[0] = (ke_tensor[0] + virial[0]) * inv_volume * nktv2p;
    vector[1] = (ke_tensor[1] + virial[1]) * inv_volume * nktv2p;
    vector[3] = (ke_tensor[3] + virial[3]) * inv_volume * nktv2p;
  }
}

/* ---------------------------------------------------------------------- */

void ComputePressure::virial_compute(int n)
{
  int i,j;
  double v[6],*vcomponent;

  for (i = 0; i < n; i++) v[i] = 0.0;

  // sum contributions to virial from forces and fixes

  for (j = 0; j < nvirial; j++) {
    vcomponent = vptr[j];
    for (i = 0; i < n; i++) v[i] += vcomponent[i];
  }

  // sum virial across procs

  MPI_Allreduce(v,virial,n,MPI_DOUBLE,MPI_SUM,world);

  // KSpace virial contribution is already summed across procs

  if (force->kspace)
    for (i = 0; i < n; i++) virial[i] += kspace_virial[i];

  // LJ long-range tail correction

  if (force->pair && force->pair->tail_flag)
    for (i = 0; i < n; i++) virial[i] += force->pair->ptail * inv_volume;
}
