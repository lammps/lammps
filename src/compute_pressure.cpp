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

  int n = strlen(arg[3]) + 1;
  id_pre = new char[n];
  strcpy(id_pre,arg[3]);

  int icompute = modify->find_compute(id_pre);

  if (icompute < 0) error->all("Could not find compute pressure temp ID");
  if (modify->compute[icompute]->tempflag == 0)
    error->all("Compute pressure temp ID does not compute temperature");

  scalar_flag = vector_flag = 1;
  size_vector = 6;
  extensive = 0;
  pressflag = 1;

  vector = new double[6];
}

/* ---------------------------------------------------------------------- */

ComputePressure::~ComputePressure()
{
  delete [] id_pre;
  delete [] vector;
}

/* ---------------------------------------------------------------------- */

void ComputePressure::init()
{
  boltz = force->boltz;
  nktv2p = force->nktv2p;

  // set temperature used by pressure

  int icompute = modify->find_compute(id_pre);
  if (icompute < 0) error->all("Could not find compute pressure temp ID");
  temperature = modify->compute[icompute];

  // set flags/ptrs for all contributions to virial
    
  pairflag = bondflag = angleflag = dihedralflag = improperflag = 0;
  kspaceflag = 0;
  shakeflag = bodyflag = rigidflag = poemsflag = 0;

  if (force->pair) {
    pairflag = 1;
    pair_virial = force->pair->virial;
  }
  if (atom->molecular) {
    if (force->bond) {
      bondflag = 1;
      bond_virial = force->bond->virial;
    }
    if (force->angle) {
      angleflag = 1;
      angle_virial = force->angle->virial;
    }
    if (force->dihedral) {
      dihedralflag = 1;
      dihedral_virial = force->dihedral->virial;
    }
    if (force->improper) {
      improperflag = 1;
      improper_virial = force->improper->virial;
    }
  }

  if (force->kspace) {
    kspaceflag = 1;
    kspace_virial = force->kspace->virial;
  }

  for (int i = 0; i < modify->nfix; i++) {
    if (strcmp(modify->fix[i]->style,"shake") == 0) {
      shakeflag = 1;
      shake_virial = modify->fix[i]->virial;
    }
    if (strcmp(modify->fix[i]->style,"rigid") == 0 ||
	strcmp(modify->fix[i]->style,"poems") == 0) {
      bodyflag = 1;
      if (strcmp(modify->fix[i]->style,"rigid") == 0) {
	rigidflag = 1;
	rigid_virial = modify->fix[i]->virial;
      } else {
	poemsflag = 1;
	poems_virial = modify->fix[i]->virial;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   compute total pressure, averaged over Pxx, Pyy, Pzz
   assume temperature has already been computed
------------------------------------------------------------------------- */

double ComputePressure::compute_scalar()
{
  inv_volume = 1.0 / (domain->xprd * domain->yprd * domain->zprd);
  virial_compute(3);
  scalar = (temperature->dof * boltz * temperature->scalar + 
	    virial[0] + virial[1] + virial[2]) / 3.0 * inv_volume * nktv2p;
  return scalar;
}

/* ----------------------------------------------------------------------
   compute pressure tensor
   assume KE tensor has already been computed
------------------------------------------------------------------------- */

void ComputePressure::compute_vector()
{
  inv_volume = 1.0 / (domain->xprd * domain->yprd * domain->zprd);
  virial_compute(6);
  double *ke_tensor = temperature->vector;
  for (int i = 0; i < 6; i++)
    vector[i] = (ke_tensor[i] + virial[i]) * inv_volume * nktv2p;
}

/* ---------------------------------------------------------------------- */

void ComputePressure::virial_compute(int n)
{
  int i;
  double v[6];

  for (i = 0; i < n; i++) v[i] = 0.0;

  // sum contributions to virial from various forces and fixes

  if (pairflag)
    for (i = 0; i < n; i++) v[i] += pair_virial[i];

  if (atom->molecular) {
    if (bondflag)
      for (i = 0; i < n; i++) v[i] += bond_virial[i];
    if (angleflag)
      for (i = 0; i < n; i++) v[i] += angle_virial[i];
    if (dihedralflag)
      for (i = 0; i < n; i++) v[i] += dihedral_virial[i];
    if (improperflag)
      for (i = 0; i < n; i++) v[i] += improper_virial[i];
    if (shakeflag)
      for (i = 0; i < n; i++) v[i] += shake_virial[i];
  }

  if (bodyflag) {
    if (rigidflag) for (i = 0; i < n; i++) v[i] += rigid_virial[i];
    if (poemsflag) for (i = 0; i < n; i++) v[i] += poems_virial[i];
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
