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

#include "stdlib.h"
#include "string.h"
#include "fix_ave_atom.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "compute.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixAveAtom::FixAveAtom(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 7) error->all("Illegal fix ave/atom command");

  nevery = atoi(arg[3]);
  nrepeat = atoi(arg[4]);
  peratom_freq = atoi(arg[5]);

  int n = strlen(arg[6]) + 1;
  id_compute = new char[n];
  strcpy(id_compute,arg[6]);

  // setup and error check

  if (nevery <= 0) error->all("Illegal fix ave/atom command");
  if (peratom_freq < nevery || peratom_freq % nevery ||
      (nrepeat-1)*nevery >= peratom_freq)
    error->all("Illegal fix ave/atom command");

  int icompute = modify->find_compute(id_compute);
  if (icompute < 0) error->all("Compute ID for fix ave/atom does not exist");

  if (modify->compute[icompute]->peratom_flag == 0)
    error->all("Fix ave/atom compute does not calculate per-atom info");

  peratom_flag = 1;

  // setup list of computes to call, including pre-computes

  ncompute = 1 + modify->compute[icompute]->npre;
  compute = new Compute*[ncompute];

  // perform initial allocation of atom-based array
  // register with Atom class

  size_peratom = modify->compute[icompute]->size_peratom;
  scalar = NULL;
  vector = NULL;
  grow_arrays(atom->nmax);
  atom->add_callback(0);

  // zero the array since dump may access it on timestep 0

  int nlocal = atom->nlocal;
  if (size_peratom == 0)
    for (int i = 0; i < nlocal; i++) scalar[i] = 0.0;
  else
    for (int i = 0; i < nlocal; i++)
      for (int m = 0; m < size_peratom; m++)
	vector[i][m] = 0.0;

  // nvalid = next step on which end_of_step does something
  // can be this timestep if multiple of peratom_freq and nrepeat = 1
  // else backup from next multiple of peratom_freq

  irepeat = 0;
  nvalid = (update->ntimestep/peratom_freq)*peratom_freq + peratom_freq;
  if (nvalid-peratom_freq == update->ntimestep && nrepeat == 1)
    nvalid = update->ntimestep;
  else
    nvalid -= (nrepeat-1)*nevery;
  if (nvalid < update->ntimestep) nvalid += peratom_freq;

  // must set timestep for all computes that store invocation times
  // since don't know a priori which are invoked by this fix
  // once in end_of_step() can just set timestep for ones actually invoked

  for (int i = 0; i < modify->ncompute; i++)
    if (modify->compute[i]->timeflag) modify->compute[i]->add_step(nvalid);
}

/* ---------------------------------------------------------------------- */

FixAveAtom::~FixAveAtom()
{
  // unregister callback to this fix from Atom class
 
  atom->delete_callback(id,0);

  delete [] id_compute;
  delete [] compute;
  memory->sfree(scalar);
  memory->destroy_2d_double_array(vector);
}

/* ---------------------------------------------------------------------- */

int FixAveAtom::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixAveAtom::init()
{
  // set ptrs to compute and its pre-computes called each end-of-step
  // put pre-computes in list before compute

  int icompute = modify->find_compute(id_compute);
  if (icompute < 0) error->all("Compute ID for fix ave/atom does not exist");
  
  ncompute = 0;
  if (modify->compute[icompute]->npre)
    for (int i = 0; i < modify->compute[icompute]->npre; i++) {
      int ic = modify->find_compute(modify->compute[icompute]->id_pre[i]);
      if (ic < 0)
	error->all("Precompute ID of compute for fix ave/atom does not exist");
      compute[ncompute++] = modify->compute[ic];
    }

  compute[ncompute++] = modify->compute[icompute];
}

/* ---------------------------------------------------------------------- */

void FixAveAtom::end_of_step()
{
  int i,m;

  // skip if not step which requires doing something

  if (update->ntimestep != nvalid) return;

  // zero if first step

  int nlocal = atom->nlocal;

  if (irepeat == 0) {
    if (size_peratom == 0)
      for (i = 0; i < nlocal; i++) scalar[i] = 0.0;
    else
      for (i = 0; i < nlocal; i++)
	for (m = 0; m < size_peratom; m++)
	  vector[i][m] = 0.0;
  }
  
  // accumulate results of compute to local copy
  
  modify->clearstep_compute();
  for (i = 0; i < ncompute; i++) compute[i]->compute_peratom();
  
  int *mask = atom->mask;

  if (size_peratom == 0) {
    double *compute_scalar = compute[ncompute-1]->scalar_atom;
    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) scalar[i] += compute_scalar[i];
  } else {
    double **compute_vector = compute[ncompute-1]->vector_atom;
    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
	for (m = 0; m < size_peratom; m++)
	  vector[i][m] += compute_vector[i][m];
  }

  // done if irepeat < nrepeat
  // else reset irepeat and nvalid

  irepeat++;
  if (irepeat < nrepeat) {
    nvalid += nevery;
    modify->addstep_compute(nvalid);
    return;
  }

  irepeat = 0;
  nvalid = update->ntimestep+peratom_freq - (nrepeat-1)*nevery;
  modify->addstep_compute(nvalid);

  // average the final result for the Nfreq timestep

  double repeat = nrepeat;
  if (size_peratom == 0)
    for (i = 0; i < nlocal; i++)
      scalar[i] /= repeat;
  else
    for (i = 0; i < nlocal; i++)
      for (m = 0; m < size_peratom; m++)
	vector[i][m] /= repeat;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixAveAtom::memory_usage()
{
  double bytes;
  if (size_peratom == 0) bytes = atom->nmax * sizeof(double);
  else bytes = atom->nmax*size_peratom * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   allocate atom-based array
------------------------------------------------------------------------- */

void FixAveAtom::grow_arrays(int nmax)
{
  if (size_peratom == 0) {
    scalar = (double *) memory->srealloc(scalar,nmax*sizeof(double),
					 "fix_ave/atom:scalar");
    scalar_atom = scalar;
  } else {
    vector = memory->grow_2d_double_array(vector,nmax,size_peratom,
					  "fix_ave/atom:vector");
    vector_atom = vector;
  }
}

/* ----------------------------------------------------------------------
   copy values within local atom-based array
------------------------------------------------------------------------- */

void FixAveAtom::copy_arrays(int i, int j)
{
  if (size_peratom == 0)
    scalar[j] = scalar[i];
  else
    for (int m = 0; m <= size_peratom; m++)
      vector[j][m] = vector[i][m];
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

int FixAveAtom::pack_exchange(int i, double *buf)
{
  if (size_peratom == 0) {
    buf[0] = scalar[i];
    return 1;
  }

  for (int m = 0; m <= size_peratom; m++) buf[m] = vector[i][m];
  return size_peratom;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
------------------------------------------------------------------------- */

int FixAveAtom::unpack_exchange(int nlocal, double *buf)
{
  if (size_peratom == 0) {
    scalar[nlocal] = buf[0];
    return 1;
  }

  for (int m = 0; m <= size_peratom; m++) vector[nlocal][m] = buf[m];
  return size_peratom;
}
