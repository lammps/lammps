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
#include <cstdlib>
#include <cstring>
#include "fix_ave_deviation.h"
#include "atom.h"
#include "domain.h"
#include "update.h"
#include "modify.h"
#include "compute.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"
#include "force.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#define INVOKED_PERATOM 8

/* ---------------------------------------------------------------------- */

FixAveDeviation::FixAveDeviation(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), nvalues(3), array(NULL)
{
  if (narg < 6) error->all(FLERR,"Illegal fix ave/deviation command");

  if (!atom->pafi_flag) error->all(FLERR,"Fix ave/deviation requires atom_style pafi or pafipath");

  nevery = force->inumeric(FLERR,arg[3]);
  nrepeat = force->inumeric(FLERR,arg[4]);
  peratom_freq = force->inumeric(FLERR,arg[5]);

  // this fix produces a per-atom array
  nvalues = 3;
  peratom_flag = 1;
  array_flag = 1;
  size_peratom_cols = 3;

  // perform initial allocation of atom-based array
  // register with Atom class
  grow_arrays(atom->nmax);
  atom->add_callback(0);

  // set to initial deviation

  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  double **x = atom->x;
  double **path = atom->path;
  double deviation[3] = {0.,0.,0.};
  for (int i = 0; i < nlocal; i++) if (mask[i] & groupbit) {
    for(int j=0;j<3;j++) deviation[j] = x[i][j]-path[i][j]; // x-path
    domain->minimum_image(deviation);
    for(int j=0;j<3;j++) array[i][j] = deviation[j];
  }

  // nvalid = next step on which end_of_step does something
  // add nvalid to all computes that store invocation times
  // since don't know a priori which are invoked by this fix
  // once in end_of_step() can set timestep for ones actually invoked

  irepeat = 0;
  nvalid_last = -1;
  nvalid = nextvalid();
  modify->addstep_compute_all(nvalid);
}

/* ---------------------------------------------------------------------- */

FixAveDeviation::~FixAveDeviation()
{
  // unregister callback to this fix from Atom class

  atom->delete_callback(id,0);
  memory->destroy(array);
}

/* ---------------------------------------------------------------------- */

int FixAveDeviation::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixAveDeviation::init()
{
  // need to reset nvalid if nvalid < ntimestep b/c minimize was performed
  if (nvalid < update->ntimestep) {
    irepeat = 0;
    nvalid = nextvalid();
    modify->addstep_compute_all(nvalid);
  }
}

/* ----------------------------------------------------------------------
   only does something if nvalid = current timestep
------------------------------------------------------------------------- */

void FixAveDeviation::setup(int vflag)
{
  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixAveDeviation::end_of_step()
{
  int i,j,m,n;

  // skip if not step which requires doing something
  // error check if timestep was reset in an invalid manner

  bigint ntimestep = update->ntimestep;
  if (ntimestep < nvalid_last || ntimestep > nvalid)
    error->all(FLERR,"Invalid timestep reset for fix ave/deviation");
  if (ntimestep != nvalid) return;
  nvalid_last = nvalid;

  // zero if first step

  int nlocal = atom->nlocal;

  if (irepeat == 0)
    for (i = 0; i < nlocal; i++)
      for (m = 0; m < nvalues; m++)
        array[i][m] = 0.0;

  // accumulate results of attributes,computes,fixes,variables to local copy
  // compute/fix/variable may invoke computes so wrap with clear/add

  //modify->clearstep_compute();

  int *mask = atom->mask;
  double **x = atom->x;
  double **path = atom->path;
  double deviation[3] = {0.,0.,0.};
  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      for(int j=0;j<3;j++) deviation[j] = x[i][j]-path[i][j]; // x-path
      domain->minimum_image(deviation);
      for(int j=0;j<3;j++) array[i][j] += deviation[j];
    }

  // done if irepeat < nrepeat
  // else reset irepeat and nvalid

  irepeat++;
  if (irepeat < nrepeat) {
    nvalid += nevery;
    //modify->addstep_compute(nvalid);
    return;
  }

  irepeat = 0;
  nvalid = ntimestep+peratom_freq - (nrepeat-1)*nevery;
  //modify->addstep_compute(nvalid);

  if (array == NULL) return;

  // average the final result for the Nfreq timestep
  double repeat = nrepeat;
  for (i = 0; i < nlocal; i++)
    for (m = 0; m < nvalues; m++)
      array[i][m] /= repeat;
}

/* ----------------------------------------------------------------------
  return array value - need to MPI_reduce??
------------------------------------------------------------------------- */

double FixAveDeviation::compute_array(int i,int j) {
  return array[i][j];
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixAveDeviation::memory_usage()
{
  double bytes;
  bytes = atom->nmax*nvalues * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   allocate atom-based array
------------------------------------------------------------------------- */

void FixAveDeviation::grow_arrays(int nmax)
{
  memory->grow(array,nmax,nvalues,"fix_ave/deviation:array");
  array_atom = array;
  if (array) vector_atom = array[0];
  else vector_atom = NULL;
}

/* ----------------------------------------------------------------------
   copy values within local atom-based array
------------------------------------------------------------------------- */

void FixAveDeviation::copy_arrays(int i, int j, int delflag)
{
  for (int m = 0; m < nvalues; m++)
    array[j][m] = array[i][m];
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

int FixAveDeviation::pack_exchange(int i, double *buf)
{
  for (int m = 0; m < nvalues; m++) buf[m] = array[i][m];
  return nvalues;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
------------------------------------------------------------------------- */

int FixAveDeviation::unpack_exchange(int nlocal, double *buf)
{
  for (int m = 0; m < nvalues; m++) array[nlocal][m] = buf[m];
  return nvalues;
}

/* ----------------------------------------------------------------------
   calculate nvalid = next step on which end_of_step does something
   can be this timestep if multiple of nfreq and nrepeat = 1
   else backup from next multiple of nfreq
------------------------------------------------------------------------- */

bigint FixAveDeviation::nextvalid()
{
  bigint nvalid = (update->ntimestep/peratom_freq)*peratom_freq + peratom_freq;
  if (nvalid-peratom_freq == update->ntimestep && nrepeat == 1)
    nvalid = update->ntimestep;
  else
    nvalid -= (nrepeat-1)*nevery;
  if (nvalid < update->ntimestep) nvalid += peratom_freq;
  return nvalid;
}
