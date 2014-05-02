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

enum{X,V,F,COMPUTE,FIX,VARIABLE};

#define INVOKED_PERATOM 8

/* ---------------------------------------------------------------------- */

FixAveAtom::FixAveAtom(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 7) error->all(FLERR,"Illegal fix ave/atom command");

  nevery = force->inumeric(FLERR,arg[3]);
  nrepeat = force->inumeric(FLERR,arg[4]);
  peratom_freq = force->inumeric(FLERR,arg[5]);

  // parse remaining values

  which = new int[narg-6];
  argindex = new int[narg-6];
  ids = new char*[narg-6];
  value2index = new int[narg-6];
  nvalues = 0;

  int iarg = 6;
  while (iarg < narg) {
    ids[nvalues] = NULL;

    if (strcmp(arg[iarg],"x") == 0) {
      which[nvalues] = X;
      argindex[nvalues++] = 0;
    } else if (strcmp(arg[iarg],"y") == 0) {
      which[nvalues] = X;
      argindex[nvalues++] = 1;
    } else if (strcmp(arg[iarg],"z") == 0) {
      which[nvalues] = X;
      argindex[nvalues++] = 2;

    } else if (strcmp(arg[iarg],"vx") == 0) {
      which[nvalues] = V;
      argindex[nvalues++] = 0;
    } else if (strcmp(arg[iarg],"vy") == 0) {
      which[nvalues] = V;
      argindex[nvalues++] = 1;
    } else if (strcmp(arg[iarg],"vz") == 0) {
      which[nvalues] = V;
      argindex[nvalues++] = 2;

    } else if (strcmp(arg[iarg],"fx") == 0) {
      which[nvalues] = F;
      argindex[nvalues++] = 0;
    } else if (strcmp(arg[iarg],"fy") == 0) {
      which[nvalues] = F;
      argindex[nvalues++] = 1;
    } else if (strcmp(arg[iarg],"fz") == 0) {
      which[nvalues] = F;
      argindex[nvalues++] = 2;

    } else if (strncmp(arg[iarg],"c_",2) == 0 ||
               strncmp(arg[iarg],"f_",2) == 0 ||
               strncmp(arg[iarg],"v_",2) == 0) {
      if (arg[iarg][0] == 'c') which[nvalues] = COMPUTE;
      else if (arg[iarg][0] == 'f') which[nvalues] = FIX;
      else if (arg[iarg][0] == 'v') which[nvalues] = VARIABLE;

      int n = strlen(arg[iarg]);
      char *suffix = new char[n];
      strcpy(suffix,&arg[iarg][2]);

      char *ptr = strchr(suffix,'[');
      if (ptr) {
        if (suffix[strlen(suffix)-1] != ']')
          error->all(FLERR,"Illegal fix ave/atom command");
        argindex[nvalues] = atoi(ptr+1);
        *ptr = '\0';
      } else argindex[nvalues] = 0;

      n = strlen(suffix) + 1;
      ids[nvalues] = new char[n];
      strcpy(ids[nvalues],suffix);
      nvalues++;
      delete [] suffix;

    } else error->all(FLERR,"Illegal fix ave/atom command");

    iarg++;
  }

  // setup and error check
  // for fix inputs, check that fix frequency is acceptable

  if (nevery <= 0 || nrepeat <= 0 || peratom_freq <= 0)
    error->all(FLERR,"Illegal fix ave/atom command");
  if (peratom_freq % nevery || (nrepeat-1)*nevery >= peratom_freq)
    error->all(FLERR,"Illegal fix ave/atom command");

  for (int i = 0; i < nvalues; i++) {
    if (which[i] == COMPUTE) {
      int icompute = modify->find_compute(ids[i]);
      if (icompute < 0)
        error->all(FLERR,"Compute ID for fix ave/atom does not exist");
      if (modify->compute[icompute]->peratom_flag == 0)
        error->all(FLERR,
                   "Fix ave/atom compute does not calculate per-atom values");
      if (argindex[i] == 0 &&
          modify->compute[icompute]->size_peratom_cols != 0)
        error->all(FLERR,"Fix ave/atom compute does not "
                   "calculate a per-atom vector");
      if (argindex[i] && modify->compute[icompute]->size_peratom_cols == 0)
        error->all(FLERR,"Fix ave/atom compute does not "
                   "calculate a per-atom array");
      if (argindex[i] &&
          argindex[i] > modify->compute[icompute]->size_peratom_cols)
        error->all(FLERR,"Fix ave/atom compute array is accessed out-of-range");

    } else if (which[i] == FIX) {
      int ifix = modify->find_fix(ids[i]);
      if (ifix < 0)
        error->all(FLERR,"Fix ID for fix ave/atom does not exist");
      if (modify->fix[ifix]->peratom_flag == 0)
        error->all(FLERR,"Fix ave/atom fix does not calculate per-atom values");
      if (argindex[i] == 0 && modify->fix[ifix]->size_peratom_cols != 0)
        error->all(FLERR,
                   "Fix ave/atom fix does not calculate a per-atom vector");
      if (argindex[i] && modify->fix[ifix]->size_peratom_cols == 0)
        error->all(FLERR,
                   "Fix ave/atom fix does not calculate a per-atom array");
      if (argindex[i] && argindex[i] > modify->fix[ifix]->size_peratom_cols)
        error->all(FLERR,"Fix ave/atom fix array is accessed out-of-range");
      if (nevery % modify->fix[ifix]->peratom_freq)
        error->all(FLERR,
                   "Fix for fix ave/atom not computed at compatible time");

    } else if (which[i] == VARIABLE) {
      int ivariable = input->variable->find(ids[i]);
      if (ivariable < 0)
        error->all(FLERR,"Variable name for fix ave/atom does not exist");
      if (input->variable->atomstyle(ivariable) == 0)
        error->all(FLERR,"Fix ave/atom variable is not atom-style variable");
    }
  }

  // this fix produces either a per-atom vector or array

  peratom_flag = 1;
  if (nvalues == 1) size_peratom_cols = 0;
  else size_peratom_cols = nvalues;

  // perform initial allocation of atom-based array
  // register with Atom class

  array = NULL;
  grow_arrays(atom->nmax);
  atom->add_callback(0);

  // zero the array since dump may access it on timestep 0
  // zero the array since a variable may access it before first run

  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++)
    for (int m = 0; m < nvalues; m++)
      array[i][m] = 0.0;

  // nvalid = next step on which end_of_step does something
  // add nvalid to all computes that store invocation times
  // since don't know a priori which are invoked by this fix
  // once in end_of_step() can set timestep for ones actually invoked

  irepeat = 0;
  nvalid = nextvalid();
  modify->addstep_compute_all(nvalid);
}

/* ---------------------------------------------------------------------- */

FixAveAtom::~FixAveAtom()
{
  // unregister callback to this fix from Atom class

  atom->delete_callback(id,0);

  delete [] which;
  delete [] argindex;
  for (int m = 0; m < nvalues; m++) delete [] ids[m];
  delete [] ids;
  delete [] value2index;

  memory->destroy(array);
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
  // set indices and check validity of all computes,fixes,variables

  for (int m = 0; m < nvalues; m++) {
    if (which[m] == COMPUTE) {
      int icompute = modify->find_compute(ids[m]);
      if (icompute < 0)
        error->all(FLERR,"Compute ID for fix ave/atom does not exist");
      value2index[m] = icompute;

    } else if (which[m] == FIX) {
      int ifix = modify->find_fix(ids[m]);
      if (ifix < 0)
        error->all(FLERR,"Fix ID for fix ave/atom does not exist");
      value2index[m] = ifix;

    } else if (which[m] == VARIABLE) {
      int ivariable = input->variable->find(ids[m]);
      if (ivariable < 0)
        error->all(FLERR,"Variable name for fix ave/atom does not exist");
      value2index[m] = ivariable;

    } else value2index[m] = -1;
  }

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

void FixAveAtom::setup(int vflag)
{
  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixAveAtom::end_of_step()
{
  int i,j,m,n;

  // skip if not step which requires doing something

  bigint ntimestep = update->ntimestep;
  if (ntimestep != nvalid) return;

  // zero if first step

  int nlocal = atom->nlocal;

  if (irepeat == 0)
    for (i = 0; i < nlocal; i++)
      for (m = 0; m < nvalues; m++)
        array[i][m] = 0.0;

  // accumulate results of attributes,computes,fixes,variables to local copy
  // compute/fix/variable may invoke computes so wrap with clear/add

  modify->clearstep_compute();

  int *mask = atom->mask;

  for (m = 0; m < nvalues; m++) {
    n = value2index[m];
    j = argindex[m];

    if (which[m] == X) {
      double **x = atom->x;
      for (i = 0; i < nlocal; i++)
        if (mask[i] & groupbit) array[i][m] += x[i][j];

    } else if (which[m] == V) {
      double **v = atom->v;
      for (i = 0; i < nlocal; i++)
        if (mask[i] & groupbit) array[i][m] += v[i][j];

    } else if (which[m] == F) {
      double **f = atom->f;
      for (i = 0; i < nlocal; i++)
        if (mask[i] & groupbit) array[i][m] += f[i][j];

    // invoke compute if not previously invoked

    } else if (which[m] == COMPUTE) {
      Compute *compute = modify->compute[n];
      if (!(compute->invoked_flag & INVOKED_PERATOM)) {
        compute->compute_peratom();
        compute->invoked_flag |= INVOKED_PERATOM;
      }

      if (j == 0) {
        double *compute_vector = compute->vector_atom;
        for (i = 0; i < nlocal; i++)
          if (mask[i] & groupbit) array[i][m] += compute_vector[i];
      } else {
        int jm1 = j - 1;
        double **compute_array = compute->array_atom;
        for (i = 0; i < nlocal; i++)
          if (mask[i] & groupbit) array[i][m] += compute_array[i][jm1];
      }

    // access fix fields, guaranteed to be ready

    } else if (which[m] == FIX) {
      if (j == 0) {
        double *fix_vector = modify->fix[n]->vector_atom;
        for (i = 0; i < nlocal; i++)
          if (mask[i] & groupbit) array[i][m] += fix_vector[i];
      } else {
        int jm1 = j - 1;
        double **fix_array = modify->fix[n]->array_atom;
        for (i = 0; i < nlocal; i++)
          if (mask[i] & groupbit) array[i][m] += fix_array[i][jm1];
      }

    // evaluate atom-style variable
    // final argument = 1 sums result to array

    } else if (which[m] == VARIABLE) {
      if (array) input->variable->compute_atom(n,igroup,&array[0][m],nvalues,1);
      else input->variable->compute_atom(n,igroup,NULL,nvalues,1);
    }
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
  nvalid = ntimestep+peratom_freq - (nrepeat-1)*nevery;
  modify->addstep_compute(nvalid);

  // average the final result for the Nfreq timestep

  double repeat = nrepeat;
  for (i = 0; i < nlocal; i++)
    for (m = 0; m < nvalues; m++)
      array[i][m] /= repeat;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixAveAtom::memory_usage()
{
  double bytes;
  bytes = atom->nmax*nvalues * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   allocate atom-based array
------------------------------------------------------------------------- */

void FixAveAtom::grow_arrays(int nmax)
{
  memory->grow(array,nmax,nvalues,"fix_ave/atom:array");
  array_atom = array;
  if (array) vector_atom = array[0];
  else vector_atom = NULL;
}

/* ----------------------------------------------------------------------
   copy values within local atom-based array
------------------------------------------------------------------------- */

void FixAveAtom::copy_arrays(int i, int j, int delflag)
{
  for (int m = 0; m < nvalues; m++)
    array[j][m] = array[i][m];
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

int FixAveAtom::pack_exchange(int i, double *buf)
{
  for (int m = 0; m < nvalues; m++) buf[m] = array[i][m];
  return nvalues;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
------------------------------------------------------------------------- */

int FixAveAtom::unpack_exchange(int nlocal, double *buf)
{
  for (int m = 0; m < nvalues; m++) array[nlocal][m] = buf[m];
  return nvalues;
}

/* ----------------------------------------------------------------------
   calculate nvalid = next step on which end_of_step does something
   can be this timestep if multiple of nfreq and nrepeat = 1
   else backup from next multiple of nfreq
------------------------------------------------------------------------- */

bigint FixAveAtom::nextvalid()
{
  bigint nvalid = (update->ntimestep/peratom_freq)*peratom_freq + peratom_freq;
  if (nvalid-peratom_freq == update->ntimestep && nrepeat == 1)
    nvalid = update->ntimestep;
  else
    nvalid -= (nrepeat-1)*nevery;
  if (nvalid < update->ntimestep) nvalid += peratom_freq;
  return nvalid;
}

/* ---------------------------------------------------------------------- */

void FixAveAtom::reset_timestep(bigint ntimestep)
{
  if (ntimestep > nvalid) error->all(FLERR,"Fix ave/atom missed timestep");
}
