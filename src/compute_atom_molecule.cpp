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

#include "string.h"
#include "stdlib.h"
#include "compute_atom_molecule.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "compute.h"
#include "fix.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

enum{COMPUTE,FIX,VARIABLE};

#define INVOKED_PERATOM 8

/* ---------------------------------------------------------------------- */

ComputeAtomMolecule::
ComputeAtomMolecule(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal compute atom/molecule command");

  if (atom->molecular == 0)
    error->all(FLERR,"Compute atom/molecule requires molecular atom style");

  // parse args

  which = new int[narg-3];
  argindex = new int[narg-3];
  ids = new char*[narg-3];
  value2index = new int[narg-3];
  nvalues = 0;

  int iarg = 3;
  while (iarg < narg) {
    if (strncmp(arg[iarg],"c_",2) == 0 ||
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
          error->all(FLERR,"Illegal compute reduce command");
        argindex[nvalues] = atoi(ptr+1);
        *ptr = '\0';
      } else argindex[nvalues] = 0;

      n = strlen(suffix) + 1;
      ids[nvalues] = new char[n];
      strcpy(ids[nvalues],suffix);
      nvalues++;
      delete [] suffix;
    } else error->all(FLERR,"Illegal compute atom/molecule command");

    iarg++;
  }

  // setup and error check

  for (int i = 0; i < nvalues; i++) {
    if (which[i] == COMPUTE) {
      int icompute = modify->find_compute(ids[i]);
      if (icompute < 0)
        error->all(FLERR,"Compute ID for compute atom/molecule does not exist");
      if (modify->compute[icompute]->peratom_flag == 0)
        error->all(FLERR,"Compute atom/molecule compute does not "
                   "calculate per-atom values");
      if (argindex[i] == 0 &&
          modify->compute[icompute]->size_peratom_cols != 0)
        error->all(FLERR,"Compute atom/molecule compute does not "
                   "calculate a per-atom vector");
      if (argindex[i] && modify->compute[icompute]->size_peratom_cols == 0)
        error->all(FLERR,"Compute atom/molecule compute does not "
                   "calculate a per-atom array");
      if (argindex[i] &&
          argindex[i] > modify->compute[icompute]->size_peratom_cols)
        error->all(FLERR,"Compute atom/molecule compute array is "
                   "accessed out-of-range");

    } else if (which[i] == FIX) {
      int ifix = modify->find_fix(ids[i]);
      if (ifix < 0)
        error->all(FLERR,"Fix ID for compute atom/molecule does not exist");
      if (modify->fix[ifix]->peratom_flag)
        error->all(FLERR,"Compute atom/molecule fix does not "
                   "calculate per-atom values");
      if (argindex[i] == 0 &&
          modify->fix[ifix]->size_peratom_cols != 0)
        error->all(FLERR,"Compute atom/molecule fix does not "
                   "calculate a per-atom vector");
      if (argindex[i] && modify->fix[ifix]->size_peratom_cols == 0)
        error->all(FLERR,"Compute atom/molecule fix does not "
                   "calculate a per-atom array");
      if (argindex[i] &&
          argindex[i] > modify->fix[ifix]->size_peratom_cols)
        error->all(FLERR,
                   "Compute atom/molecule fix array is accessed out-of-range");

    } else if (which[i] == VARIABLE) {
      int ivariable = input->variable->find(ids[i]);
      if (ivariable < 0)
        error->all(FLERR,
                   "Variable name for compute atom/molecule does not exist");
      if (input->variable->atomstyle(ivariable) == 0)
        error->all(FLERR,"Compute atom/molecule variable is not "
                   "atom-style variable");
    }
  }

  // setup molecule-based data

  nmolecules = molecules_in_group(idlo,idhi);

  vone = vector = NULL;
  aone = array = NULL;

  if (nvalues == 1) {
    memory->create(vone,nmolecules,"atom/molecule:vone");
    memory->create(vector,nmolecules,"atom/molecule:vector");
    vector_flag = 1;
    size_vector = nmolecules;
    extvector = 0;
  } else {
    memory->create(aone,nmolecules,nvalues,"atom/molecule:aone");
    memory->create(array,nmolecules,nvalues,"atom/molecule:array");
    array_flag = 1;
    size_array_rows = nmolecules;
    size_array_cols = nvalues;
    extarray = 0;
  }

  maxatom = 0;
  scratch = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeAtomMolecule::~ComputeAtomMolecule()
{
  delete [] which;
  delete [] argindex;
  for (int m = 0; m < nvalues; m++) delete [] ids[m];
  delete [] ids;
  delete [] value2index;

  memory->destroy(vone);
  memory->destroy(vector);
  memory->destroy(aone);
  memory->destroy(array);
  memory->destroy(scratch);
}

/* ---------------------------------------------------------------------- */

void ComputeAtomMolecule::init()
{
  int ntmp = molecules_in_group(idlo,idhi);
  if (ntmp != nmolecules)
    error->all(FLERR,"Molecule count changed in compute atom/molecule");

  // set indices and check validity of all computes,fixes,variables

  for (int m = 0; m < nvalues; m++) {
    if (which[m] == COMPUTE) {
      int icompute = modify->find_compute(ids[m]);
      if (icompute < 0)
        error->all(FLERR,"Compute ID for compute atom/molecule does not exist");
      value2index[m] = icompute;

    } else if (which[m] == FIX) {
      int ifix = modify->find_fix(ids[m]);
      if (ifix < 0)
        error->all(FLERR,"Fix ID for compute atom/molecule does not exist");
      value2index[m] = ifix;

    } else if (which[m] == VARIABLE) {
      int ivariable = input->variable->find(ids[m]);
      if (ivariable < 0)
        error->all(FLERR,
                   "Variable name for compute atom/molecule does not exist");
      value2index[m] = ivariable;

    } else value2index[m] = -1;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeAtomMolecule::compute_vector()
{
  int i,j,n,imol;

  invoked_vector = update->ntimestep;

  for (n = 0; n < nmolecules; n++) vone[n] = 0.0;
  compute_one(0);

  int *molecule = atom->molecule;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  j = 0;
  for (i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      imol = molecule[i];
      if (molmap) imol = molmap[imol-idlo];
      else imol--;
      vone[imol] += peratom[j];
    }
    j += nstride;
  }

  int me;
  MPI_Comm_rank(world,&me);
  MPI_Allreduce(vone,vector,nmolecules,MPI_DOUBLE,MPI_SUM,world);
}

/* ---------------------------------------------------------------------- */

void ComputeAtomMolecule::compute_array()
{
  int i,j,m,n,imol;

  invoked_array = update->ntimestep;

  int *molecule = atom->molecule;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (m = 0; m < nvalues; m++) {
    for (n = 0; n < nmolecules; n++) aone[n][m] = 0.0;
    compute_one(m);

    j = 0;
    for (i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        imol = molecule[i];
        if (molmap) imol = molmap[imol-idlo];
        else imol--;
        aone[imol][m] += peratom[j];
      }
      j += nstride;
    }
  }

  if (array)
    MPI_Allreduce(&aone[0][0],&array[0][0],nvalues*nmolecules,
                  MPI_DOUBLE,MPI_SUM,world);
}

/* ----------------------------------------------------------------------
   calculate per-atom values for one input M
   invoke the appropriate compute,fix,variable
   reallocate scratch if necessary for per-atom variable scratch space
------------------------------------------------------------------------- */

void ComputeAtomMolecule::compute_one(int m)
{
  int vidx = value2index[m];
  int aidx = argindex[m];

  // invoke compute if not previously invoked

  if (which[m] == COMPUTE) {
    Compute *compute = modify->compute[vidx];

    if (!(compute->invoked_flag & INVOKED_PERATOM)) {
      compute->compute_peratom();
      compute->invoked_flag |= INVOKED_PERATOM;
    }

    if (aidx == 0) {
      peratom = compute->vector_atom;
      nstride = 1;
    } else {
      if (compute->array_atom) peratom = &compute->array_atom[0][aidx-1];
      else peratom = NULL;
      nstride = compute->size_array_cols;
    }

  // access fix fields, check if fix frequency is a match

  } else if (which[m] == FIX) {
    if (update->ntimestep % modify->fix[vidx]->peratom_freq)
      error->all(FLERR,"Fix used in compute atom/molecule not computed "
                 "at compatible time");
    Fix *fix = modify->fix[vidx];

    if (aidx == 0) {
      peratom = fix->vector_atom;
      nstride = 1;
    } else {
      peratom = &fix->array_atom[0][aidx-1];
      nstride = fix->size_array_cols;
    }

  // evaluate atom-style variable

  } else if (which[m] == VARIABLE) {
    if (atom->nlocal > maxatom) {
      maxatom = atom->nmax;
      memory->destroy(scratch);
      memory->create(scratch,maxatom,"atom/molecule:scratch");
      peratom = scratch;
    }

    input->variable->compute_atom(vidx,igroup,peratom,1,0);
    nstride = 1;
  }
}

/* ----------------------------------------------------------------------
   memory usage of local data
------------------------------------------------------------------------- */

double ComputeAtomMolecule::memory_usage()
{
  double bytes = 2*nmolecules*nvalues * sizeof(double);
  if (molmap) bytes += (idhi-idlo+1) * sizeof(int);
  bytes += maxatom * sizeof(double);
  return bytes;
}
