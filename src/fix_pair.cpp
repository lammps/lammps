// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "fix_pair.h"

#include "atom.h"
#include "dump.h"
#include "error.h"
#include "force.h"
#include "fix.h"
#include "input.h"
#include "memory.h"
#include "pair.h"
#include "output.h"
#include "variable.h"
#include "update.h"
#include "variable.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixPair::FixPair(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 7) error->all(FLERR,"Illegal fix pair command");

  nevery = utils::inumeric(FLERR,arg[3],false,lmp);
  if (nevery < 1) error->all(FLERR,"Illegal fix pair command");

  int n = strlen(arg[4]) + 1;
  pairname = new char[n];
  strcpy(pairname,arg[4]);
  pstyle = force->pair_match(pairname,1,0);
  if (pstyle == nullptr) error->all(FLERR,"Fix pair pair style not found");

  nfield = (narg-5) / 2;
  fieldname = new char*[nfield];
  trigger = new int[nfield];

  nfield = 0;
  int iarg = 5;

  while (iarg < narg) {
    if (iarg+2 > narg) error->all(FLERR,"Illegal fix pair command");
    n = strlen(arg[iarg]) + 1;
    fieldname[nfield] = new char[n];
    strcpy(fieldname[nfield],arg[iarg]);
    int flag = utils::inumeric(FLERR,arg[iarg+1],true,lmp);
    if (flag == 0) trigger[nfield] = 0;
    else if (flag == 1) trigger[nfield] = 1;
    else error->all(FLERR,"Illegal fix pair command");
    nfield++;
    iarg += 2;
  }

  // set trigger names = fieldname + "_flag"

  triggername = new char*[nfield];

  for (int ifield = 0; ifield < nfield; ifield++) {
    if (trigger[ifield] == 0) triggername[ifield] = nullptr;
    else {
      n = strlen(fieldname[ifield]) + 6;
      triggername[ifield] = new char[n];
      sprintf(triggername[ifield],"%s_flag",fieldname[ifield]);
    }
  }

  // extract all fields just to get number of per-atom values, ptrs may be NULL
  // also check if triggername can be extracted as a scalar value

  triggerptr = new int*[nfield];

  int columns,dim;
  ncols = 0;

  for (int ifield = 0; ifield < nfield; ifield++) {
    void *tmp = pstyle->extract(fieldname[ifield],columns);
    //if (!tmp) 
    //  error->all(FLERR,"Fix pair pair style cannot extract {}",fieldname[ifield]);
    if (columns) ncols += columns;
    else ncols++;
    if (trigger[ifield]) {
      triggerptr[ifield] = (int *) pstyle->extract(triggername[ifield],dim);
      if (!triggerptr[ifield]) 
        error->all(FLERR,"Fix pair pair style cannot extract {}",
                   triggername[ifield]);
      if (dim) error->all(FLERR,"Fix pair pair style {} is not a scalar",
                          triggername[ifield]);
    }
  }

  // settings
  // freq = 1 since vector/array always have values
  // allows user to specify nevery = 0 for a dump
  //   which this fix outputs whenever it wants

  peratom_flag = 1;
  if (ncols == 1) size_peratom_cols = 0;
  else size_peratom_cols = ncols;
  peratom_freq = nevery;

  // perform initial allocation of atom-based array
  // register with Atom class

  vector = nullptr;
  array = nullptr;
  grow_arrays(atom->nmax);
  atom->add_callback(Atom::GROW);

  // zero the vector/array since dump may access it on timestep 0
  // zero the vector/array since a variable may access it before first run

  int nlocal = atom->nlocal;

  if (ncols == 0) {
    for (int i = 0; i < nlocal; i++)
      vector[i] = 0.0;
  } else {
    for (int i = 0; i < nlocal; i++)
      for (int m = 0; m < ncols; m++)
        array[i][m] = 0.0;
  }
}

/* ---------------------------------------------------------------------- */

FixPair::~FixPair()
{
  // unregister callbacks to this fix from Atom class

  atom->delete_callback(id,Atom::GROW);

  delete [] pairname;
  for (int ifield = 0; ifield < nfield; ifield++) {
    delete [] fieldname[ifield];
    delete [] triggername[ifield];
  }

  delete [] fieldname;
  delete [] trigger;
  delete [] triggername;
  delete [] triggerptr;
  
  if (ncols == 0) memory->destroy(vector);
  else memory->destroy(array);
}

/* ---------------------------------------------------------------------- */

int FixPair::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixPair::init()
{
  // make sure pair style still exists

  pstyle = force->pair_match(pairname,1,0);
  if (pstyle == nullptr) error->all(FLERR,"Fix pair pair style not found");
}

/* ----------------------------------------------------------------------
   trigger pair style computation on steps which are multiples of Nevery
------------------------------------------------------------------------- */

void FixPair::pre_force(int /*vflag*/)
{
  if (update->ntimestep % nevery) return;

  // set pair style triggers

  for (int ifield = 0; ifield < nfield; ifield++)
    if (trigger[ifield]) *(triggerptr[ifield]) = 1;
}

/* ----------------------------------------------------------------------
   extract results from pair style
------------------------------------------------------------------------- */

void FixPair::post_force(int /*vflag*/)
{
  if (update->ntimestep % nevery) return;

  // extract pair style fields one by one 
  // store their values in this fix

  int nlocal = atom->nlocal;

  int icol = 0;
  int columns;

  for (int ifield = 0; ifield < nfield; ifield++) {
    void *pvoid = pstyle->extract(fieldname[ifield],columns);
    if (pvoid == nullptr) 
      error->all(FLERR,"Fix pair pair style cannot extract {}",fieldname[ifield]);

    if (columns == 0) {
      double *pvector = (double *) pvoid;
      if (ncols == 0) { 
        for (int i = 0; i < nlocal; i++)
          vector[i] = pvector[i];
      } else {
        for (int i = 0; i < nlocal; i++)
          array[i][icol] = pvector[i];
      }
      icol++;

    } else {
      double **parray = (double **) pvoid;
      for (int i = 0; i < nlocal; i++)
        for (int m = 0; m < columns; m++) {
          array[i][icol] = parray[i][m];
          icol++;
        }
    }
  }

  // unset pair style triggers

  for (int ifield = 0; ifield < nfield; ifield++)
    if (trigger[ifield]) *(triggerptr[ifield]) = 0;
}

/* ----------------------------------------------------------------------
   allocate atom-based vector or array
------------------------------------------------------------------------- */

void FixPair::grow_arrays(int nmax)
{
  if (ncols == 0) {
    memory->grow(vector,nmax,"store/state:vector");
    vector_atom = vector;
  } else {
    memory->grow(array,nmax,ncols,"store/state:array");
    array_atom = array;
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based vector or array
------------------------------------------------------------------------- */

double FixPair::memory_usage()
{
  double bytes = 0.0;
  if (ncols == 0) bytes += (double)atom->nmax * sizeof(double);
  else bytes += (double)atom->nmax*ncols * sizeof(double);
  return bytes;
}
