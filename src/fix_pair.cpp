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
  if (narg < 7) utils::missing_cmd_args(FLERR, "fix pair", error);

  nevery = utils::inumeric(FLERR,arg[3],false,lmp);
  if (nevery < 1) error->all(FLERR,"Illegal fix pair every value: {}", nevery);

  pairname = utils::strdup(arg[4]);
  pstyle = force->pair_match(pairname,1,0);
  if (pstyle == nullptr) error->all(FLERR,"Pair style {} for fix pair not found", pairname);

  nfield = (narg-5) / 2;
  fieldname = new char*[nfield];
  trigger = new int[nfield];

  nfield = 0;
  int iarg = 5;

  while (iarg < narg) {
    if (iarg+2 > narg) utils::missing_cmd_args(FLERR, fmt::format("fix pair {}", arg[iarg]), error);
    fieldname[nfield] = utils::strdup(arg[iarg]);
    int flag = utils::inumeric(FLERR,arg[iarg+1],true,lmp);
    if (flag == 0) trigger[nfield] = 0;
    else if (flag == 1) trigger[nfield] = 1;
    else error->all(FLERR,"Illegal fix pair {} command flag: {}", arg[iarg], arg[iarg+1]);
    nfield++;
    iarg += 2;
  }

  // set trigger names = fieldname + "_flag"

  triggername = new char*[nfield];

  for (int ifield = 0; ifield < nfield; ifield++) {
    if (trigger[ifield] == 0) triggername[ifield] = nullptr;
    else triggername[ifield] = utils::strdup(fmt::format("{}_flag", fieldname[ifield]));
  }

  // extract all fields just to get number of per-atom values
  //   returned data ptr may be NULL, if pair style has not allocated field yet
  //   check for recognized field cannot be done until post_force()
  // also check if triggername can be extracted as a scalar value

  triggerptr = new int*[nfield];

  ncols = 0;

  for (int ifield = 0; ifield < nfield; ifield++) {
    int columns = 0;         // set in case fieldname not recognized by pstyle
    void *pvoid = pstyle->extract_peratom(fieldname[ifield],columns);
    if (columns) ncols += columns;
    else ncols++;
    if (trigger[ifield]) {
      int dim;
      triggerptr[ifield] = (int *) pstyle->extract(triggername[ifield],dim);
      if (!triggerptr[ifield])
        error->all(FLERR,"Fix pair pair style cannot extract {}", triggername[ifield]);
      if (dim)
        error->all(FLERR,"Fix pair pair style {} trigger {} is not a scalar",
                   pairname, triggername[ifield]);
    }
  }

  // if set peratom_freq = Nevery, then cannot access the per-atom
  //   values as part of thermo output during minimiziation
  //   at different frequency or on last step of minimization
  // instead set peratom_freq = 1
  //   ok, since vector/array always have values
  //   but requires the vector/array be persisted between Nevery steps
  //     since it may be accessed

  peratom_flag = 1;
  if (ncols == 1) size_peratom_cols = 0;
  else size_peratom_cols = ncols;
  peratom_freq = 1;

  // perform initial allocation of atom-based array
  // register with Atom class

  vector = nullptr;
  array = nullptr;
  grow_arrays(atom->nmax);
  atom->add_callback(Atom::GROW);

  // zero the vector/array since dump may access it on timestep 0
  // zero the vector/array since a variable may access it before first run

  int nlocal = atom->nlocal;

  if (ncols == 1) {
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

  delete[] pairname;
  for (int ifield = 0; ifield < nfield; ifield++) {
    delete[] fieldname[ifield];
    delete[] triggername[ifield];
  }

  delete[] fieldname;
  delete[] trigger;
  delete[] triggername;
  delete[] triggerptr;

  if (ncols == 1) memory->destroy(vector);
  else memory->destroy(array);
}

/* ---------------------------------------------------------------------- */

int FixPair::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  mask |= MIN_PRE_FORCE;
  mask |= POST_FORCE;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixPair::init()
{
  // insure pair style still exists

  pstyle = force->pair_match(pairname,1,0);
  if (pstyle == nullptr) error->all(FLERR,"Pair style {} for fix pair not found", pairname);
}

/* ---------------------------------------------------------------------- */

void FixPair::setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixPair::setup_pre_force(int vflag)
{
  pre_force(vflag);
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

/* ---------------------------------------------------------------------- */

void FixPair::min_pre_force(int vflag)
{
  pre_force(vflag);
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
    void *pvoid = pstyle->extract_peratom(fieldname[ifield],columns);
    if (pvoid == nullptr)
      error->all(FLERR,"Fix pair pair style cannot extract {}",fieldname[ifield]);

    if (columns == 0) {
      double *pvector = (double *) pvoid;
      if (ncols == 1) {
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

/* ---------------------------------------------------------------------- */

void FixPair::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   allocate atom-based vector or array
------------------------------------------------------------------------- */

void FixPair::grow_arrays(int nmax)
{
  if (ncols == 1) {
    memory->grow(vector,nmax,"store/state:vector");
    vector_atom = vector;
  } else {
    memory->grow(array,nmax,ncols,"store/state:array");
    array_atom = array;
  }
}

/* ----------------------------------------------------------------------
   copy values within local atom-based array
------------------------------------------------------------------------- */

void FixPair::copy_arrays(int i, int j, int /*delflag*/)
{
  if (ncols == 1) {
    vector[j] = vector[i];
  } else {
    for (int m = 0; m < ncols; m++)
      array[j][m] = array[i][m];
  }
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

int FixPair::pack_exchange(int i, double *buf)
{
  if (ncols == 1) {
    buf[0] = vector[i];
  } else {
    for (int m = 0; m < ncols; m++)
      buf[m] = array[i][m];
  }

  return ncols;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
------------------------------------------------------------------------- */

int FixPair::unpack_exchange(int nlocal, double *buf)
{
  if (ncols == 1) {
    vector[nlocal] = buf[0];
  } else {
    for (int m = 0; m < ncols; m++)
      array[nlocal][m] = buf[m];
  }

  return ncols;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based vector or array
------------------------------------------------------------------------- */

double FixPair::memory_usage()
{
  double bytes = 0.0;
  if (ncols == 1) bytes += (double)atom->nmax * sizeof(double);
  else bytes += (double)atom->nmax*ncols * sizeof(double);
  return bytes;
}
