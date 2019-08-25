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

#include <cstring>
#include <cstdlib>
#include "compute_global_atom.h"
#include "atom.h"
#include "update.h"
#include "domain.h"
#include "modify.h"
#include "fix.h"
#include "force.h"
#include "comm.h"
#include "group.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

enum{COMPUTE,FIX,VARIABLE};
enum{VECTOR,ARRAY};

#define INVOKED_VECTOR 2
#define INVOKED_ARRAY 4
#define INVOKED_PERATOM 8

#define BIG 1.0e20

/* ---------------------------------------------------------------------- */

ComputeGlobalAtom::ComputeGlobalAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  idref(NULL), which(NULL), argindex(NULL), value2index(NULL), ids(NULL),
  indices(NULL), varatom(NULL), vecglobal(NULL)
{
  if (narg < 5) error->all(FLERR,"Illegal compute global/atom command");

  peratom_flag = 1;

  // process index arg

  int iarg = 3;

  if (strncmp(arg[iarg],"c_",2) == 0 ||
      strncmp(arg[iarg],"f_",2) == 0 ||
      strncmp(arg[iarg],"v_",2) == 0) {
    if (arg[iarg][0] == 'c') whichref = COMPUTE;
    else if (arg[iarg][0] == 'f') whichref = FIX;
    else if (arg[iarg][0] == 'v') whichref = VARIABLE;

    int n = strlen(arg[iarg]);
    char *suffix = new char[n];
    strcpy(suffix,&arg[iarg][2]);

    char *ptr = strchr(suffix,'[');
    if (ptr) {
      if (suffix[strlen(suffix)-1] != ']')
        error->all(FLERR,"Illegal compute global/atom command");
      indexref = atoi(ptr+1);
      *ptr = '\0';
    } else indexref = 0;

    n = strlen(suffix) + 1;
    idref = new char[n];
    strcpy(idref,suffix);
    delete [] suffix;
  } else error->all(FLERR,"Illegal compute global/atom command");

  iarg++;

  // expand args if any have wildcard character "*"

  int expand = 0;
  char **earg;
  int nargnew = input->expand_args(narg-iarg,&arg[iarg],1,earg);

  if (earg != &arg[iarg]) expand = 1;
  arg = earg;

  // parse values until one isn't recognized

  which = new int[nargnew];
  argindex = new int[nargnew];
  ids = new char*[nargnew];
  value2index = new int[nargnew];
  nvalues = 0;

  iarg = 0;
  while (iarg < nargnew) {
    ids[nvalues] = NULL;
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
          error->all(FLERR,"Illegal compute global/atom command");
        argindex[nvalues] = atoi(ptr+1);
        *ptr = '\0';
      } else argindex[nvalues] = 0;

      n = strlen(suffix) + 1;
      ids[nvalues] = new char[n];
      strcpy(ids[nvalues],suffix);
      nvalues++;
      delete [] suffix;

    } else error->all(FLERR,"Illegal compute global/atom command");

    iarg++;
  }

  // if wildcard expansion occurred, free earg memory from expand_args()

  if (expand) {
    for (int i = 0; i < nargnew; i++) delete [] earg[i];
    memory->sfree(earg);
  }

  // setup and error check both index arg and values

  if (whichref == COMPUTE) {
    int icompute = modify->find_compute(idref);
    if (icompute < 0)
      error->all(FLERR,"Compute ID for compute global/atom does not exist");

    if (!modify->compute[icompute]->peratom_flag)
      error->all(FLERR,"Compute global/atom compute does not "
                 "calculate a per-atom vector or array");
    if (indexref == 0 &&
        modify->compute[icompute]->size_peratom_cols != 0)
      error->all(FLERR,"Compute global/atom compute does not "
                 "calculate a per-atom vector");
    if (indexref && modify->compute[icompute]->size_peratom_cols == 0)
      error->all(FLERR,"Compute global/atom compute does not "
                 "calculate a per-atom array");
    if (indexref && indexref > modify->compute[icompute]->size_peratom_cols)
      error->all(FLERR,
                 "Compute global/atom compute array is accessed out-of-range");

  } else if (whichref == FIX) {
    int ifix = modify->find_fix(idref);
    if (ifix < 0)
      error->all(FLERR,"Fix ID for compute global/atom does not exist");
    if (!modify->fix[ifix]->peratom_flag)
      error->all(FLERR,"Compute global/atom fix does not "
                 "calculate a per-atom vector or array");
    if (indexref == 0 &&
        modify->fix[ifix]->size_peratom_cols != 0)
      error->all(FLERR,"Compute global/atom fix does not "
                 "calculate a per-atom vector");
    if (indexref && modify->fix[ifix]->size_peratom_cols == 0)
      error->all(FLERR,"Compute global/atom fix does not "
                 "calculate a per-atom array");
    if (indexref && indexref > modify->fix[ifix]->size_peratom_cols)
      error->all(FLERR,
                 "Compute global/atom fix array is accessed out-of-range");

  } else if (whichref == VARIABLE) {
    int ivariable = input->variable->find(idref);
    if (ivariable < 0)
      error->all(FLERR,"Variable name for compute global/atom does not exist");
    if (input->variable->atomstyle(ivariable) == 0)
      error->all(FLERR,"Compute global/atom variable is not "
                 "atom-style variable");
  }

  for (int i = 0; i < nvalues; i++) {
    if (which[i] == COMPUTE) {
      int icompute = modify->find_compute(ids[i]);
      if (icompute < 0)
        error->all(FLERR,"Compute ID for compute global/atom does not exist");
      if (argindex[i] == 0) {
        if (!modify->compute[icompute]->vector_flag)
          error->all(FLERR,"Compute global/atom compute does not "
                     "calculate a global vector");
      } else {
        if (!modify->compute[icompute]->array_flag)
          error->all(FLERR,"Compute global/atom compute does not "
                     "calculate a global array");
        if (argindex[i] > modify->compute[icompute]->size_array_cols)
          error->all(FLERR,"Compute global/atom compute array is "
                     "accessed out-of-range");
      }

    } else if (which[i] == FIX) {
      int ifix = modify->find_fix(ids[i]);
      if (ifix < 0)
        error->all(FLERR,"Fix ID for compute global/atom does not exist");
      if (argindex[i] == 0) {
        if (!modify->fix[ifix]->vector_flag)
          error->all(FLERR,"Compute global/atom fix does not "
                     "calculate a global vector");
      } else {
        if (!modify->fix[ifix]->array_flag)
          error->all(FLERR,"Compute global/atom fix does not "
                     "calculate a global array");
        if (argindex[i] > modify->fix[ifix]->size_array_cols)
          error->all(FLERR,"Compute global/atom fix array is "
                     "accessed out-of-range");
      }

    } else if (which[i] == VARIABLE) {
      int ivariable = input->variable->find(ids[i]);
      if (ivariable < 0)
        error->all(FLERR,"Variable name for compute global/atom "
                   "does not exist");
      if (input->variable->vectorstyle(ivariable) == 0)
        error->all(FLERR,"Compute global/atom variable is not "
                   "vector-style variable");
    }
  }

  // this compute produces either a peratom vector or array

  if (nvalues == 1) size_peratom_cols = 0;
  else size_peratom_cols = nvalues;

  nmax = maxvector = 0;
  indices = NULL;
  varatom = NULL;
  vecglobal = NULL;

  vector_atom = NULL;
  array_atom = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeGlobalAtom::~ComputeGlobalAtom()
{
  delete [] idref;

  delete [] which;
  delete [] argindex;
  for (int m = 0; m < nvalues; m++) delete [] ids[m];
  delete [] ids;
  delete [] value2index;

  memory->destroy(indices);
  memory->destroy(varatom);
  memory->destroy(vecglobal);
  memory->destroy(vector_atom);
  memory->destroy(array_atom);
}

/* ---------------------------------------------------------------------- */

void ComputeGlobalAtom::init()
{
  // set indices of all computes,fixes,variables

  if (whichref == COMPUTE) {
    int icompute = modify->find_compute(idref);
    if (icompute < 0)
      error->all(FLERR,"Compute ID for compute global/atom does not exist");
    ref2index = icompute;
  } else if (whichref == FIX) {
    int ifix = modify->find_fix(idref);
    if (ifix < 0)
      error->all(FLERR,"Fix ID for compute global/atom does not exist");
    ref2index = ifix;
  } else if (whichref == VARIABLE) {
    int ivariable = input->variable->find(idref);
    if (ivariable < 0)
      error->all(FLERR,"Variable name for compute global/atom does not exist");
    ref2index = ivariable;
  }

  for (int m = 0; m < nvalues; m++) {
    if (which[m] == COMPUTE) {
      int icompute = modify->find_compute(ids[m]);
      if (icompute < 0)
        error->all(FLERR,"Compute ID for compute global/atom does not exist");
      value2index[m] = icompute;

    } else if (which[m] == FIX) {
      int ifix = modify->find_fix(ids[m]);
      if (ifix < 0)
        error->all(FLERR,"Fix ID for compute global/atom does not exist");
      value2index[m] = ifix;

    } else if (which[m] == VARIABLE) {
      int ivariable = input->variable->find(ids[m]);
      if (ivariable < 0)
        error->all(FLERR,"Variable name for compute global/atom "
                   "does not exist");
      value2index[m] = ivariable;
    }
  }
}

/* ---------------------------------------------------------------------- */

void ComputeGlobalAtom::compute_peratom()
{
  int i,j;

  invoked_peratom = update->ntimestep;

  // grow indices and output vector or array if necessary

  if (atom->nmax > nmax) {
    nmax = atom->nmax;
    memory->destroy(indices);
    memory->create(indices,nmax,"global/atom:indices");
    if (nvalues == 1) {
      memory->destroy(vector_atom);
      memory->create(vector_atom,nmax,"global/atom:vector_atom");
    } else {
      memory->destroy(array_atom);
      memory->create(array_atom,nmax,nvalues,"global/atom:array_atom");
    }
  }

  // setup current peratom indices
  // integer indices are rounded down from double values
  // indices are decremented from 1 to N -> 0 to N-1

  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  if (whichref == COMPUTE) {
    Compute *compute = modify->compute[ref2index];

    if (!(compute->invoked_flag & INVOKED_PERATOM)) {
      compute->compute_peratom();
      compute->invoked_flag |= INVOKED_PERATOM;
    }

    if (indexref == 0) {
      double *compute_vector = compute->vector_atom;
      for (i = 0; i < nlocal; i++)
        if (mask[i] & groupbit)
          indices[i] = static_cast<int> (compute_vector[i]) - 1;
    } else {
      double **compute_array = compute->array_atom;
      int im1 = indexref - 1;
      for (i = 0; i < nlocal; i++)
        if (mask[i] & groupbit)
          indices[i] = static_cast<int> (compute_array[i][im1]) - 1;
    }

  } else if (whichref == FIX) {
    if (update->ntimestep % modify->fix[ref2index]->peratom_freq)
      error->all(FLERR,"Fix used in compute global/atom not "
                 "computed at compatible time");
    Fix *fix = modify->fix[ref2index];

    if (indexref == 0) {
      double *fix_vector = fix->vector_atom;
      for (i = 0; i < nlocal; i++)
        if (mask[i] & groupbit)
          indices[i] = static_cast<int> (fix_vector[i]) - 1;
    } else {
      double **fix_array = fix->array_atom;
      int im1 = indexref - 1;
      for (i = 0; i < nlocal; i++)
        if (mask[i] & groupbit)
          indices[i] = static_cast<int> (fix_array[i][im1]) - 1;
    }

  } else if (whichref == VARIABLE) {
    if (atom->nmax > nmax) {
      nmax = atom->nmax;
      memory->destroy(varatom);
      memory->create(varatom,nmax,"global/atom:varatom");
    }

    input->variable->compute_atom(ref2index,igroup,varatom,1,0);
    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
        indices[i] = static_cast<int> (varatom[i]) - 1;
  }

  // loop over values to fill output vector or array

  for (int m = 0; m < nvalues; m++) {

    // output = vector

    if (argindex[m] == 0) {
      int vmax;
      double *source;

      if (which[m] == COMPUTE) {
        Compute *compute = modify->compute[value2index[m]];

        if (!(compute->invoked_flag & INVOKED_VECTOR)) {
          compute->compute_vector();
          compute->invoked_flag |= INVOKED_VECTOR;
        }

        source = compute->vector;
        vmax = compute->size_vector;

      } else if (which[m] == FIX) {
        if (update->ntimestep % modify->fix[value2index[m]]->peratom_freq)
          error->all(FLERR,"Fix used in compute global/atom not "
                     "computed at compatible time");
        Fix *fix = modify->fix[value2index[m]];
        vmax = fix->size_vector;

        if (vmax > maxvector) {
          maxvector = vmax;
          memory->destroy(vecglobal);
          memory->create(vecglobal,maxvector,"global/atom:vecglobal");
        }

        for (i = 0; i < vmax; i++)
          vecglobal[i] = fix->compute_vector(i);

        source = vecglobal;

      } else if (which[m] == VARIABLE) {
        vmax = input->variable->compute_vector(value2index[m],&source);
      }

      if (nvalues == 1) {
        for (i = 0; i < nlocal; i++) {
          vector_atom[i] = 0.0;
          if (mask[i] & groupbit) {
            j = indices[i];
            if (j >= 0 && j < vmax) vector_atom[i] = source[j];
          }
        }
      } else {
        for (i = 0; i < nlocal; i++) {
          array_atom[i][m] = 0.0;
          if (mask[i] & groupbit) {
            j = indices[i];
            if (j >= 0 && j < vmax) array_atom[i][m] = source[j];
          }
        }
      }

    // output = array

    } else {
      int vmax;
      double *source;
      int col = argindex[m] - 1;

      if (which[m] == COMPUTE) {
        Compute *compute = modify->compute[value2index[m]];

        if (!(compute->invoked_flag & INVOKED_ARRAY)) {
          compute->compute_array();
          compute->invoked_flag |= INVOKED_ARRAY;
        }

        double **compute_array = compute->array;
        vmax = compute->size_array_rows;

        if (vmax > maxvector) {
          maxvector = vmax;
          memory->destroy(vecglobal);
          memory->create(vecglobal,maxvector,"global/atom:vecglobal");
        }

        for (i = 0; i < vmax; i++)
          vecglobal[i] = compute_array[i][col];

        source = vecglobal;

      } else if (which[m] == FIX) {
        if (update->ntimestep % modify->fix[value2index[m]]->peratom_freq)
          error->all(FLERR,"Fix used in compute global/atom not "
                     "computed at compatible time");
        Fix *fix = modify->fix[value2index[m]];
        vmax = fix->size_array_rows;

        if (vmax > maxvector) {
          maxvector = vmax;
          memory->destroy(vecglobal);
          memory->create(vecglobal,maxvector,"global/atom:vecglobal");
        }

        for (i = 0; i < vmax; i++)
          vecglobal[i] = fix->compute_array(i,col);

        source = vecglobal;

      } else if (which[m] == VARIABLE) {
        vmax = input->variable->compute_vector(value2index[m],&source);
      }

      if (nvalues == 1) {
        for (i = 0; i < nlocal; i++) {
          vector_atom[i] = 0.0;
          if (mask[i] & groupbit) {
            j = indices[i];
            if (j >= 0 && j < vmax) vector_atom[i] = source[j];
          }
        }
      } else {
        for (i = 0; i < nlocal; i++) {
          array_atom[i][m] = 0.0;
          if (mask[i] & groupbit) {
            j = indices[i];
            if (j >= 0 && j < vmax) array_atom[i][m] = source[j];
          }
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeGlobalAtom::memory_usage()
{
  double bytes = nmax*nvalues * sizeof(double);
  bytes += nmax * sizeof(int);                    // indices
  if (varatom) bytes += nmax * sizeof(double);    // varatom
  bytes += maxvector * sizeof(double);            // vecglobal
  return bytes;
}
