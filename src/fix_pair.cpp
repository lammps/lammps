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

enum {PERSTEP, PERTIME, EXTERNAL};     // same as Output

/* ---------------------------------------------------------------------- */

FixPair::FixPair(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 6) error->all(FLERR,"Illegal fix pair command");

  nevery = utils::inumeric(FLERR,arg[3],false,lmp);
  if (nevery < 1) error->all(FLERR,"Illegal fix pair command");
  
  pstyle = force->pair_match(arg[4],1,0);
  if (pstyle == nullptr) error->all(FLERR,"Fix pair pair style not found");

  int n = strlen(arg[5]) + 1;
  field = new char[n];
  strcpy(field,arg[5]);

  // optional args

  ndump = 0;
  dumpID = nullptr;
  dumptrigger = nullptr;
  dumpindex = nullptr;
  varindex = nullptr;

  int iarg = 6;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"dump") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix pair command");
      dumpID = (char **) 
        memory->srealloc(dumpID,(ndump+1)*sizeof(char *),"pair:dumpID");
      dumptrigger = (char **) 
        memory->srealloc(dumptrigger,(ndump+1)*sizeof(char *),"pair:dumptrigger");
      dumpindex = (int *) 
        memory->srealloc(dumpindex,(ndump+1)*sizeof(int),"pair:dumpindex");
      varindex = (int *) 
        memory->srealloc(varindex,(ndump+1)*sizeof(int),"pair:varindex");
      int n = strlen(arg[iarg+1]) + 1;
      dumpID[ndump] = new char[n];
      strcpy(dumpID[ndump],arg[iarg+1]);
      n = strlen(arg[iarg+2]) + 1;
      dumptrigger[ndump] = new char[n];
      strcpy(dumptrigger[ndump],arg[iarg+2]);
      ndump++;
      iarg += 3;
    } else error->all(FLERR,"Illegal fix pair command");
  }

  // register request with pair style
  // pair style returns size of data
  
  //int ierr = pstyle->request(field,ncols);
  //if (ierr < 0) 
  //  error->all(FLERR,"Fix pair pair style does not recognize request");

  // DEBUG for now
  ncols = 3;

  // settings
  // freq = 1 since vector/array always have values
  // allows user to specify nevery = 0 for a dump
  //   which this fix outputs whenever it wants

  peratom_flag = 1;
  peratom_freq = 1;

  if (ncols == 0) size_peratom_cols = 0;
  else size_peratom_cols = ncols;

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

  for (int m = 0; m < ndump; m++) {
    delete [] dumpID[m];
    delete [] dumptrigger[m];
  }

  memory->sfree(dumpID);
  memory->sfree(dumptrigger);
  memory->sfree(dumpindex);
  memory->sfree(varindex);

  if (ncols == 0) memory->destroy(vector);
  else memory->destroy(array);
}

/* ---------------------------------------------------------------------- */

int FixPair::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  mask |= POST_FORCE;    // should this be post-reverse instead ?
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixPair::init()
{
  // setup dump requsets

  for (int m = 0; m < ndump; m++) {
    int idump = output->find_dump(dumpID[m]);
    if (idump < 0) 
      error->all(FLERR, "Dump {} for fix pair does not exist", dumpID[m]);
    if (output->mode_dump[m] != EXTERNAL)
      error->all(FLERR, "Dump {} for fix pair is not external", dumpID[m]);
    dumpindex[m] = idump;

    int ivar = input->variable->find(dumptrigger[m]);
    if (ivar < 0) 
      error->all(FLERR, "Variable {} for fix pair does not exist", 
                 dumptrigger[m]);
    if (!input->variable->equalstyle(ivar))
      error->all(FLERR, "Variable {} for fix pair is not equal style", 
                 dumptrigger[m]);
    varindex[m] = ivar;
  }
}

/* ----------------------------------------------------------------------
   trigger pair style computation on steps which are multiples of Nevery
------------------------------------------------------------------------- */

void FixPair::pre_force(int /*vflag*/)
{
  if (update->ntimestep % nevery) return;

  //pstyle->request_set();
}

/* ----------------------------------------------------------------------
   extract results from pair style
------------------------------------------------------------------------- */

void FixPair::post_force(int /*vflag*/)
{
  if (update->ntimestep % nevery) return;

  // extract the pair style field and store values in this fix
  // NOTE: do this with memcpy ?
  
  if (ncols == 0) {
    int dim;
    double *pvector = (double *) pstyle->extract(field,dim);

    int nlocal = atom->nlocal;
    for (int i = 0; i < nlocal; i++)
      vector[i] = pvector[i];
  } else {
    int dim;
    double **parray = (double **) pstyle->extract(field,dim);

    int nlocal = atom->nlocal;
    for (int i = 0; i < nlocal; i++)
      for (int m = 0; m < ncols; m++)
    // DEBUG
        array[i][m] = atom->x[i][m];
    //        array[i][m] = parray[i][m];
  }

  // trigger any dumps that were requested
  // reset next_dump and next_time_dump, 1 arg for write()

  for (int m = 0; m < ndump; m++) {
    double value = input->variable->compute_equal(varindex[m]);
    if (value == 0.0) continue;

    output->dump[dumpindex[m]]->write();
  }

  // unset the pair style trigger

  //pstyle->request_unset();
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

