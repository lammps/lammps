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

#include <cmath>
#include <cstring>
#include <cctype>
#include "improper_hybrid.h"
#include "atom.h"
#include "neighbor.h"
#include "domain.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define EXTRA 1000

/* ---------------------------------------------------------------------- */

ImproperHybrid::ImproperHybrid(LAMMPS *lmp) : Improper(lmp)
{
  nstyles = 0;
}

/* ---------------------------------------------------------------------- */

ImproperHybrid::~ImproperHybrid()
{
  if (nstyles) {
    for (int i = 0; i < nstyles; i++) delete styles[i];
    delete [] styles;
    for (int i = 0; i < nstyles; i++) delete [] keywords[i];
    delete [] keywords;
  }

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(map);
    delete [] nimproperlist;
    delete [] maximproper;
    for (int i = 0; i < nstyles; i++)
      memory->destroy(improperlist[i]);
    delete [] improperlist;
  }
}

/* ---------------------------------------------------------------------- */

void ImproperHybrid::compute(int eflag, int vflag)
{
  int i,j,m,n;

  // save ptrs to original improperlist

  int nimproperlist_orig = neighbor->nimproperlist;
  int **improperlist_orig = neighbor->improperlist;

  // if this is re-neighbor step, create sub-style improperlists
  // nimproperlist[] = length of each sub-style list
  // realloc sub-style improperlist if necessary
  // load sub-style improperlist with 5 values from original improperlist

  if (neighbor->ago == 0) {
    for (m = 0; m < nstyles; m++) nimproperlist[m] = 0;
    for (i = 0; i < nimproperlist_orig; i++) {
      m = map[improperlist_orig[i][4]];
      if (m >= 0) nimproperlist[m]++;
    }
    for (m = 0; m < nstyles; m++) {
      if (nimproperlist[m] > maximproper[m]) {
        memory->destroy(improperlist[m]);
        maximproper[m] = nimproperlist[m] + EXTRA;
        memory->create(improperlist[m],maximproper[m],5,
                       "improper_hybrid:improperlist");
      }
      nimproperlist[m] = 0;
    }
    for (i = 0; i < nimproperlist_orig; i++) {
      m = map[improperlist_orig[i][4]];
      if (m < 0) continue;
      n = nimproperlist[m];
      improperlist[m][n][0] = improperlist_orig[i][0];
      improperlist[m][n][1] = improperlist_orig[i][1];
      improperlist[m][n][2] = improperlist_orig[i][2];
      improperlist[m][n][3] = improperlist_orig[i][3];
      improperlist[m][n][4] = improperlist_orig[i][4];
      nimproperlist[m]++;
    }
  }

  // call each sub-style's compute function
  // set neighbor->improperlist to sub-style improperlist before call
  // accumulate sub-style global/peratom energy/virial in hybrid

  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = eflag_global = vflag_global = eflag_atom = vflag_atom = 0;

  for (m = 0; m < nstyles; m++) {
    neighbor->nimproperlist = nimproperlist[m];
    neighbor->improperlist = improperlist[m];

    styles[m]->compute(eflag,vflag);

    if (eflag_global) energy += styles[m]->energy;
    if (vflag_global)
      for (n = 0; n < 6; n++) virial[n] += styles[m]->virial[n];
    if (eflag_atom) {
      n = atom->nlocal;
      if (force->newton_bond) n += atom->nghost;
      double *eatom_substyle = styles[m]->eatom;
      for (i = 0; i < n; i++) eatom[i] += eatom_substyle[i];
    }
    if (vflag_atom) {
      n = atom->nlocal;
      if (force->newton_bond) n += atom->nghost;
      double **vatom_substyle = styles[m]->vatom;
      for (i = 0; i < n; i++)
        for (j = 0; j < 6; j++)
          vatom[i][j] += vatom_substyle[i][j];
    }
  }

  // restore ptrs to original improperlist

  neighbor->nimproperlist = nimproperlist_orig;
  neighbor->improperlist = improperlist_orig;
}

/* ---------------------------------------------------------------------- */

void ImproperHybrid::allocate()
{
  allocated = 1;
  int n = atom->nimpropertypes;

  memory->create(map,n+1,"improper:map");
  memory->create(setflag,n+1,"improper:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;

  nimproperlist = new int[nstyles];
  maximproper = new int[nstyles];
  improperlist = new int**[nstyles];
  for (int m = 0; m < nstyles; m++) maximproper[m] = 0;
  for (int m = 0; m < nstyles; m++) improperlist[m] = NULL;
}

/* ---------------------------------------------------------------------- */

void ImproperHybrid::init_style()
{
    for (int i = 0; i < nstyles; i++)
      styles[i]->init_style();
}


/* ----------------------------------------------------------------------
   create one improper style for each arg in list
------------------------------------------------------------------------- */

void ImproperHybrid::settings(int narg, char **arg)
{
  int i,m,istyle;

  if (narg < 1) error->all(FLERR,"Illegal improper_style command");

  // delete old lists, since cannot just change settings

  if (nstyles) {
    for (int i = 0; i < nstyles; i++) delete styles[i];
    delete [] styles;
    for (int i = 0; i < nstyles; i++) delete [] keywords[i];
    delete [] keywords;
  }

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(map);
    delete [] nimproperlist;
    delete [] maximproper;
    for (int i = 0; i < nstyles; i++)
      memory->destroy(improperlist[i]);
    delete [] improperlist;
  }
  allocated = 0;

  // count sub-styles by skipping numeric args
  // one exception is 1st arg of style "table", which is non-numeric word
  // need a better way to skip these exceptions

  nstyles = 0;
  i = 0;
  while (i < narg) {
    if (strcmp(arg[i],"table") == 0) i++;
    i++;
    while (i < narg && !isalpha(arg[i][0])) i++;
    nstyles++;
  }

  // allocate list of sub-styles

  styles = new Improper*[nstyles];
  keywords = new char*[nstyles];

  // allocate each sub-style and call its settings() with subset of args
  // allocate uses suffix, but don't store suffix version in keywords,
  //   else syntax in coeff() will not match
  // define subset of args for a sub-style by skipping numeric args
  // one exception is 1st arg of style "table", which is non-numeric
  // need a better way to skip these exceptions

  int dummy;
  nstyles = 0;
  i = 0;

  while (i < narg) {
    for (m = 0; m < nstyles; m++)
      if (strcmp(arg[i],keywords[m]) == 0)
        error->all(FLERR,"Improper style hybrid cannot use "
                   "same improper style twice");
    if (strcmp(arg[i],"hybrid") == 0)
      error->all(FLERR,
                 "Improper style hybrid cannot have hybrid as an argument");
    if (strcmp(arg[i],"none") == 0)
      error->all(FLERR,"Improper style hybrid cannot have none as an argument");

    styles[nstyles] = force->new_improper(arg[i],1,dummy);
    force->store_style(keywords[nstyles],arg[i],0);

    istyle = i;
    if (strcmp(arg[i],"table") == 0) i++;
    i++;
    while (i < narg && !isalpha(arg[i][0])) i++;
    styles[nstyles]->settings(i-istyle-1,&arg[istyle+1]);
    nstyles++;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one type
---------------------------------------------------------------------- */

void ImproperHybrid::coeff(int narg, char **arg)
{
  if (!allocated) allocate();

  int ilo,ihi;
  force->bounds(FLERR,arg[0],atom->nimpropertypes,ilo,ihi);

  // 2nd arg = improper sub-style name
  // allow for "none" as valid sub-style name

  int m;
  for (m = 0; m < nstyles; m++)
    if (strcmp(arg[1],keywords[m]) == 0) break;

  int none = 0;
  if (m == nstyles) {
    if (strcmp(arg[1],"none") == 0) none = 1;
    else error->all(FLERR,"Improper coeff for hybrid has invalid style");
  }

  // move 1st arg to 2nd arg
  // just copy ptrs, since arg[] points into original input line

  arg[1] = arg[0];

  // invoke sub-style coeff() starting with 1st arg

  if (!none) styles[m]->coeff(narg-1,&arg[1]);

  // set setflag and which type maps to which sub-style
  // if sub-style is none: set hybrid setflag, wipe out map

  for (int i = ilo; i <= ihi; i++) {
    if (none) {
      setflag[i] = 1;
      map[i] = -1;
    } else {
      setflag[i] = styles[m]->setflag[i];
      map[i] = m;
    }
  }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void ImproperHybrid::write_restart(FILE *fp)
{
  fwrite(&nstyles,sizeof(int),1,fp);

  int n;
  for (int m = 0; m < nstyles; m++) {
    n = strlen(keywords[m]) + 1;
    fwrite(&n,sizeof(int),1,fp);
    fwrite(keywords[m],sizeof(char),n,fp);
  }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void ImproperHybrid::read_restart(FILE *fp)
{
  int me = comm->me;
  if (me == 0) fread(&nstyles,sizeof(int),1,fp);
  MPI_Bcast(&nstyles,1,MPI_INT,0,world);
  styles = new Improper*[nstyles];
  keywords = new char*[nstyles];

  allocate();

  int n,dummy;
  for (int m = 0; m < nstyles; m++) {
    if (me == 0) fread(&n,sizeof(int),1,fp);
    MPI_Bcast(&n,1,MPI_INT,0,world);
    keywords[m] = new char[n];
    if (me == 0) fread(keywords[m],sizeof(char),n,fp);
    MPI_Bcast(keywords[m],n,MPI_CHAR,0,world);
    styles[m] = force->new_improper(keywords[m],0,dummy);
  }
}

/* ----------------------------------------------------------------------
   memory usage
------------------------------------------------------------------------- */

double ImproperHybrid::memory_usage()
{
  double bytes = maxeatom * sizeof(double);
  bytes += maxvatom*6 * sizeof(double);
  for (int m = 0; m < nstyles; m++) bytes += maximproper[m]*5 * sizeof(int);
  for (int m = 0; m < nstyles; m++)
    if (styles[m]) bytes += styles[m]->memory_usage();
  return bytes;
}
