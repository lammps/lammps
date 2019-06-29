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

#include <mpi.h>
#include <cstring>
#include <cctype>
#include "dihedral_hybrid.h"
#include "atom.h"
#include "neighbor.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define EXTRA 1000

/* ---------------------------------------------------------------------- */

DihedralHybrid::DihedralHybrid(LAMMPS *lmp) : Dihedral(lmp)
{
  nstyles = 0;
}

/* ---------------------------------------------------------------------- */

DihedralHybrid::~DihedralHybrid()
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
    delete [] ndihedrallist;
    delete [] maxdihedral;
    for (int i = 0; i < nstyles; i++)
      memory->destroy(dihedrallist[i]);
    delete [] dihedrallist;
  }
}

/* ---------------------------------------------------------------------- */

void DihedralHybrid::compute(int eflag, int vflag)
{
  int i,j,m,n;

  // save ptrs to original dihedrallist

  int ndihedrallist_orig = neighbor->ndihedrallist;
  int **dihedrallist_orig = neighbor->dihedrallist;

  // if this is re-neighbor step, create sub-style dihedrallists
  // ndihedrallist[] = length of each sub-style list
  // realloc sub-style dihedrallist if necessary
  // load sub-style dihedrallist with 5 values from original dihedrallist

  if (neighbor->ago == 0) {
    for (m = 0; m < nstyles; m++) ndihedrallist[m] = 0;
    for (i = 0; i < ndihedrallist_orig; i++) {
      m = map[dihedrallist_orig[i][4]];
      if (m >= 0) ndihedrallist[m]++;
    }
    for (m = 0; m < nstyles; m++) {
      if (ndihedrallist[m] > maxdihedral[m]) {
        memory->destroy(dihedrallist[m]);
        maxdihedral[m] = ndihedrallist[m] + EXTRA;
        memory->create(dihedrallist[m],maxdihedral[m],5,
                       "dihedral_hybrid:dihedrallist");
      }
      ndihedrallist[m] = 0;
    }
    for (i = 0; i < ndihedrallist_orig; i++) {
      m = map[dihedrallist_orig[i][4]];
      if (m < 0) continue;
      n = ndihedrallist[m];
      dihedrallist[m][n][0] = dihedrallist_orig[i][0];
      dihedrallist[m][n][1] = dihedrallist_orig[i][1];
      dihedrallist[m][n][2] = dihedrallist_orig[i][2];
      dihedrallist[m][n][3] = dihedrallist_orig[i][3];
      dihedrallist[m][n][4] = dihedrallist_orig[i][4];
      ndihedrallist[m]++;
    }
  }

  // call each sub-style's compute function
  // set neighbor->dihedrallist to sub-style dihedrallist before call
  // accumulate sub-style global/peratom energy/virial in hybrid

  ev_init(eflag,vflag);

  for (m = 0; m < nstyles; m++) {
    neighbor->ndihedrallist = ndihedrallist[m];
    neighbor->dihedrallist = dihedrallist[m];

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

  // restore ptrs to original dihedrallist

  neighbor->ndihedrallist = ndihedrallist_orig;
  neighbor->dihedrallist = dihedrallist_orig;
}

/* ---------------------------------------------------------------------- */

void DihedralHybrid::allocate()
{
  allocated = 1;
  int n = atom->ndihedraltypes;

  memory->create(map,n+1,"dihedral:map");
  memory->create(setflag,n+1,"dihedral:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;

  ndihedrallist = new int[nstyles];
  maxdihedral = new int[nstyles];
  dihedrallist = new int**[nstyles];
  for (int m = 0; m < nstyles; m++) maxdihedral[m] = 0;
  for (int m = 0; m < nstyles; m++) dihedrallist[m] = NULL;
}

/* ----------------------------------------------------------------------
   create one dihedral style for each arg in list
------------------------------------------------------------------------- */

void DihedralHybrid::settings(int narg, char **arg)
{
  int i,m,istyle;

  if (narg < 1) error->all(FLERR,"Illegal dihedral_style command");

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
    delete [] ndihedrallist;
    delete [] maxdihedral;
    for (int i = 0; i < nstyles; i++)
      memory->destroy(dihedrallist[i]);
    delete [] dihedrallist;
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

  styles = new Dihedral*[nstyles];
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
        error->all(FLERR,"Dihedral style hybrid cannot use "
                   "same dihedral style twice");
    if (strcmp(arg[i],"hybrid") == 0)
      error->all(FLERR,
                 "Dihedral style hybrid cannot have hybrid as an argument");
    if (strcmp(arg[i],"none") == 0)
      error->all(FLERR,"Dihedral style hybrid cannot have none as an argument");

    styles[nstyles] = force->new_dihedral(arg[i],1,dummy);
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

void DihedralHybrid::coeff(int narg, char **arg)
{
  if (!allocated) allocate();

  int ilo,ihi;
  force->bounds(FLERR,arg[0],atom->ndihedraltypes,ilo,ihi);

  // 2nd arg = dihedral sub-style name
  // allow for "none" or "skip" as valid sub-style name

  int m;
  for (m = 0; m < nstyles; m++)
    if (strcmp(arg[1],keywords[m]) == 0) break;

  int none = 0;
  int skip = 0;
  if (m == nstyles) {
    if (strcmp(arg[1],"none") == 0) none = 1;
    else if (strcmp(arg[1],"skip") == 0) none = skip = 1;
    else error->all(FLERR,"Dihedral coeff for hybrid has invalid style");
  }

  // move 1st arg to 2nd arg
  // just copy ptrs, since arg[] points into original input line

  arg[1] = arg[0];

  // invoke sub-style coeff() starting with 1st arg

  if (!none) styles[m]->coeff(narg-1,&arg[1]);

  // set setflag and which type maps to which sub-style
  // if sub-style is skip: auxiliary class2 setting in data file so ignore
  // if sub-style is none and not skip: set hybrid setflag, wipe out map

  for (int i = ilo; i <= ihi; i++) {
    if (skip) continue;
    else if (none) {
      setflag[i] = 1;
      map[i] = -1;
    } else {
      setflag[i] = styles[m]->setflag[i];
      map[i] = m;
    }
  }
}

/* ---------------------------------------------------------------------- */

void DihedralHybrid::init_style()
{
  for (int m = 0; m < nstyles; m++)
    if (styles[m]) styles[m]->init_style();
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void DihedralHybrid::write_restart(FILE *fp)
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

void DihedralHybrid::read_restart(FILE *fp)
{
  int me = comm->me;
  if (me == 0) fread(&nstyles,sizeof(int),1,fp);
  MPI_Bcast(&nstyles,1,MPI_INT,0,world);
  styles = new Dihedral*[nstyles];
  keywords = new char*[nstyles];

  allocate();

  int n,dummy;
  for (int m = 0; m < nstyles; m++) {
    if (me == 0) fread(&n,sizeof(int),1,fp);
    MPI_Bcast(&n,1,MPI_INT,0,world);
    keywords[m] = new char[n];
    if (me == 0) fread(keywords[m],sizeof(char),n,fp);
    MPI_Bcast(keywords[m],n,MPI_CHAR,0,world);
    styles[m] = force->new_dihedral(keywords[m],0,dummy);
  }
}

/* ----------------------------------------------------------------------
   memory usage
------------------------------------------------------------------------- */

double DihedralHybrid::memory_usage()
{
  double bytes = maxeatom * sizeof(double);
  bytes += maxvatom*6 * sizeof(double);
  for (int m = 0; m < nstyles; m++) bytes += maxdihedral[m]*5 * sizeof(int);
  for (int m = 0; m < nstyles; m++)
    if (styles[m]) bytes += styles[m]->memory_usage();
  return bytes;
}
