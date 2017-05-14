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

/* ----------------------------------------------------------------------
   Contributing author: Axel Kohlmeyer and Richard Berger (Temple U)
------------------------------------------------------------------------- */

#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pair_python.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "python.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairPython::PairPython(LAMMPS *lmp) : Pair(lmp) {
  respa_enable = 0;
  single_enable = 0;
  writedata = 0;
  restartinfo = 0;
  one_coeff = 1;
  reinitflag = 0;

  python->init();
}

/* ---------------------------------------------------------------------- */

PairPython::~PairPython()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
  }
}

/* ---------------------------------------------------------------------- */

void PairPython::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,factor_lj;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        const double r = sqrt(rsq);
        printf("compute f at r=%g for types %d,%d with factor %g\n",
               r,itype,jtype,factor_lj);
        fpair = 0.0;
        fpair *= factor_lj;

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        if (eflag) {
          printf("compute e at r=%g for types %d,%d with factor %g\n",
                 r,itype,jtype,factor_lj);
          evdwl = 0.0;
          evdwl *= factor_lj;
        }

        if (evflag) ev_tally(i,j,nlocal,newton_pair,
                             evdwl,0.0,fpair,delx,dely,delz);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairPython::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairPython::settings(int narg, char **arg)
{
  if (narg != 1)
    error->all(FLERR,"Illegal pair_style command");

  cut_global = force->numeric(FLERR,arg[0]);
}

/* ----------------------------------------------------------------------
   set coeffs for all type pairs
------------------------------------------------------------------------- */

void PairPython::coeff(int narg, char **arg)
{
  const int ntypes = atom->ntypes;

  if (narg != 3+ntypes)
    error->all(FLERR,"Incorrect args for pair coefficients");

  if (!allocated) allocate();

  // make sure I,J args are * *

  if (strcmp(arg[0],"*") != 0 || strcmp(arg[1],"*") != 0)
    error->all(FLERR,"Incorrect args for pair coefficients");

  // check if potential python file exists

  FILE *fp = fopen(arg[2],"r");
  if (fp == NULL)
    error->all(FLERR,"Cannot open python pair potential class file");

  PyGILState_STATE gstate = PyGILState_Ensure();
  PyObject *pModule = PyImport_AddModule("__main__");
  if (!pModule) error->all(FLERR,"Could not initialize embedded Python");

  int err = PyRun_SimpleFile(fp,arg[2]);
  if (err) error->all(FLERR,"Loading python pair style class failure");
  fclose(fp);

  PyObject *py_pair_instance =
    PyObject_GetAttrString(pModule,"lammps_pair_style");

  if (!py_pair_instance) {
    PyGILState_Release(gstate);
    error->all(FLERR,"Could not find 'lammps_pair_style instance'");
  }


  for (int i = 1; i <= ntypes ; i++) {
    for (int j = i; j <= ntypes ; j++) {
      if (strcmp(arg[2+i],"NULL") != 0) {
        setflag[i][j] = 1;
        cutsq[i][j] = cut_global*cut_global;
      }
    }
  }
  PyGILState_Release(gstate);
}

/* ---------------------------------------------------------------------- */

double PairPython::init_one(int, int)
{
  return cut_global;
}

