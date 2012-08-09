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

#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "compute_coord_atom.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "force.h"
#include "pair.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeCoordAtom::ComputeCoordAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal compute coord/atom command");

  double cutoff = atof(arg[3]);
  cutsq = cutoff*cutoff;

  ncol = narg-4 + 1;
  int ntypes = atom->ntypes;
  typelo = new int[ncol];
  typehi = new int[ncol];

  if (narg == 4) {
    ncol = 1;
    typelo[0] = 1;
    typehi[0] = ntypes;
  } else {
    ncol = 0;
    int iarg = 4;
    while (iarg < narg) {
      force->bounds(arg[iarg],ntypes,typelo[ncol],typehi[ncol]);
      if (typelo[ncol] > typehi[ncol])
        error->all(FLERR,"Illegal compute coord/atom command");
      ncol++;
      iarg++;
    }
  }

  peratom_flag = 1;
  if (ncol == 1) size_peratom_cols = 0;
  else size_peratom_cols = ncol;

  nmax = 0;
  cvec = NULL;
  carray = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeCoordAtom::~ComputeCoordAtom()
{
  delete [] typelo;
  delete [] typehi;
  memory->destroy(cvec);
  memory->destroy(carray);
}

/* ---------------------------------------------------------------------- */

void ComputeCoordAtom::init()
{
  if (force->pair == NULL)
    error->all(FLERR,"Compute coord/atom requires a pair style be defined");
  if (sqrt(cutsq) > force->pair->cutforce)
    error->all(FLERR,
               "Compute coord/atom cutoff is longer than pairwise cutoff");

  // need an occasional full neighbor list

  int irequest = neighbor->request((void *) this);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->compute = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->occasional = 1;

  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"coord/atom") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute coord/atom");
}

/* ---------------------------------------------------------------------- */

void ComputeCoordAtom::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeCoordAtom::compute_peratom()
{
  int i,j,m,ii,jj,inum,jnum,jtype,n;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double *count;

  invoked_peratom = update->ntimestep;

  // grow coordination array if necessary

  if (atom->nlocal > nmax) {
    if (ncol == 1) {
      memory->destroy(cvec);
      nmax = atom->nmax;
      memory->create(cvec,nmax,"coord/atom:cvec");
      vector_atom = cvec;
    } else {
      memory->destroy(carray);
      nmax = atom->nmax;
      memory->create(carray,nmax,ncol,"coord/atom:carray");
      array_atom = carray;
    }
  }

  // invoke full neighbor list (will copy or build if necessary)

  neighbor->build_one(list->index);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // compute coordination number(s) for each atom in group
  // use full neighbor list to count atoms less than cutoff

  double **x = atom->x;
  int *type = atom->type;
  int *mask = atom->mask;

  if (ncol == 1) {
    for (ii = 0; ii < inum; ii++) {
      i = ilist[ii];
      if (mask[i] & groupbit) {
        xtmp = x[i][0];
        ytmp = x[i][1];
        ztmp = x[i][2];
        jlist = firstneigh[i];
        jnum = numneigh[i];

        n = 0;
        for (jj = 0; jj < jnum; jj++) {
          j = jlist[jj];
          j &= NEIGHMASK;
          
          jtype = type[j];
          delx = xtmp - x[j][0];
          dely = ytmp - x[j][1];
          delz = ztmp - x[j][2];
          rsq = delx*delx + dely*dely + delz*delz;
          if (rsq < cutsq && jtype >= typelo[0] && jtype <= typehi[0]) n++;
        }
        
        cvec[i] = n;
      } else cvec[i] = 0.0;
    }

  } else {
    for (ii = 0; ii < inum; ii++) {
      i = ilist[ii];
      count = carray[i];
      for (m = 0; m < ncol; m++) count[m] = 0.0;

      if (mask[i] & groupbit) {
        xtmp = x[i][0];
        ytmp = x[i][1];
        ztmp = x[i][2];
        jlist = firstneigh[i];
        jnum = numneigh[i];


        for (jj = 0; jj < jnum; jj++) {
          j = jlist[jj];
          j &= NEIGHMASK;
          
          jtype = type[j];
          delx = xtmp - x[j][0];
          dely = ytmp - x[j][1];
          delz = ztmp - x[j][2];
          rsq = delx*delx + dely*dely + delz*delz;
          if (rsq < cutsq) {
            for (m = 0; m < ncol; m++)
              if (jtype >= typelo[m] && jtype <= typehi[m])
                count[m] += 1.0;
          }
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeCoordAtom::memory_usage()
{
  double bytes = ncol*nmax * sizeof(double);
  return bytes;
}
