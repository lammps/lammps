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
   Contributing author: Pablo Piaggi (EPFL Lausanne)
------------------------------------------------------------------------- */

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "compute_pair_entropy_atom.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "force.h"
#include "pair.h"
#include "comm.h"
#include "math_extra.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

ComputePairEntropyAtom::ComputePairEntropyAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  gofr(NULL), pair_entropy(NULL)
{
  if (narg < 5 || narg > 7)
    error->all(FLERR,"Illegal compute pentropy/atom command");

  sigma = force->numeric(FLERR,arg[3]);
  cutoff = force->numeric(FLERR,arg[4]);
  if (cutoff < 0.0) error->all(FLERR,"Illegal compute pentropy/atom command; negative cutoff");

  avg_flag = 0;

  // optional keywords

  int iarg = 5;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"avg") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal compute pentropy/atom missing arguments after avg");
      if (strcmp(arg[iarg+1],"yes") == 0) avg_flag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) avg_flag = 0;
      else error->all(FLERR,"Illegal compute centro/atom argument after avg should be yes or no");
      cutoff2 = force->numeric(FLERR,arg[iarg+2]);
      if (cutoff2 < 0.0) error->all(FLERR,"Illegal compute pentropy/atom command; negative cutoff2");
      cutsq2 = cutoff2*cutoff2;
      iarg += 2;
    } else error->all(FLERR,"Illegal compute centro/atom argument after sigma and cutoff should be avg");
  }

  peratom_flag = 1;

  cutsq = cutoff*cutoff;
  nbin = static_cast<int>(cutoff / sigma) + 1;
  nmax = 0;
  maxneigh = 0;
  deltabin = 2;
  deltar = sigma; 
  rbin = new double(nbin);
  rbinsq = new double(nbin);
  for (int i = 0; i < nbin; i++) {
    rbin[i] = i*deltar;
    rbinsq[i] = rbin[i]*rbin[i];
  }
}

/* ---------------------------------------------------------------------- */

ComputePairEntropyAtom::~ComputePairEntropyAtom()
{
  memory->destroy(pair_entropy);
  memory->destroy(gofr);
}

/* ---------------------------------------------------------------------- */

void ComputePairEntropyAtom::init()
{
  if (force->pair == NULL)
    error->all(FLERR,"Compute centro/atom requires a pair style be defined");

  double largest_cutsq;
  largest_cutsq = cutsq;
  if (cutsq2 > cutsq) largest_cutsq = cutsq2;

  if (sqrt(largest_cutsq) > force->pair->cutforce)
    error->all(FLERR,"Compute pentropy/atom cutoff is longer than pairwise cutoff");

  if (2.0*sqrt(largest_cutsq) > force->pair->cutforce + neighbor->skin &&
      comm->me == 0)
    error->warning(FLERR,"Compute pentropy/atom cutoff may be too large to find "
                   "ghost atom neighbors");

  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"pentropy/atom") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute pentropy/atom");

  // need an occasional full neighbor list

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->compute = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->occasional = 1;
}

/* ---------------------------------------------------------------------- */

void ComputePairEntropyAtom::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputePairEntropyAtom::compute_peratom()
{
  int i,j,k,ii,jj,kk,n,inum,jnum;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq,value;
  int *ilist,*jlist,*numneigh,**firstneigh;

  invoked_peratom = update->ntimestep;

  // grow pair_entropy array if necessary

  if (atom->nmax > nmax) {
    memory->destroy(pair_entropy);
    nmax = atom->nmax;
    memory->create(pair_entropy,nmax,"pentropy/atom:pair_entropy");
    vector_atom = pair_entropy;
  }

  // invoke full neighbor list (will copy or build if necessary)

  neighbor->build_one(list);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // compute pair entropy for each atom in group
  // use full neighbor list

  double **x = atom->x;
  int *mask = atom->mask;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (mask[i] & groupbit) {
      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];
      jlist = firstneigh[i];
      jnum = numneigh[i];

      
      // calculate kernel normalization
      double nlist_cutoff = force->pair->cutforce;
      double density = jnum/((4./3.)*MY_PI*cutoff*nlist_cutoff*nlist_cutoff*nlist_cutoff);
      double normConstantBase = 2*MY_PI*density; // Normalization of g(r)
      normConstantBase *= sqrt(2.*MY_PI)*sigma; // Normalization of gaussian
      double invNormConstantBase = 1./normConstantBase;
      double sigmasq2=2*sigma*sigma;

      // loop over list of all neighbors within force cutoff

      // re-initialize gofr
      delete [] gofr;
      gofr = new double(nbin);

      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        j &= NEIGHMASK;

        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
        rsq = delx*delx + dely*dely + delz*delz;
        if (rsq < cutsq) {
          // contribute to gofr
          double r=sqrt(rsq);
          int bin=floor(r/deltar);
          int minbin, maxbin; // These cannot be unsigned
          // Only consider contributions to g(r) of atoms less than n*sigma bins apart from the actual distance
          minbin=bin - deltabin;
          if (minbin < 0) minbin=0;
          if (minbin > (nbin-1)) minbin=nbin-1;
          maxbin=bin +  deltabin;
          if (maxbin > (nbin-1)) maxbin=nbin-1;
          for(int k=minbin;k<maxbin+1;k++) {
            double invNormKernel=invNormConstantBase/rbinsq[k];
            double distance = r - rbin[k];
            gofr[k] += invNormKernel*exp(-distance*distance/sigmasq2) ; 
          }
        }
      }

      //value = 0.0;
      //for (j = 0; j < nhalf; j++) value += pairs[j];
      pair_entropy[i] = gofr[0];

    }
  }

  //delete [] pairs;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (mask[i] & groupbit)
      array_atom[i][0] = pair_entropy[i];
  }
}


/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputePairEntropyAtom::memory_usage()
{
  double bytes = nmax * sizeof(double);
  return bytes;
}
