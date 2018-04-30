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
#include "domain.h"


using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

ComputePairEntropyAtom::ComputePairEntropyAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  pair_entropy(NULL), pair_entropy_avg(NULL)
{
  if (narg < 5 || narg > 8)
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
      else error->all(FLERR,"Illegal compute pentropy/atom argument after avg should be yes or no");
      cutoff2 = force->numeric(FLERR,arg[iarg+2]);
      if (cutoff2 < 0.0) error->all(FLERR,"Illegal compute pentropy/atom command; negative cutoff2");
      cutsq2 = cutoff2*cutoff2;
      iarg += 3;
    } else error->all(FLERR,"Illegal compute pentropy/atom argument after sigma and cutoff should be avg");
  }


  cutsq = cutoff*cutoff;
  nbin = static_cast<int>(cutoff / sigma) + 1;
  nmax = 0;
  maxneigh = 0;
  deltabin = 2;
  deltar = sigma; 
  peratom_flag = 1;
  size_peratom_cols = 0;
}

/* ---------------------------------------------------------------------- */

ComputePairEntropyAtom::~ComputePairEntropyAtom()
{
  memory->destroy(pair_entropy);
  if (avg_flag) memory->destroy(pair_entropy_avg);
}

/* ---------------------------------------------------------------------- */

void ComputePairEntropyAtom::init()
{
  if (force->pair == NULL)
    error->all(FLERR,"Compute centro/atom requires a pair style be defined");

  //double largest_cutsq;
  //largest_cutsq = cutsq;
  //if (cutsq2 > cutsq) largest_cutsq = cutsq2;

  if ( (cutoff+cutoff2) > (force->pair->cutforce  + neighbor->skin) ) 
    {
	//fprintf(screen, "%f %f %f %f \n", cutoff, cutoff2, (cutoff+cutoff2) , force->pair->cutforce  + neighbor->skin );
    	error->all(FLERR,"Compute pentropy/atom cutoff is longer than pairwise cutoff. Increase the neighbor list skin distance.");
    }

  /*
  if (2.0*sqrt(largest_cutsq) > force->pair->cutforce + neighbor->skin &&
      comm->me == 0)
    error->warning(FLERR,"Compute pentropy/atom cutoff may be too large to find "
                   "ghost atom neighbors");
  */

  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"pentropy/atom") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute pentropy/atom");

  // need a full neighbor list with neighbors of the ghost atoms

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->compute = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->occasional = 0;
  neighbor->requests[irequest]->ghost = 1;

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
  double rbin[nbin], rbinsq[nbin];

  invoked_peratom = update->ntimestep;

  // Initialize distance vectors
  for (int i = 0; i < nbin; i++) {
    rbin[i] = i*deltar;
    rbinsq[i] = rbin[i]*rbin[i];
  }

  // grow pair_entropy and pair_entropy_avg array if necessary

  if (atom->nmax > nmax) {
    if (!avg_flag) {
      memory->destroy(pair_entropy);
      nmax = atom->nmax;
      memory->create(pair_entropy,nmax,"pentropy/atom:pair_entropy");
      vector_atom = pair_entropy;
    } else {
      memory->destroy(pair_entropy);
      memory->destroy(pair_entropy_avg);
      nmax = atom->nmax;
      memory->create(pair_entropy,nmax,"pentropy/atom:pair_entropy");
      memory->create(pair_entropy_avg,nmax,"pentropy/atom:pair_entropy_avg");
      vector_atom = pair_entropy_avg;
    }
  }

  // invoke full neighbor list (will copy or build if necessary)

  //neighbor->build_one(list);

  inum = list->inum +  list->gnum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // Compute some constants
  double nlist_cutoff = force->pair->cutforce;
  double sigmasq2=2*sigma*sigma;
  double volume = domain->xprd * domain->yprd * domain->zprd;
  double density = atom->natoms / volume;

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
      double normConstantBase = 4*MY_PI*density; // Normalization of g(r)
      normConstantBase *= sqrt(2.*MY_PI)*sigma; // Normalization of gaussian
      double invNormConstantBase = 1./normConstantBase;

      // loop over list of all neighbors within force cutoff

      // initialize gofr
      double gofr[nbin];
      for(int k=0;k<nbin;++k) gofr[k]=0.;

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
            double invNormKernel=invNormConstantBase/rbinsq[bin];
            double distance = r - rbin[k];
            gofr[k] += invNormKernel*exp(-distance*distance/sigmasq2) ; 
          }
        }
      }

      /*
      if (ii==500) {
        for(int k=0;k<nbin;++k) {
          fprintf(screen,"%f %f \n",rbin[k], gofr[k]);
        }
      }
      */

      // Calculate integrand
      double integrand[nbin];
      for(int k=0;k<nbin;++k){
        if (gofr[k]<1.e-10) {
          integrand[k] = rbinsq[k];
        } else {
          integrand[k] = (gofr[k]*log(gofr[k])-gofr[k]+1)*rbinsq[k];
        }
      }

      // Integrate with trapezoid rule
      double value = 0.;
      for(int k=1;k<nbin-1;++k){
        value += integrand[k];
      }
      value += 0.5*integrand[0];
      value += 0.5*integrand[nbin-1];
      value *= deltar;

      pair_entropy[i] = -2*MY_PI*density*value;

    }
  }


  if (avg_flag) {
    for (ii = 0; ii < inum; ii++) {
      i = ilist[ii];
      if (mask[i] & groupbit) {
        xtmp = x[i][0];
        ytmp = x[i][1];
        ztmp = x[i][2];
        jlist = firstneigh[i];
        jnum = numneigh[i];
 
        pair_entropy_avg[i] = pair_entropy[i];
        double counter = 1;
        // loop over list of all neighbors within force cutoff
 
        for (jj = 0; jj < jnum; jj++) {
          j = jlist[jj];
          j &= NEIGHMASK;
 
          delx = xtmp - x[j][0];
          dely = ytmp - x[j][1];
          delz = ztmp - x[j][2];
          rsq = delx*delx + dely*dely + delz*delz;
          if (rsq < cutsq2) {
            pair_entropy_avg[i] += pair_entropy[j];
            counter += 1;
          }
        }
        pair_entropy_avg[i] /= counter;
      }
    }
  }

}


/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputePairEntropyAtom::memory_usage()
{
  double bytes;
  if (!avg_flag) {
    bytes = nmax * sizeof(double);
  } else {
    bytes = 2 * nmax * sizeof(double);
  }
  return bytes;
}
