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

#include "compute_entropy_atom.h"
#include <cmath>
#include <cstring>
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "force.h"
#include "pair.h"
#include "comm.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "domain.h"

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

ComputeEntropyAtom::
ComputeEntropyAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  pair_entropy(NULL), pair_entropy_avg(NULL)
{
  if (narg < 5 || narg > 10)
    error->all(FLERR,"Illegal compute entropy/atom command; wrong number"
               " of arguments");

  // Arguments are: sigma cutoff avg yes/no cutoff2 local yes/no
  //   sigma is the gaussian width
  //   cutoff is the cutoff for the calculation of g(r)
  //   avg is optional and allows averaging the pair entropy over neighbors
  //   the next argument should be yes or no
  //   cutoff2 is the cutoff for the averaging
  //   local is optional and allows using the local density to normalize
  //     the g(r)

  sigma = force->numeric(FLERR,arg[3]);
  if (sigma <= 0.0) error->all(FLERR,"Illegal compute entropy/atom"
                              " command; sigma must be positive");
  cutoff = force->numeric(FLERR,arg[4]);
  if (cutoff <= 0.0) error->all(FLERR,"Illegal compute entropy/atom"
                               " command; cutoff must be positive");

  cutoff2 = 0.;
  avg_flag = 0;
  local_flag = 0;

  // optional keywords
  int iarg = 5;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"avg") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal compute entropy/atom;"
                   " missing arguments after avg");
      if (strcmp(arg[iarg+1],"yes") == 0) avg_flag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) avg_flag = 0;
      else error->all(FLERR,"Illegal compute entropy/atom;"
                      " argument after avg should be yes or no");
      cutoff2 = force->numeric(FLERR,arg[iarg+2]);
      if (cutoff2 < 0.0) error->all(FLERR,"Illegal compute entropy/atom"
                                    " command; negative cutoff2");
      cutsq2 = cutoff2*cutoff2;
      iarg += 3;
    } else if (strcmp(arg[iarg],"local") == 0) {
      if (strcmp(arg[iarg+1],"yes") == 0) local_flag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) local_flag = 0;
      else error->all(FLERR,"Illegal compute entropy/atom;"
                      " argument after local should be yes or no");
      iarg += 2;
    } else error->all(FLERR,"Illegal compute entropy/atom; argument after"
                      " sigma and cutoff should be avg or local");
  }


  cutsq = cutoff*cutoff;
  nbin = static_cast<int>(cutoff / sigma) + 1;
  nmax = 0;
  maxneigh = 0;
  // Number of bins above and below the central one that will be
  // considered as affected by the gaussian kernel
  // 3 seems a good compromise between speed and good mollification
  deltabin = 3;
  deltar = sigma;
  peratom_flag = 1;
  size_peratom_cols = 0;
}

/* ---------------------------------------------------------------------- */

ComputeEntropyAtom::~ComputeEntropyAtom()
{
  memory->destroy(pair_entropy);
  if (avg_flag) memory->destroy(pair_entropy_avg);
}

/* ---------------------------------------------------------------------- */

void ComputeEntropyAtom::init()
{
  if (force->pair == NULL)
    error->all(FLERR,"Compute entropy/atom requires a pair style be"
               " defined");

  if ( (cutoff+cutoff2) > (force->pair->cutforce  + neighbor->skin) )
    {
        error->all(FLERR,"Compute entropy/atom cutoff is longer than the"
                   " pairwise cutoff. Increase the neighbor list skin"
                   " distance.");
    }

  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"entropy/atom") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute entropy/atom");

  // Request neighbor list
  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->compute = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->occasional = 1;
  if (avg_flag) {
    // need a full neighbor list with neighbors of the ghost atoms
    neighbor->requests[irequest]->ghost = 1;
  } else {
    // need a full neighbor list
    neighbor->requests[irequest]->ghost = 0;
  }

}

/* ---------------------------------------------------------------------- */

void ComputeEntropyAtom::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeEntropyAtom::compute_peratom()
{
  int i,j,ii,jj,inum,jnum;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double *rbin = new double[nbin];
  double *rbinsq = new double[nbin];

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
      memory->create(pair_entropy,nmax,"entropy/atom:pair_entropy");
      vector_atom = pair_entropy;
    } else {
      memory->destroy(pair_entropy);
      memory->destroy(pair_entropy_avg);
      nmax = atom->nmax;
      memory->create(pair_entropy,nmax,"entropy/atom:pair_entropy");
      memory->create(pair_entropy_avg,nmax,
                     "entropy/atom:pair_entropy_avg");
      vector_atom = pair_entropy_avg;
    }
  }

  // invoke full neighbor list (will copy or build if necessary)

  neighbor->build_one(list);

  inum = list->inum + list->gnum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // Compute some constants
  double sigmasq2=2*sigma*sigma;
  double volume = domain->xprd * domain->yprd * domain->zprd;
  double density = atom->natoms / volume;

  // compute pair entropy for each atom in group
  // use full neighbor list

  double **x = atom->x;
  int *mask = atom->mask;
  double *gofr = new double[nbin];
  double *integrand = new double[nbin];

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (mask[i] & groupbit) {
      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];
      jlist = firstneigh[i];
      jnum = numneigh[i];

      // If local density is used, calculate it
      if (local_flag) {
        double neigh_cutoff = force->pair->cutforce  + neighbor->skin;
        double volume =
               (4./3.)*MY_PI*neigh_cutoff*neigh_cutoff*neigh_cutoff;
        density = jnum / volume;
      }

      // calculate kernel normalization
      // Normalization of g(r)
      double normConstantBase = 4*MY_PI*density;
      // Normalization of gaussian
      normConstantBase *= sqrt(2.*MY_PI)*sigma;
      double invNormConstantBase = 1./normConstantBase;

      // loop over list of all neighbors within force cutoff

      // initialize gofr
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
          int minbin, maxbin;
          minbin=bin - deltabin;
          if (minbin < 0) minbin=0;
          if (minbin > (nbin-1)) minbin=nbin-1;
          maxbin=bin +  deltabin;
          if (maxbin > (nbin-1)) maxbin=nbin-1;
          for(int k=minbin;k<maxbin+1;k++) {
            double invNormKernel=invNormConstantBase/rbinsq[k];
            double distance = r - rbin[k];
            gofr[k] += invNormKernel*exp(-distance*distance/sigmasq2);
          }
        }
      }

      // Calculate integrand
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
  delete [] gofr;
  delete [] integrand;

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
  delete [] rbin;
  delete [] rbinsq;
}


/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeEntropyAtom::memory_usage()
{
  double bytes;
  if (!avg_flag) {
    bytes = nmax * sizeof(double);
  } else {
    bytes = 2 * nmax * sizeof(double);
  }
  return bytes;
}
