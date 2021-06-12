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

// Marsaglia random number generator
// see RANMAR in F James, Comp Phys Comm, 60, 329 (1990)

#include "random_mars.h"
#include <cmath>
#include <cstring>
#include "error.h"
#include "math_const.h"

using namespace LAMMPS_NS;

enum{ADD,SUBTRACT};

/* ---------------------------------------------------------------------- */

RanMars::RanMars(LAMMPS *lmp, int seed) : Pointers(lmp),
  u(nullptr)
{
  int ij,kl,i,j,k,l,ii,jj,m;
  double s,t;

  if (seed <= 0 || seed > 900000000)
    error->one(FLERR,"Invalid seed for Marsaglia random # generator");

  save = 0;
  u = new double[97+1];
  memset(u,0,98*sizeof(double));

  ij = (seed-1)/30082;
  kl = (seed-1) - 30082*ij;
  i = (ij/177) % 177 + 2;
  j = ij %177 + 2;
  k = (kl/169) % 178 + 1;
  l = kl % 169;
  for (ii = 1; ii <= 97; ii++) {
    s = 0.0;
    t = 0.5;
    for (jj = 1; jj <= 24; jj++) {
      m = ((i*j) % 179)*k % 179;
      i = j;
      j = k;
      k = m;
      l = (53*l+1) % 169;
      if ((l*m) % 64 >= 32) s = s + t;
      t = 0.5*t;
    }
    u[ii] = s;
  }
  c = 362436.0 / 16777216.0;
  cd = 7654321.0 / 16777216.0;
  cm = 16777213.0 / 16777216.0;
  i97 = 97;
  j97 = 33;
  uniform();
}

/* ---------------------------------------------------------------------- */

RanMars::~RanMars()
{
  delete [] u;
}

/* ----------------------------------------------------------------------
   uniform RN
------------------------------------------------------------------------- */

double RanMars::uniform()
{
  double uni = u[i97] - u[j97];
  if (uni < 0.0) uni += 1.0;
  u[i97] = uni;
  i97--;
  if (i97 == 0) i97 = 97;
  j97--;
  if (j97 == 0) j97 = 97;
  c -= cd;
  if (c < 0.0) c += cm;
  uni -= c;
  if (uni < 0.0) uni += 1.0;
  return uni;
}


/* ----------------------------------------------------------------------
   gaussian RN
------------------------------------------------------------------------- */

double RanMars::gaussian()
{
  double first,v1,v2,rsq,fac;

  if (!save) {
    do {
      v1 = 2.0*uniform()-1.0;
      v2 = 2.0*uniform()-1.0;
      rsq = v1*v1 + v2*v2;
    } while ((rsq >= 1.0) || (rsq == 0.0));
    fac = sqrt(-2.0*log(rsq)/rsq);
    second = v1*fac;
    first = v2*fac;
    save = 1;
  } else {
    first = second;
    save = 0;
  }
  return first;
}

/* ----------------------------------------------------------------------
   Gaussian RN
------------------------------------------------------------------------- */

double RanMars::gaussian(double mu, double sigma)
{
  double v1;
  v1 = mu+sigma*gaussian();
  return v1;
}

/* ----------------------------------------------------------------------
   Rayleigh RN
------------------------------------------------------------------------- */

double RanMars::rayleigh(double sigma)
{
  double v1;

  if (sigma <= 0.0) error->all(FLERR,"Invalid Rayleigh parameter");

  v1 = uniform();
  // avoid a floating point exception due to log(0.0)
  // and just return a very big number
  if (v1 == 0.0) return 1.0e300;

  return sigma*sqrt(-2.0*log(v1));
}

/* ----------------------------------------------------------------------
   Bessel exponential RN
------------------------------------------------------------------------- */

double RanMars::besselexp(double theta, double alpha, double cp)
{
  double first,v1,v2;

  if (theta < 0.0 || alpha < 0.0 || alpha > 1.0)
    error->all(FLERR,"Invalid Bessel exponential distribution parameters");

  v1 = uniform();
  v2 = uniform();

  if (cp < 0.0)
    first = sqrt((1.0-alpha)*cp*cp - 2.0*alpha*theta*log(v1) +
                 2.0*sqrt(-2.0*theta*(1.0-alpha)*alpha*log(v1)) *
                 cos(2.0*MathConst::MY_PI*v2)*cp);
  else
    first = -sqrt((1.0-alpha)*cp*cp - 2.0*alpha*theta*log(v1) -
                  2.0*sqrt(-2.0*theta*(1.0-alpha)*alpha*log(v1)) *
                  cos(2.0*MathConst::MY_PI*v2)*cp);

  return first;
}

/* ----------------------------------------------------------------------
   select random subset of size Ntarget out of Ntotal items
   Ntotal = sum of Nmine across all procs
   mark,next = vectors of length Nmine
   return mark = 0 for unselected item, 1 for selected item
   next = work vector used to store linked lists for 2 active sets of items
   IMPORTANT: this method must be called simultaneously by all procs
------------------------------------------------------------------------- */

void RanMars::select_subset(bigint ntarget, int nmine, int *mark, int *next)
{
  int mode,index,oldindex,newvalue,nflip,which,niter;
  int active[2],first[2];
  int newactive[2],newfirst[2],newlast[2];
  bigint nmark,nflipall;
  bigint activeall[2],bsum[3],bsumall[3];
  double thresh;

  active[0] = nmine;
  active[1] = 0;
  first[0] = 0;
  first[1] = -1;

  bigint bnmine = nmine;
  bigint bnall;
  MPI_Allreduce(&bnmine,&bnall,1,MPI_LMP_BIGINT,MPI_SUM,world);
  activeall[0] = bnall;

  for (int i = 0; i < nmine; i++) mark[i] = 0;
  for (int i = 0; i < nmine; i++) next[i] = i+1;
  if (nmine > 0) next[nmine-1] = -1;

  nmark = 0;
  niter = 0;

  while (nmark != ntarget) {

    // choose to ADD or SUBTRACT from current nmark
    // thresh = desired flips / size of active set

    if (ntarget-nmark > 0) {
      mode = ADD;
      thresh = 1.0 * (ntarget-nmark) / activeall[mode];
    } else {
      mode = SUBTRACT;
      thresh = 1.0 * (nmark-ntarget) / activeall[mode];
    }

    // bound thresh for RNG accuracy

    thresh = MAX(thresh,0.01);
    thresh = MIN(thresh,0.99);

    // new empty active sets for next iteration

    newactive[0] = newactive[1] = 0;
    newfirst[0] = newfirst[1] = -1;
    newlast[0] = newlast[1] = -1;

    // index = first value in ADD or SUBTRACT set

    if (mode == ADD) newvalue = 1;
    else if (mode == SUBTRACT) newvalue = 0;
    index = first[mode];

    // flip marks from 0 -> 1 (ADD) or 1 -> 0 (SUBTRACT)
    // loop over active set via next vector = linked list
    // flip each value based on RN < thresh

    nflip = 0;
    while ((nmine > 0) && (index >= 0)) {
      if (uniform() < thresh) {
        mark[index] = newvalue;
        nflip++;
      }
      oldindex = index;
      index = next[index];

      // oldindex can now be appended to a new active set
      // which = which of two new active sets to append to

      which = mark[oldindex];
      newactive[which]++;
      if (newfirst[which] < 0) newfirst[which] = oldindex;
      if (newlast[which] >= 0) next[newlast[which]] = oldindex;
      newlast[which] = oldindex;
      next[oldindex] = -1;

      // set active sets for next iteration to the new ones
      // next vector is already updated

      active[0] = newactive[0];
      active[1] = newactive[1];
      first[0] = newfirst[0];
      first[1] = newfirst[1];
    }

    // update nmark and activeall

    bsum[0] = nflip;
    bsum[1] = active[0];
    bsum[2] = active[1];
    MPI_Allreduce(&bsum,&bsumall,3,MPI_LMP_BIGINT,MPI_SUM,world);
    nflipall = bsumall[0];
    activeall[0] = bsumall[1];
    activeall[1] = bsumall[2];

    if (mode == ADD) nmark += nflipall;
    else if (mode == SUBTRACT) nmark -= nflipall;

    niter++;

    // DEBUG output of stats

    //if (comm->me == 0) printf("%d %ld %ld %g %ld\n",
    //                          niter,nmark,nactiveall,thresh,nflipall);
  }
}

/* ----------------------------------------------------------------------
   store state in buffer
------------------------------------------------------------------------- */

void RanMars::get_state(double *state)
{
  for (int i=0; i < 98; ++i) state[i] = u[i];
  state[98] = i97;
  state[99] = j97;
  state[100]= c;
  state[101]= cd;
  state[102]= cm;
}

/* ----------------------------------------------------------------------
   restore state from buffer
------------------------------------------------------------------------- */

void RanMars::set_state(double *state)
{
  for (int i=0; i < 98; ++i) u[i] = state[i];
  i97 = state[98];
  j97 = state[99];
  c   = state[100];
  cd  = state[101];
  cm  = state[102];
}
