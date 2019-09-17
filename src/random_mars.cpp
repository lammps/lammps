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

// Marsaglia random number generator
// see RANMAR in F James, Comp Phys Comm, 60, 329 (1990)

#include "random_mars.h"
#include <cmath>
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

RanMars::RanMars(LAMMPS *lmp, int seed) : Pointers(lmp),
  u(NULL)
{
  int ij,kl,i,j,k,l,ii,jj,m;
  double s,t;

  if (seed <= 0 || seed > 900000000)
    error->one(FLERR,"Invalid seed for Marsaglia random # generator");

  save = 0;
  u = new double[97+1];

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
/* ---------------------------------------------------------------------- */


/* ----------------------------------------------------------------------
   Rayleigh RN
------------------------------------------------------------------------- */

double RanMars::rayleigh(double sigma)
{
  double first,v1;
 
   if (sigma <= 0)
    error->all(FLERR,"Invalid Rayleigh parameter");
   else {
    v1 = uniform();
    first = sigma*sqrt(-2.0*log(v1));
    return first;
   }
}


/* ----------------------------------------------------------------------
   Bessel exponential RN
------------------------------------------------------------------------- */

double RanMars::besselexp(double theta, double alpha, double cp)
{
  double first,v1,v2;


   if ((theta < 0) || (alpha < 0) || (alpha >1))
    error->all(FLERR,"Invalid Bessel exponential distribution parameters");
   else {
    v1 = uniform();
    v2 = uniform();
    if (cp < 0) 
    first = sqrt((1-alpha)*cp*cp-2*alpha*theta*log(v1)+2*sqrt(-2*theta*(1-alpha)*alpha*log(v1))*cos(2*M_PI*v2)*cp);
    else {
    first = - sqrt((1-alpha)*cp*cp-2*alpha*theta*log(v1)-2*sqrt(-2*theta*(1-alpha)*alpha*log(v1))*cos(2*M_PI*v2)*cp) ;}
        return first;
   }
}