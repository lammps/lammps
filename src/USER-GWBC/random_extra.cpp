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

// Park/Miller RNG

#include "random_extra.h"
#include <cmath>
#include "error.h"

using namespace LAMMPS_NS;

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836

/* ---------------------------------------------------------------------- */

RanExtra::RanExtra(LAMMPS *lmp, int seed_init) : Pointers(lmp)
{
  if (seed_init <= 0)
    error->one(FLERR,"Invalid seed for Park random # generator (random extra)");
  seed = seed_init;
  save = 0;
}

/* ----------------------------------------------------------------------
   uniform RN
------------------------------------------------------------------------- */

double RanExtra::uniform()
{
  int k = seed/IQ;
  seed = IA*(seed-k*IQ) - IR*k;
  if (seed < 0) seed += IM;
  double ans = AM*seed;
  return ans;
}

/* ----------------------------------------------------------------------
   Rayleigh RN
------------------------------------------------------------------------- */

double RanExtra::rayleigh(double sigma)
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
/* ---------------------------------------------------------------------- */


/* ----------------------------------------------------------------------
   Gaussian RN
------------------------------------------------------------------------- */

double RanExtra::gaussian(double mu, double sigma)
{
  double first,v1,v2,rsq,fac;

    do {
      v1 = 2.0*uniform()-1.0;
      v2 = 2.0*uniform()-1.0;
      rsq = v1*v1 + v2*v2;
    } while ((rsq >= 1.0) || (rsq == 0.0));
    fac = sqrt(-2.0*log(rsq)/rsq);

    first = sigma*v1*fac+mu;
    return first;
}
/* ---------------------------------------------------------------------- */


/* ----------------------------------------------------------------------
   Bessel exponential RN
------------------------------------------------------------------------- */

double RanExtra::besselexp(double theta, double alpha, double cp)
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
/* ---------------------------------------------------------------------- */


/* ---------------------------------------------------------------------- */

void RanExtra::reset(int seed_init)
{
  if (seed_init <= 0)
    error->all(FLERR,"Invalid seed for Park random # generator (random extra)");
  seed = seed_init;
  save = 0;
}

/* ---------------------------------------------------------------------- */

int RanExtra::state()
{
  return seed;
}
