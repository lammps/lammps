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

#include "random_park.h"
#include <cmath>
#include "error.h"
#include "math_const.h"

using namespace LAMMPS_NS;

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836

/* ---------------------------------------------------------------------- */

RanPark::RanPark(LAMMPS *lmp, int seed_init) : Pointers(lmp)
{
  if (seed_init <= 0)
    error->one(FLERR,"Invalid seed for Park random # generator");
  seed = seed_init;
  save = 0;
}

/* ----------------------------------------------------------------------
   uniform RN
------------------------------------------------------------------------- */

double RanPark::uniform()
{
  int k = seed/IQ;
  seed = IA*(seed-k*IQ) - IR*k;
  if (seed < 0) seed += IM;
  double ans = AM*seed;
  return ans;
}

/* ----------------------------------------------------------------------
   gaussian RN
------------------------------------------------------------------------- */

double RanPark::gaussian()
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

/* ---------------------------------------------------------------------- */

void RanPark::reset(int seed_init)
{
  if (seed_init <= 0)
    error->all(FLERR,"Invalid seed for Park random # generator");
  seed = seed_init;
  save = 0;
}

/* ----------------------------------------------------------------------
   reset the seed based on atom coords and ibase = caller seed
   use hash function, treating user seed and coords as sequence of input ints
   this is Jenkins One-at-a-time hash, see Wikipedia entry on hash tables
------------------------------------------------------------------------- */

void RanPark::reset(int ibase, double *coord)
{
  int i;

  char *str = (char *) &ibase;
  int n = sizeof(int);

  unsigned int hash = 0;
  for (i = 0; i < n; i++) {
    hash += str[i];
    hash += (hash << 10);
    hash ^= (hash >> 6);
  }

  str = (char *) coord;
  n = 3 * sizeof(double);
  for (i = 0; i < n; i++) {
    hash += str[i];
    hash += (hash << 10);
    hash ^= (hash >> 6);
  }

  hash += (hash << 3);
  hash ^= (hash >> 11);
  hash += (hash << 15);

  // keep 31 bits of unsigned int as new seed
  // do not allow seed = 0, since will cause hang in gaussian()

  seed = hash & 0x7ffffff;
  if (!seed) seed = 1;

  // warm up the RNG

  for (i = 0; i < 5; i++) uniform();
  save = 0;
}

/* ---------------------------------------------------------------------- */

int RanPark::state()
{
  return seed;
}
/* ----------------------------------------------------------------------
   Gaussian RN
------------------------------------------------------------------------- */

double RanPark::gaussian(double mu, double sigma)
{
  double v1;
  v1 = mu+sigma*gaussian();
  return v1;
}

/* ----------------------------------------------------------------------
   Rayleigh RN
------------------------------------------------------------------------- */

double RanPark::rayleigh(double sigma)
{
  double first,v1;

  if (sigma <= 0) error->all(FLERR,"Invalid Rayleigh parameter");

  v1 = uniform();
  first = sigma*sqrt(-2.0*log(v1));
  return first;
}

/* ----------------------------------------------------------------------
   Bessel exponential RN
------------------------------------------------------------------------- */

double RanPark::besselexp(double theta, double alpha, double cp)
{
  double first,v1,v2;

  if (theta < 0.0 || alpha < 0.0 || alpha < 1.0)
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
