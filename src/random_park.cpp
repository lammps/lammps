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

#include "math.h"
#include "random_park.h"
#include "domain.h"
#include "error.h"

using namespace LAMMPS_NS;

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836

#define IA1 1366
#define IC1 150889
#define IM1 714025
#define IA2 8121
#define IC2 28411
#define IM2 134456
#define IA3 7141
#define IC3 54773
#define IM3 259200

/* ---------------------------------------------------------------------- */

RanPark::RanPark(LAMMPS *lmp, int seed_init) : Pointers(lmp)
{
  if (seed_init <= 0) error->all("Invalid seed for Park random # generator");
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
    int again = 1;
    while (again) {
      v1 = 2.0*uniform()-1.0;
      v2 = 2.0*uniform()-1.0;
      rsq = v1*v1 + v2*v2;
      if (rsq < 1.0 && rsq != 0.0) again = 0;
    }
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
  if (seed_init <= 0) error->all("Invalid seed for Park random # generator");
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

  seed = hash & 0x7ffffff;

  // warm up the RNG

  for (i = 0; i < 5; i++) uniform();
  save = 0;
}

/* ---------------------------------------------------------------------- */

int RanPark::state()
{
  return seed;
}
