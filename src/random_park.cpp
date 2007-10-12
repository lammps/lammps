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
   reset the seed based on atom position within box and ibase = caller seed
   combine 3 RNGs based on fractional position in box into one new seed
------------------------------------------------------------------------- */

void RanPark::reset(int ibase, double *coord)
{
  // for orthogonal box, lamda = fractional position in box
  // for triclinic box, convert to lamda coords

  double lamda[3];

  if (domain->triclinic == 0) {
    lamda[0] = (coord[0] - domain->boxlo[0]) / domain->prd[0];
    lamda[1] = (coord[1] - domain->boxlo[1]) / domain->prd[1];
    lamda[2] = (coord[2] - domain->boxlo[2]) / domain->prd[2];
  } else domain->x2lamda(coord,lamda);

  // seed 1,2,3 = combination of atom coord in each dim and user-input seed
  // map geometric extent into range of each of 3 RNGs
  // warm-up each RNG by calling it twice

  int seed1,seed2,seed3;

  seed1 = static_cast<int> (lamda[0] * IM1);
  seed1 = (seed1+ibase) % IM1;
  seed1 = (seed1*IA1+IC1) % IM1;
  seed1 = (seed1*IA1+IC1) % IM1;

  seed2 = static_cast<int> (lamda[1] * IM2);
  seed2 = (seed2+ibase) % IM2;
  seed2 = (seed2*IA2+IC2) % IM2;
  seed2 = (seed2*IA2+IC2) % IM2;

  seed3 = static_cast<int> (lamda[2] * IM3);
  seed3 = (seed3+ibase) % IM3;
  seed3 = (seed3*IA3+IC3) % IM3;
  seed3 = (seed3*IA3+IC3) % IM3;

  // fraction = 0-1 with giving each dim an equal weighting
  // use fraction to reset Park/Miller RNG seed
  // warm-up master RNG with new seed by calling it twice

  double fraction = 1.0*seed1/(3*IM1) + 1.0*seed2/(3*IM2) + 1.0*seed3/(3*IM3);
  seed = static_cast<int> (fraction*IM) + 1;
  if (seed >= IM) seed = IM-1;

  uniform();
  uniform();
}

/* ---------------------------------------------------------------------- */

int RanPark::state()
{
  return seed;
}
