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
   Contributing author: Axel Kohlmeyer (ICTP, Italy)
------------------------------------------------------------------------- */

#include "string.h"
#include "fix_oneway.h"
#include "atom.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "region.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{NONE=-1,X=0,Y=1,Z=2,XYZMASK=3,MINUS=4,PLUS=0};

/* ---------------------------------------------------------------------- */

FixOneWay::FixOneWay(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  direction = NONE;
  regionidx = 0;
  regionstr = NULL;

  if (narg < 6) error->all(FLERR,"Illegal fix oneway command");

  nevery = force->inumeric(FLERR,arg[3]);
  if (nevery < 1) error->all(FLERR,"Illegal nevery value for fix oneway");

  int len = strlen(arg[4]);
  regionstr = new char[len];
  strcpy(regionstr,arg[4]);

  if (strcmp(arg[5], "x") == 0) direction = X|PLUS;
  if (strcmp(arg[5], "X") == 0) direction = X|PLUS;
  if (strcmp(arg[5], "y") == 0) direction = Y|PLUS;
  if (strcmp(arg[5], "Y") == 0) direction = Y|PLUS;
  if (strcmp(arg[5], "z") == 0) direction = Z|PLUS;
  if (strcmp(arg[5], "Z") == 0) direction = Z|PLUS;
  if (strcmp(arg[5],"-x") == 0) direction = X|MINUS;
  if (strcmp(arg[5],"-X") == 0) direction = X|MINUS;
  if (strcmp(arg[5],"-y") == 0) direction = Y|MINUS;
  if (strcmp(arg[5],"-Y") == 0) direction = Y|MINUS;
  if (strcmp(arg[5],"-z") == 0) direction = Z|MINUS;
  if (strcmp(arg[5],"-Z") == 0) direction = Z|MINUS;

  global_freq = nevery;
}

/* ---------------------------------------------------------------------- */

FixOneWay::~FixOneWay()
{
  if (regionstr) delete[] regionstr;
}

/* ---------------------------------------------------------------------- */

int FixOneWay::setmask()
{
  return END_OF_STEP;
}

/* ---------------------------------------------------------------------- */

void FixOneWay::init()
{
  regionidx = domain->find_region(regionstr);
  if (regionidx < 0)
    error->warning(FLERR,"Region for fix oneway does not exist");
}

/* ---------------------------------------------------------------------- */

void FixOneWay::end_of_step()
{
  Region *region = domain->regions[regionidx];
  region->prematch();

  const int idx = direction & XYZMASK;
  const double * const * const x = atom->x;
  double * const * const v = atom->v;
  const int *mask = atom->mask;
  const int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; ++i) {
    if ((mask[i] & groupbit) && region->match(x[i][0],x[i][1],x[i][2])) {
      if (direction & MINUS) {
        if (v[i][idx] > 0.0) v[i][idx] = -v[i][idx];
      } else {
        if (v[i][idx] < 0.0) v[i][idx] = -v[i][idx];
      }
    }
  }
}

