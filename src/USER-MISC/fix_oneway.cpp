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

#include "fix_oneway.h"
#include "atom.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "region.h"

#include <string.h>
#include <stdlib.h>

using namespace LAMMPS_NS;
using namespace FixConst;

enum{NONE=-1,X=0,Y=1,Z=2,XYZMASK=3,MINUS=4,PLUS=0};

/* ---------------------------------------------------------------------- */

FixOneWay::FixOneWay(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  direction = NONE;
  regionidx = 0;
  regionstr = NULL;

  if (narg < 6)
    error->all(FLERR,"Illegal fix oneway command");

  int len = strlen(arg[3]);
  regionstr = new char[len];
  strcpy(regionstr,arg[3]);

  if (strcmp(arg[4], "x") == 0) direction = X|PLUS;
  if (strcmp(arg[4], "X") == 0) direction = X|PLUS;
  if (strcmp(arg[4], "y") == 0) direction = Y|PLUS;
  if (strcmp(arg[4], "Y") == 0) direction = Y|PLUS;
  if (strcmp(arg[4], "z") == 0) direction = Z|PLUS;
  if (strcmp(arg[4], "Z") == 0) direction = Z|PLUS;
  if (strcmp(arg[4],"-x") == 0) direction = X|MINUS;
  if (strcmp(arg[4],"-X") == 0) direction = X|MINUS;
  if (strcmp(arg[4],"-y") == 0) direction = Y|MINUS;
  if (strcmp(arg[4],"-Y") == 0) direction = Y|MINUS;
  if (strcmp(arg[4],"-z") == 0) direction = Z|MINUS;
  if (strcmp(arg[4],"-Z") == 0) direction = Z|MINUS;

  nevery = force->inumeric(FLERR,arg[5]);
  if (nevery < 1)
    error->all(FLERR,"Illegal nevery value for fix oneway");
  global_freq = nevery;
}

/* ---------------------------------------------------------------------- */

FixOneWay::~FixOneWay()
{
  if (regionstr) 
    delete[] regionstr;
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
  // nothing to do for non-existing region
  if (regionidx < 0) return;

  const double * const * const x = atom->x;
  double * const * const v = atom->v;
  const int *mask = atom->mask;
  Region *region = domain->regions[regionidx];
  const int nlocal = atom->nlocal;
  const int idx = direction & XYZMASK;

  for (int i = 0; i < nlocal; ++i) {
    if ((mask[i] & groupbit) 
        && region->match(x[i][0],x[i][1],x[i][2])) {
      if (direction & MINUS) {
        if (v[i][idx] > 0.0) v[i][idx] = -v[i][idx];
      } else {
        if (v[i][idx] < 0.0) v[i][idx] = -v[i][idx];
      }
    }
  }
}

