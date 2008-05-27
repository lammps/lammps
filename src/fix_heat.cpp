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
   Contributing author: Paul Crozier (SNL)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "fix_heat.h"
#include "atom.h"
#include "domain.h"
#include "group.h"
#include "force.h"
#include "update.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixHeat::FixHeat(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  if (narg < 4) error->all("Illegal fix heat command");

  nevery = atoi(arg[3]);
  if (nevery <= 0) error->all("Illegal fix heat command");

  heat = atof(arg[4]);
  heat *= nevery*update->dt*force->ftm2v;

  // cannot have 0 atoms in group

  if (group->count(igroup) == 0.0) error->all("Fix heat group has no atoms");
}

/* ---------------------------------------------------------------------- */

int FixHeat::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixHeat::init()
{
  masstotal = group->mass(igroup);
}

/* ---------------------------------------------------------------------- */

void FixHeat::end_of_step()
{ 
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double vsub[3],vcm[3];

  double ke = group->ke(igroup)*force->ftm2v;
  group->vcm(igroup,masstotal,vcm);
  double vcmsq = vcm[0]*vcm[0] + vcm[1]*vcm[1] + vcm[2]*vcm[2];
  double escale = (ke + heat - 0.5*vcmsq*masstotal)/(ke - 0.5*vcmsq*masstotal);
  if (escale < 0.0) error->all("Fix heat kinetic energy went negative");
  double r = sqrt(escale);

  vsub[0] = (r-1.0) * vcm[0];
  vsub[1] = (r-1.0) * vcm[1];
  vsub[2] = (r-1.0) * vcm[2];
  
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      v[i][0] = r*v[i][0] - vsub[0];
      v[i][1] = r*v[i][1] - vsub[1];
      v[i][2] = r*v[i][2] - vsub[2];
    }
}
