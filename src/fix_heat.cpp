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
#include "region.h"
#include "group.h"
#include "force.h"
#include "update.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixHeat::FixHeat(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal fix heat command");

  scalar_flag = 1;
  global_freq = 1;
  extscalar = 0;

  nevery = atoi(arg[3]);
  if (nevery <= 0) error->all(FLERR,"Illegal fix heat command");

  heat_input = atof(arg[4]);

  // optional args

  iregion = -1;
  idregion = NULL;

  int iarg = 5;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"region") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix heat command");
      iregion = domain->find_region(arg[iarg+1]);
      if (iregion == -1) error->all(FLERR,"Region ID for fix heat does not exist");
      int n = strlen(arg[iarg+1]) + 1;
      idregion = new char[n];
      strcpy(idregion,arg[iarg+1]);
      iarg += 2;
    } else error->all(FLERR,"Illegal fix heat command");
  }

  scale = 1.0;
}

/* ---------------------------------------------------------------------- */

FixHeat::~FixHeat()
{
  delete [] idregion;
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
  // set index and check validity of region

  if (iregion >= 0) {
    iregion = domain->find_region(idregion);
    if (iregion == -1) error->all(FLERR,"Region ID for fix heat does not exist");
  }

  // cannot have 0 atoms in group

  if (group->count(igroup) == 0) error->all(FLERR,"Fix heat group has no atoms");
  masstotal = group->mass(igroup);
}

/* ---------------------------------------------------------------------- */

void FixHeat::end_of_step()
{ 
  double heat,ke;
  double vsub[3],vcm[3];
  Region *region = NULL;
  if (iregion >= 0) region = domain->regions[iregion];
  
  if (iregion < 0) {
    heat = heat_input*nevery*update->dt*force->ftm2v;
    ke = group->ke(igroup)*force->ftm2v;
    group->vcm(igroup,masstotal,vcm);
  } else {
    masstotal = group->mass(igroup,iregion);
    if (masstotal == 0.0) error->all(FLERR,"Fix heat group has no atoms");
    heat = heat_input*nevery*update->dt*force->ftm2v;
    ke = group->ke(igroup,iregion)*force->ftm2v;
    group->vcm(igroup,masstotal,vcm,iregion);
  }

  double vcmsq = vcm[0]*vcm[0] + vcm[1]*vcm[1] + vcm[2]*vcm[2];
  double escale = (ke + heat - 0.5*vcmsq*masstotal)/(ke - 0.5*vcmsq*masstotal);
  if (escale < 0.0) error->all(FLERR,"Fix heat kinetic energy went negative");
  scale = sqrt(escale);

  vsub[0] = (scale-1.0) * vcm[0];
  vsub[1] = (scale-1.0) * vcm[1];
  vsub[2] = (scale-1.0) * vcm[2];

  double **x = atom->x;
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  if (iregion < 0) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
	v[i][0] = scale*v[i][0] - vsub[0];
	v[i][1] = scale*v[i][1] - vsub[1];
	v[i][2] = scale*v[i][2] - vsub[2];
      }
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2])) {
	v[i][0] = scale*v[i][0] - vsub[0];
	v[i][1] = scale*v[i][1] - vsub[1];
	v[i][2] = scale*v[i][2] - vsub[2];
      }
  }
}

/* ---------------------------------------------------------------------- */

double FixHeat::compute_scalar()
{
  return scale;
}
