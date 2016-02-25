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
   Contributing author: James Larentzos (U.S. Army Research Laboratory)
------------------------------------------------------------------------- */

#include <stdlib.h>
#include <string.h>
#include "fix_eos_cv.h"
#include "atom.h"
#include "error.h"
#include "force.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixEOScv::FixEOScv(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 4) error->all(FLERR,"Illegal fix eos/cv command");
  cvEOS = force->numeric(FLERR,arg[3]);
  if(cvEOS <= double(0.0)) error->all(FLERR,"EOS cv must be > 0.0");

  restart_peratom = 1;
  nevery = 1;
}

/* ---------------------------------------------------------------------- */

int FixEOScv::setmask()
{
  int mask = 0;
  mask |= POST_INTEGRATE;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixEOScv::init()
{
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  double *uCond = atom->uCond;
  double *uMech = atom->uMech;
  double *dpdTheta = atom->dpdTheta;

  if(this->restart_reset){
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
	dpdTheta[i] = (uCond[i]+uMech[i])/cvEOS;
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
	uCond[i] = double(0.5)*cvEOS*dpdTheta[i];
	uMech[i] = double(0.5)*cvEOS*dpdTheta[i];
      }
  }
}

/* ---------------------------------------------------------------------- */

void FixEOScv::post_integrate()
{
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  double *uCond = atom->uCond;
  double *uMech = atom->uMech;
  double *dpdTheta = atom->dpdTheta;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit){
      dpdTheta[i] = (uCond[i]+uMech[i])/cvEOS;
      if(dpdTheta[i] <= double(0.0)) 
	error->one(FLERR,"Internal temperature < zero");
    }
}

/* ---------------------------------------------------------------------- */

void FixEOScv::end_of_step()
{
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  double *uCond = atom->uCond;
  double *uMech = atom->uMech;
  double *dpdTheta = atom->dpdTheta;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit){
      dpdTheta[i] = (uCond[i]+uMech[i])/cvEOS;
      if(dpdTheta[i] <= double(0.0)) 
	error->one(FLERR,"Internal temperature < zero");
    }
}
