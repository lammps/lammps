// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
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

#include "fix_eos_cv.h"

#include "atom.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixEOScv::FixEOScv(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 4) error->all(FLERR,"Illegal fix eos/cv command");
  cvEOS = utils::numeric(FLERR,arg[3],false,lmp);
  if (cvEOS <= 0.0) error->all(FLERR,"EOS cv must be > 0.0");

  nevery = 1;

  if (atom->dpd_flag != 1)
    error->all(FLERR,"FixEOScv requires atom_style with internal temperature and energies (e.g. dpd)");
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

  if (this->restart_reset) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
        dpdTheta[i] = (uCond[i]+uMech[i])/cvEOS;
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        if (dpdTheta[i] <= 0.0)
          error->one(FLERR,"Internal temperature <= zero");
        uCond[i] = 0.0;
        uMech[i] = cvEOS*dpdTheta[i];
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
    if (mask[i] & groupbit) {
      dpdTheta[i] = (uCond[i]+uMech[i])/cvEOS;
      if (dpdTheta[i] <= 0.0)
        error->one(FLERR,"Internal temperature <= zero");
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
    if (mask[i] & groupbit) {
      dpdTheta[i] = (uCond[i]+uMech[i])/cvEOS;
      if (dpdTheta[i] <= 0.0)
        error->one(FLERR,"Internal temperature <= zero");
    }
}
