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

#include <cstdlib>
#include <cstring>
#include "fix_tdpd_source.h"
#include "atom.h"
#include "comm.h"
#include "update.h"
#include "modify.h"
#include "domain.h"
#include "lattice.h"
#include "input.h"
#include "variable.h"
#include "error.h"
#include "force.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixTDPDSource::FixTDPDSource(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (strcmp(style,"tdpd/source") != 0 && narg < 4)
    error->all(FLERR,"Illegal fix tdpd/source command");

  int iarg = 3;
  cc_index = force->inumeric(FLERR,arg[iarg++]);

  if (strcmp(arg[iarg],"sphere") == 0) option = 0;
  else if (strcmp(arg[iarg],"cuboid") == 0) option = 1;
  else error->all(FLERR,"Illegal fix tdpd/source command");
  iarg++;

  if(option == 0){
    if (narg != 10 ) error->all(FLERR,"Illegal fix tdpd/source command (5 args for sphere)");
    center[0] = force->numeric(FLERR,arg[iarg++]);
    center[1] = force->numeric(FLERR,arg[iarg++]);
    center[2] = force->numeric(FLERR,arg[iarg++]);
    radius  = force->numeric(FLERR,arg[iarg++]);
    value   = force->numeric(FLERR,arg[iarg++]);
  }
  else if(option == 1){
    if (narg != 12 ) error->all(FLERR,"Illegal fix tdpd/edpd command (7 args for cuboid)");
    center[0] = force->numeric(FLERR,arg[iarg++]);
    center[1] = force->numeric(FLERR,arg[iarg++]);
    center[2] = force->numeric(FLERR,arg[iarg++]);
    dLx = force->numeric(FLERR,arg[iarg++]);
    dLy = force->numeric(FLERR,arg[iarg++]);
    dLz = force->numeric(FLERR,arg[iarg++]);
    value = force->numeric(FLERR,arg[iarg++]);
  }
  else error->all(FLERR,"Illegal fix tdpd/source command");
}

/* ---------------------------------------------------------------------- */

FixTDPDSource::~FixTDPDSource()
{
}

/* ---------------------------------------------------------------------- */

int FixTDPDSource::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixTDPDSource::init()
{
}

/* ---------------------------------------------------------------------- */

void FixTDPDSource::post_force(int /*vflag*/)
{
  double **x = atom->x;
  double **cc_flux = atom->cc_flux;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double drx, dry, drz, rsq;
  double radius_sq = radius*radius*radius;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      if(option == 0){
        drx = x[i][0] - center[0];
        dry = x[i][1] - center[1];
        drz = x[i][2] - center[2];
        rsq = drx*drx + dry*dry + drz*drz;
        if(rsq < radius_sq)
          cc_flux[i][cc_index-1] += value;
      }
      else if(option == 1){
        drx = x[i][0] - center[0];
        dry = x[i][1] - center[1];
        drz = x[i][2] - center[2];
        if(fabs(drx) <= 0.5*dLx && fabs(dry) <= 0.5*dLy && fabs(drz) <= 0.5*dLz)
          cc_flux[i][cc_index-1] += value;
      }
    }
  }
}
