/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "fix_graphics_periodic.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "graphics.h"
#include "group.h"
#include "memory.h"
#include "modify.h"
#include "respa.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixGraphicsPeriodic::FixGraphicsPeriodic(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), imgobjs(nullptr), imgparms(nullptr)
{
  if (narg < 4) utils::missing_cmd_args(FLERR, "fix graphics/periodic", error);

  // parse mandatory arg

  nevery = utils::inumeric(FLERR, arg[3], false, lmp);
  if (nevery <= 0) error->all(FLERR, 3, "Illegal fix graphics/periodic nevery value");
  global_freq = nevery;
  dynamic_group_allow = 1;

  // initialize and set defaults
  atomflag = true;
  bondflag = (atom->molecular == Atom::MOLECULAR);
  numobjs = 0;
  radius = -1.0;
  pxlo = pxhi = pylo = pyhi = pzlo = pzhi = 0;

  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "xlo") == 0) {
      pxlo = 1;
      ++iarg;
    } else if (strcmp(arg[iarg], "xhi") == 0) {
      pxhi = 1;
      ++iarg;
    } else if (strcmp(arg[iarg], "ylo") == 0) {
      pylo = 1;
      ++iarg;
    } else if (strcmp(arg[iarg], "yhi") == 0) {
      pyhi = 1;
      ++iarg;
    } else if (strcmp(arg[iarg], "zlo") == 0) {
      pzlo = 1;
      ++iarg;
    } else if (strcmp(arg[iarg], "zhi") == 0) {
      pzhi = 1;
      ++iarg;
    } else if (strcmp(arg[iarg], "radius") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix graphics/periodic radius", error);
      if (strcmp(arg[iarg + 1], "auto") == 0) {
        radius = -1.0;
      } else {
        radius = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
        if (radius <= 0.0) error->all(FLERR, iarg, "Fix graphics/periodic radius must be > 0");
      }
      iarg += 2;
    } else if (strcmp(arg[iarg], "atoms") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix graphics/periodic atoms", error);
      atomflag = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "bonds") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix graphics/periodic bonds", error);
      bondflag = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else {
      error->all(FLERR, iarg, "Unknown fix graphics/periodic keyword: {}", arg[iarg]);
    }
  }

  // error checks
  if (bondflag && (atom->molecular != Atom::MOLECULAR))
    error->all(FLERR, Error::NOLASTLINE,
               "Cannot display periodic images of bonds with non-molecular system");
}

/* ---------------------------------------------------------------------- */

FixGraphicsPeriodic::~FixGraphicsPeriodic()
{
  memory->destroy(imgobjs);
  memory->destroy(imgparms);
}

/* ---------------------------------------------------------------------- */

int FixGraphicsPeriodic::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixGraphicsPeriodic::setup(int)
{
  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixGraphicsPeriodic::end_of_step()
{
  // determine number of replica needed
  int nx = 1 + pxlo + pxhi;
  int ny = 1 + pylo + pyhi;
  int nz = 1 + pzlo + pzhi;
  int nrep = nx * ny * nz - 1;

  const auto *const *const x = atom->x;
  const auto *const mask = atom->mask;
  const auto *const type = atom->type;
  const auto *const num_bond = atom->num_bond;
  const auto *const *const bond_atom = atom->bond_atom;

  const auto nlocal = atom->nlocal;
  const auto *const prd = domain->prd;

  // count number of replica objects needed
  int n = 0;
  for (int i = 0; i < nlocal; ++i) {
    if (mask[i] & groupbit) {
      if (atomflag) { ++n; }
      if (bondflag) {
        for (int j = 0; j < num_bond[i]; ++j) {
          int m = atom->map(bond_atom[i][j]);
          m = domain->closest_image(i, m);
          if (m < 0) continue;
          if (mask[m] & groupbit) ++n;
        }
      }
    }
  }

  numobjs = n * nrep;
  memory->destroy(imgobjs);
  memory->destroy(imgparms);
  memory->create(imgobjs, numobjs, "fix_graphics_periodic:imgobjs");
  memory->create(imgparms, numobjs, 8, "fix_graphics_periodic:imgparms");

  n = 0;
  for (int i = 0; i < nlocal; ++i) {
    if (mask[i] & groupbit) {
      for (int ix = pxlo ? -1 : 0; ix <= pxhi; ++ix) {
        for (int iy = pylo ? -1 : 0; iy <= pyhi; ++iy) {
          for (int iz = pzlo ? -1 : 0; iz <= pzhi; ++iz) {
            if ((ix == 0) && (iy == 0) && (iz == 0)) continue;
            if (atomflag) {
              imgobjs[n] = Graphics::SPHERE;
              imgparms[n][0] = type[i];
              imgparms[n][1] = x[i][0] + ix * prd[0];
              imgparms[n][2] = x[i][1] + iy * prd[1];
              imgparms[n][3] = x[i][2] + iz * prd[2];
              imgparms[n][4] = 2.0 * radius;
              ++n;
            }
            if (bondflag) {
              for (int j = 0; j < num_bond[i]; ++j) {
                int m = atom->map(bond_atom[i][j]);
                m = domain->closest_image(i, m);
                if (m < 0) continue;
                if (mask[m] & groupbit) {
                  imgobjs[n] = Graphics::BOND;
                  imgparms[n][0] = type[i];
                  imgparms[n][1] = type[m];
                  imgparms[n][2] = x[i][0] + ix * prd[0];
                  imgparms[n][3] = x[i][1] + iy * prd[1];
                  imgparms[n][4] = x[i][2] + iz * prd[2];
                  imgparms[n][5] = x[m][0] + ix * prd[0];
                  imgparms[n][6] = x[m][1] + iy * prd[1];
                  imgparms[n][7] = x[m][2] + iz * prd[2];
                  ++n;
                }
              }
            }
          }
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   provide graphics information to dump image
------------------------------------------------------------------------- */

int FixGraphicsPeriodic::image(int *&objs, double **&parms)
{
  objs = imgobjs;
  parms = imgparms;
  return numobjs;
}
