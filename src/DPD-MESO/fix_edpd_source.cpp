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

#include "fix_edpd_source.h"

#include "atom.h"
#include "domain.h"
#include "error.h"
#include "region.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

enum { SPHERE, CUBOID, REGION };

/* ---------------------------------------------------------------------- */

FixEDPDSource::FixEDPDSource(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), idregion(nullptr), region(nullptr)
{
  if (narg < 4) utils::missing_cmd_args(FLERR, "fix edpd/source", error);

  int iarg = 3;
  if (strcmp(arg[iarg], "sphere") == 0) {
    option = SPHERE;
    if (narg != 9) error->all(FLERR, "Illegal fix edpd/source command (4 args for sphere)");
    ++iarg;
    center[0] = utils::numeric(FLERR, arg[iarg++], false, lmp);
    center[1] = utils::numeric(FLERR, arg[iarg++], false, lmp);
    center[2] = utils::numeric(FLERR, arg[iarg++], false, lmp);
    radius = utils::numeric(FLERR, arg[iarg++], false, lmp);
    value = utils::numeric(FLERR, arg[iarg++], false, lmp);
  } else if (strcmp(arg[iarg], "cuboid") == 0) {
    option = CUBOID;
    if (narg != 11) error->all(FLERR, "Illegal fix edpd/source command (6 args for cuboid)");
    ++iarg;
    center[0] = utils::numeric(FLERR, arg[iarg++], false, lmp);
    center[1] = utils::numeric(FLERR, arg[iarg++], false, lmp);
    center[2] = utils::numeric(FLERR, arg[iarg++], false, lmp);
    dLx = utils::numeric(FLERR, arg[iarg++], false, lmp);
    dLy = utils::numeric(FLERR, arg[iarg++], false, lmp);
    dLz = utils::numeric(FLERR, arg[iarg++], false, lmp);
    value = utils::numeric(FLERR, arg[iarg++], false, lmp);
  } else if (strcmp(arg[iarg], "region") == 0) {
    option = REGION;
    if (narg != 6) error->all(FLERR, "Illegal fix edpd/source command (2 args for region)");
    ++iarg;
    if (!domain->get_region_by_id(arg[iarg]))
      error->all(FLERR, "Region {} for fix edpd/source does not exist", arg[iarg]);
    idregion = utils::strdup(arg[iarg++]);
    value = utils::numeric(FLERR, arg[iarg++], false, lmp);
  } else
    error->all(FLERR, "Illegal fix edpd/source command option {}", arg[iarg]);
}

/* ---------------------------------------------------------------------- */

FixEDPDSource::~FixEDPDSource()
{
  delete[] idregion;
}

/* ---------------------------------------------------------------------- */

int FixEDPDSource::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixEDPDSource::init()
{
  // set index and check validity of region

  if (idregion) {
    region = domain->get_region_by_id(idregion);
    if (!region) error->all(FLERR, "Region {} for fix edpd/source does not exist", idregion);
  }
}

/* ---------------------------------------------------------------------- */

void FixEDPDSource::post_force(int /*vflag*/)
{
  double **x = atom->x;
  double *edpd_flux = atom->edpd_flux;
  double *edpd_cv = atom->edpd_cv;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double drx, dry, drz, rsq;
  double radius_sq = radius * radius * radius;

  if (region) region->prematch();

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      if (option == SPHERE) {
        drx = x[i][0] - center[0];
        dry = x[i][1] - center[1];
        drz = x[i][2] - center[2];
        rsq = drx * drx + dry * dry + drz * drz;
        if (rsq < radius_sq) edpd_flux[i] += value * edpd_cv[i];
      } else if (option == CUBOID) {
        drx = x[i][0] - center[0];
        dry = x[i][1] - center[1];
        drz = x[i][2] - center[2];
        if (fabs(drx) <= 0.5 * dLx && fabs(dry) <= 0.5 * dLy && fabs(drz) <= 0.5 * dLz)
          edpd_flux[i] += value * edpd_cv[i];
      } else if (option == REGION) {
        if (region->match(x[i][0], x[i][1], x[i][2])) edpd_flux[i] += value * edpd_cv[i];
      }
    }
  }
}
