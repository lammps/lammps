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

#include "fix_nve_noforce.h"

#include "atom.h"
#include "respa.h"
#include "update.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixNVENoforce::FixNVENoforce(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  if (narg != 3) utils::missing_cmd_args(FLERR, "fix nve/noforce", error);

  time_integrate = 1;
}

/* ---------------------------------------------------------------------- */

int FixNVENoforce::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= INITIAL_INTEGRATE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixNVENoforce::init()
{
  dtv = update->dt;

  if (utils::strmatch(update->integrate_style, "^respa"))
    step_respa = (dynamic_cast<Respa *>(update->integrate))->step;
}

/* ---------------------------------------------------------------------- */

void FixNVENoforce::initial_integrate(int /*vflag*/)
{
  double **x = atom->x;
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      x[i][0] += dtv * v[i][0];
      x[i][1] += dtv * v[i][1];
      x[i][2] += dtv * v[i][2];
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixNVENoforce::initial_integrate_respa(int vflag, int ilevel, int flag)
{
  if (flag) return;    // only used by NPT,NPH

  dtv = step_respa[ilevel];

  if (ilevel == 0) initial_integrate(vflag);
}

/* ---------------------------------------------------------------------- */

void FixNVENoforce::reset_dt()
{
  dtv = update->dt;
}
