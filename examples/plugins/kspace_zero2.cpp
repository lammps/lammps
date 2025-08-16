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

/* ----------------------------------------------------------------------
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "kspace_zero2.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "pair.h"

#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

KSpaceZero2::KSpaceZero2(LAMMPS *lmp) : KSpace(lmp)
{
  ewaldflag = 1;
  pppmflag = 1;
  msmflag = 1;
  dispersionflag = 1;
  tip4pflag = 1;
  dipoleflag = 1;
  spinflag = 1;
}

/* ---------------------------------------------------------------------- */

void KSpaceZero2::settings(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR, "Illegal kspace_style {} command", force->kspace_style);

  accuracy_relative = fabs(utils::numeric(FLERR, arg[0], false, lmp));
  if (accuracy_relative > 1.0)
    error->all(FLERR, "Invalid relative accuracy {:g} for kspace_style {}", accuracy_relative,
               force->kspace_style);
  if ((narg != 0) && (narg != 1)) error->all(FLERR, "Illegal kspace_style command");
}

/* ---------------------------------------------------------------------- */

void KSpaceZero2::init()
{
  if (comm->me == 0) utils::logmesg(lmp, "Dummy KSpace initialization ...\n");

  // error checks

  if (force->pair == nullptr) error->all(FLERR, "KSpace solver requires a pair style");
  if (!atom->q_flag) error->all(FLERR, "KSpace style zero2 requires atom attribute q");

  // compute two charge force

  two_charge();

  int itmp;
  auto p_cutoff = (double *) force->pair->extract("cut_coul", itmp);
  if (p_cutoff == nullptr) error->all(FLERR, "KSpace style is incompatible with Pair style");
  double cutoff = *p_cutoff;

  qsum_qsq();

  accuracy = accuracy_relative * two_charge_force;

  // make initial g_ewald estimate
  // based on desired accuracy and real space cutoff
  // fluid-occupied volume used to estimate real-space error
  // zprd used rather than zprd_slab

  if (!gewaldflag) {
    if (accuracy <= 0.0) error->all(FLERR, "KSpace accuracy must be > 0");
    if (q2 == 0.0) error->all(FLERR, "Must use 'kspace_modify gewald' for uncharged system");
    g_ewald = accuracy * sqrt(atom->natoms * cutoff * domain->xprd * domain->yprd * domain->zprd) /
        (2.0 * q2);
    if (g_ewald >= 1.0)
      g_ewald = (1.35 - 0.15 * log(accuracy)) / cutoff;
    else
      g_ewald = sqrt(-log(g_ewald)) / cutoff;
  }

  if (comm->me == 0) utils::logmesg(lmp, "  G vector (1/distance) = {:.8g}\n", g_ewald);
}

/* ---------------------------------------------------------------------- */

void KSpaceZero2::setup()
{
  if (comm->me == 0) utils::logmesg(lmp, "Dummy KSpace setup\n");
}

/* ---------------------------------------------------------------------- */

void KSpaceZero2::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);
}
