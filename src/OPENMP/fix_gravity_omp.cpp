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
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "fix_gravity_omp.h"

#include "atom.h"
#include "update.h"
#include "input.h"
#include "modify.h"
#include "variable.h"

#include "omp_compat.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{CHUTE,SPHERICAL,GRADIENT,VECTOR};
enum{CONSTANT,EQUAL};

/* ---------------------------------------------------------------------- */

FixGravityOMP::FixGravityOMP(LAMMPS *lmp, int narg, char **arg) :
  FixGravity(lmp, narg, arg) { }

/* ---------------------------------------------------------------------- */

void FixGravityOMP::post_force(int /* vflag */)
{
  // update gravity due to variables

  if (varflag != CONSTANT) {
    modify->clearstep_compute();
    if (mstyle == EQUAL) magnitude = input->variable->compute_equal(mvar);
    if (vstyle == EQUAL) magnitude = input->variable->compute_equal(vvar);
    if (pstyle == EQUAL) magnitude = input->variable->compute_equal(pvar);
    if (tstyle == EQUAL) magnitude = input->variable->compute_equal(tvar);
    if (xstyle == EQUAL) magnitude = input->variable->compute_equal(xvar);
    if (ystyle == EQUAL) magnitude = input->variable->compute_equal(yvar);
    if (zstyle == EQUAL) magnitude = input->variable->compute_equal(zvar);
    modify->addstep_compute(update->ntimestep + 1);

    set_acceleration();
  }

  const double * const * const x = atom->x;
  double * const * const f = atom->f;
  double * const rmass = atom->rmass;
  double * const mass = atom->mass;
  int * const mask = atom->mask;
  int * const type = atom->type;
  const int nlocal = atom->nlocal;
  const double xacc_thr = xacc;
  const double yacc_thr = yacc;
  const double zacc_thr = zacc;

  eflag = 0;
  double grav = 0.0;

  if (rmass) {
#if defined(_OPENMP)
#pragma omp parallel for LMP_DEFAULT_NONE reduction(-:grav)
#endif
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        const double massone = rmass[i];
        f[i][0] += massone*xacc_thr;
        f[i][1] += massone*yacc_thr;
        f[i][2] += massone*zacc_thr;
        grav -= massone * (xacc_thr*x[i][0] + yacc_thr*x[i][1] + zacc_thr*x[i][2]);
      }
  } else {
#if defined(_OPENMP)
#pragma omp parallel for LMP_DEFAULT_NONE reduction(-:grav)
#endif
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        const double massone = mass[type[i]];
        f[i][0] += massone*xacc_thr;
        f[i][1] += massone*yacc_thr;
        f[i][2] += massone*zacc_thr;
        grav -= massone * (xacc_thr*x[i][0] + yacc_thr*x[i][1] + zacc_thr*x[i][2]);
      }
  }
  egrav = grav;
}

/* ---------------------------------------------------------------------- */

void FixGravityOMP::post_force_respa(int vflag, int ilevel, int /* iloop */)
{
  if (ilevel == ilevel_respa) post_force(vflag);
}
