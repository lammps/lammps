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

#include "fix_nve_omp.h"
#include "atom.h"
#include "force.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixNVEOMP::FixNVEOMP(LAMMPS *lmp, int narg, char **arg) :
  FixNVE(lmp, narg, arg) { }

/* ----------------------------------------------------------------------
   allow for both per-type and per-atom mass
------------------------------------------------------------------------- */

void FixNVEOMP::initial_integrate(int vflag)
{
  // update v and x of atoms in group

  double * const * const x = atom->x;
  double * const * const v = atom->v;
  const double * const * const f = atom->f;
  const int * const mask = atom->mask;
  const int nlocal = (igroup == atom->firstgroup) ? atom->nfirst : atom->nlocal;
  int i;

  if (atom->rmass) {
    const double * const rmass = atom->rmass;
#if defined (_OPENMP)
#pragma omp parallel for private(i) default(none) schedule(static)
#endif
    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        const double dtfm = dtf / rmass[i];
        v[i][0] += dtfm * f[i][0];
        v[i][1] += dtfm * f[i][1];
        v[i][2] += dtfm * f[i][2];
        x[i][0] += dtv * v[i][0];
        x[i][1] += dtv * v[i][1];
        x[i][2] += dtv * v[i][2];
      }

  } else {
    const double * const mass = atom->mass;
    const int * const type = atom->type;
#if defined (_OPENMP)
#pragma omp parallel for private(i) default(none) schedule(static)
#endif
    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        const double dtfm = dtf / mass[type[i]];
        v[i][0] += dtfm * f[i][0];
        v[i][1] += dtfm * f[i][1];
        v[i][2] += dtfm * f[i][2];
        x[i][0] += dtv * v[i][0];
        x[i][1] += dtv * v[i][1];
        x[i][2] += dtv * v[i][2];
      }
  }
}

/* ---------------------------------------------------------------------- */

void FixNVEOMP::final_integrate()
{
  // update v of atoms in group

  double * const * const v = atom->v;
  const double * const * const f = atom->f;
  const int * const mask = atom->mask;
  const int nlocal = (igroup == atom->firstgroup) ? atom->nfirst : atom->nlocal;
  int i;

  if (atom->rmass) {
    const double * const rmass = atom->rmass;
    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        const double dtfm = dtf / rmass[i];
        v[i][0] += dtfm * f[i][0];
        v[i][1] += dtfm * f[i][1];
        v[i][2] += dtfm * f[i][2];
      }

  } else {
    const double * const mass = atom->mass;
    const int * const type = atom->type;
    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        const double dtfm = dtf / mass[type[i]];
        v[i][0] += dtfm * f[i][0];
        v[i][1] += dtfm * f[i][1];
        v[i][2] += dtfm * f[i][2];
      }
  }
}

