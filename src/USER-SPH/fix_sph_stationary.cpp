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

#include "fix_sph_stationary.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixSPHStationary::FixSPHStationary(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg) {

  if ((atom->esph_flag != 1) || (atom->rho_flag != 1))
    error->all(FLERR,
        "Fix sph/stationary command requires atom_style with both energy and density, e.g. meso");

  if (narg != 3)
    error->all(FLERR,"Illegal number of arguments for fix sph/stationary command");

  time_integrate = 0;
}

/* ---------------------------------------------------------------------- */

int FixSPHStationary::setmask() {
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSPHStationary::init() {
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
}

/* ----------------------------------------------------------------------
 allow for both per-type and per-atom mass
 ------------------------------------------------------------------------- */

void FixSPHStationary::initial_integrate(int /*vflag*/) {

  double *rho = atom->rho;
  double *drho = atom->drho;
  double *esph = atom->esph;
  double *desph = atom->desph;

  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int i;

  if (igroup == atom->firstgroup)
    nlocal = atom->nfirst;

  for (i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      esph[i] += dtf * desph[i]; // half-step update of particle internal energy
      rho[i] += dtf * drho[i]; // ... and density
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixSPHStationary::final_integrate() {

  double *esph = atom->esph;
  double *desph = atom->desph;
  double *rho = atom->rho;
  double *drho = atom->drho;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup)
    nlocal = atom->nfirst;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      esph[i] += dtf * desph[i];
      rho[i] += dtf * drho[i];
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixSPHStationary::reset_dt() {
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
}
