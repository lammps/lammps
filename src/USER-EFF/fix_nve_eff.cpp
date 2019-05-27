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
   Contributing author: Andres Jaramillo-Botero (Caltech)
------------------------------------------------------------------------- */

#include <cmath>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include "fix_nve_eff.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "domain.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixNVEEff::FixNVEEff(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (!atom->electron_flag)
    error->all(FLERR,"Fix nve/eff requires atom style electron");

  time_integrate = 1;
}

/* ---------------------------------------------------------------------- */

int FixNVEEff::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  mask |= INITIAL_INTEGRATE_RESPA;
  mask |= FINAL_INTEGRATE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixNVEEff::init()
{
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;

  if (strstr(update->integrate_style,"respa"))
    step_respa = ((Respa *) update->integrate)->step;
}

/* ----------------------------------------------------------------------
   allow for both per-type and per-atom mass
------------------------------------------------------------------------- */

void FixNVEEff::initial_integrate(int /*vflag*/)
{
  double dtfm;

  // update v,vr and x,radius of atoms in group

  double **x = atom->x;
  double *eradius = atom->eradius;
  double **v = atom->v;
  double *ervel = atom->ervel;
  double **f = atom->f;
  double *erforce = atom->erforce;
  double *mass = atom->mass;
  int *spin = atom->spin;
  double mefactor = domain->dimension/4.0;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  // x + dt * [v + 0.5 * dt * (f / m)];

  if (mass) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        dtfm = dtf / mass[type[i]];
        v[i][0] += dtfm * f[i][0];
        v[i][1] += dtfm * f[i][1];
        v[i][2] += dtfm * f[i][2];
        x[i][0] += dtv * v[i][0];
        x[i][1] += dtv * v[i][1];
        x[i][2] += dtv * v[i][2];
        if (abs(spin[i])==1) {
          ervel[i] += dtfm * erforce[i] / mefactor;
          eradius[i] += dtv * ervel[i];
        }
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixNVEEff::final_integrate()
{
  double dtfm;

  double **v = atom->v;
  double *ervel = atom->ervel;
  double *erforce = atom->erforce;
  double **f = atom->f;
  double *mass = atom->mass;
  int *spin = atom->spin;
  double mefactor = domain->dimension/4.0;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  // dyn_v[i] += m * dt * dyn_f[i];

  if (mass) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        dtfm = dtf / mass[type[i]];
        v[i][0] += dtfm * f[i][0];
        v[i][1] += dtfm * f[i][1];
        v[i][2] += dtfm * f[i][2];
        if (abs(spin[i])==1)
          ervel[i] += dtfm * erforce[i] / mefactor;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixNVEEff::initial_integrate_respa(int vflag, int ilevel, int /*iloop*/)
{
  dtv = step_respa[ilevel];
  dtf = 0.5 * step_respa[ilevel] * force->ftm2v;

  // innermost level - NVE update of v and x
  // all other levels - NVE update of v

  if (ilevel == 0) initial_integrate(vflag);
  else final_integrate();
}

/* ---------------------------------------------------------------------- */

void FixNVEEff::final_integrate_respa(int ilevel, int /*iloop*/)
{
  dtf = 0.5 * step_respa[ilevel] * force->ftm2v;
  final_integrate();
}

/* ---------------------------------------------------------------------- */

void FixNVEEff::reset_dt()
{
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
}
