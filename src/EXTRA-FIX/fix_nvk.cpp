// clang-format off
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
     Contributing author: Efrem Braun (UC Berkeley)
------------------------------------------------------------------------- */

#include "fix_nvk.h"

#include "atom.h"
#include "error.h"
#include "force.h"
#include "math_extra.h"
#include "respa.h"
#include "update.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace FixConst;

FixNVK::FixNVK(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 3)
    error->all(FLERR,"Illegal fix nvk command");
  if (igroup) error->all(FLERR,"Fix nvk only supports group all");

  dynamic_group_allow = 1;
  time_integrate = 1;
}

/* ---------------------------------------------------------------------- */

int FixNVK::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  mask |= INITIAL_INTEGRATE_RESPA;
  mask |= FINAL_INTEGRATE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixNVK::init()
{
  dtv = update->dt;
  dtf = 0.5 * update->dt;

  if (utils::strmatch(update->integrate_style,"^respa")) {
    error->all(FLERR,"Fix nvk not yet enabled for RESPA");
    step_respa = (dynamic_cast<Respa *>(update->integrate))->step;
  }

  // compute initial kinetic energy
  // make better by calling compute_ke instead of copy/pasting code from compute_ke.cpp
  double pfactor = 0.5 * force->mvv2e;
  double **v = atom->v;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double ke = 0.0;
  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
        ke += rmass[i] * (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]);
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
        ke += mass[type[i]] *
          (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]);
  }
  MPI_Allreduce(&ke,&K_target,1,MPI_DOUBLE,MPI_SUM,world);
  K_target *= pfactor;
}

/* ----------------------------------------------------------------------
   allow for both per-type and per-atom mass
------------------------------------------------------------------------- */

void FixNVK::initial_integrate(int /*vflag*/)
{
  double sm;
  double a,b,sqtb,s,sdot;

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  // calculate s and sdot from Minary 2003, equations 4.12 and 4.13
  double a_local = 0.0;
  double b_local = 0.0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      a_local += MathExtra::dot3(f[i], v[i]);
      if (rmass) b_local += MathExtra::dot3(f[i], f[i]) / rmass[i];
      else b_local += MathExtra::dot3(f[i], f[i]) / mass[type[i]];
    }
  MPI_Allreduce(&a_local,&a,1,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(&b_local,&b,1,MPI_DOUBLE,MPI_SUM,world);
  a /= (2.0*K_target); // units of inverse time
  b /= (2.0*K_target * force->mvv2e); // units of inverse time squared
  sqtb = sqrt(b);
  s = a/b * (cosh(dtf*sqtb) - 1.0) + sinh(dtf*sqtb) / sqtb;
  sdot = a/b * sqtb * sinh(dtf*sqtb) + cosh(dtf*sqtb);

  // update v and x of atoms in group per Minary 2003, equations 4.15-4.17
  // note that equation 4.15, 4.17 should read p = (p+F*s/m)/sdot
  // note that equation 4.16 should read r = r + delt*p/m
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      if (rmass) sm = s / rmass[i];
      else sm = s / mass[type[i]];
      v[i][0] = (v[i][0] + f[i][0] * sm * force->ftm2v) / sdot;
      v[i][1] = (v[i][1] + f[i][1] * sm * force->ftm2v) / sdot;
      v[i][2] = (v[i][2] + f[i][2] * sm * force->ftm2v) / sdot;
      x[i][0] += dtv * v[i][0];
      x[i][1] += dtv * v[i][1];
      x[i][2] += dtv * v[i][2];
    }
}

/* ---------------------------------------------------------------------- */

void FixNVK::final_integrate()
{
  double sm;
  double a,b,sqtb,s,sdot;

  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  // calculate s and sdot from Minary 2003, equations 4.12 and 4.13
  double a_local = 0.0;
  double b_local = 0.0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      a_local += MathExtra::dot3(f[i], v[i]);
      if (rmass) b_local += MathExtra::dot3(f[i], f[i]) / rmass[i];
      else b_local += MathExtra::dot3(f[i], f[i]) / mass[type[i]];
    }
  MPI_Allreduce(&a_local,&a,1,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(&b_local,&b,1,MPI_DOUBLE,MPI_SUM,world);
  a /= (2.0*K_target); // units of inverse time
  b /= (2.0*K_target * force->mvv2e); // units of inverse time squared
  sqtb = sqrt(b);
  s = a/b * (cosh(dtf*sqtb) - 1.0) + sinh(dtf*sqtb) / sqtb;
  sdot = a/b * sqtb * sinh(dtf*sqtb) + cosh(dtf*sqtb);

  // update v and x of atoms in group per Minary 2003, equations 4.15-4.17
  // note that equation 4.15, 4.17 should read p = (p+F*s/m)/sdot
  // note that equation 4.16 should read r = r + delt*p/m
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      if (rmass) sm = s / rmass[i];
      else sm = s / mass[type[i]];
      v[i][0] = (v[i][0] + f[i][0] * sm * force->ftm2v) / sdot;
      v[i][1] = (v[i][1] + f[i][1] * sm * force->ftm2v) / sdot;
      v[i][2] = (v[i][2] + f[i][2] * sm * force->ftm2v) / sdot;
    }
}

/* ---------------------------------------------------------------------- */

void FixNVK::initial_integrate_respa(int vflag, int ilevel, int /*iloop*/)
{
  dtv = step_respa[ilevel];
  dtf = 0.5 * step_respa[ilevel];

  // innermost level - NVK update of v and x
  // all other levels - NVK update of v

  if (ilevel == 0) initial_integrate(vflag);
  else final_integrate();
}

/* ---------------------------------------------------------------------- */

void FixNVK::final_integrate_respa(int ilevel, int /*iloop*/)
{
  dtf = 0.5 * step_respa[ilevel];
  final_integrate();
}

/* ---------------------------------------------------------------------- */

void FixNVK::reset_dt()
{
  dtv = update->dt;
  dtf = 0.5 * update->dt;
}
