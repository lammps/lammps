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

#include "stdio.h"
#include "string.h"
#include "fix_nve_gran.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "domain.h"
#include "error.h"

using namespace LAMMPS_NS;

// moments of inertia for sphere and disk

#define INERTIA3D 0.4
#define INERTIA2D 0.5

/* ---------------------------------------------------------------------- */

FixNVEGran::FixNVEGran(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 3) error->all("Illegal fix nve/gran command");

  if (!atom->xorient_flag || !atom->omega_flag || !atom->torque_flag)
    error->all("Fix nve/gran requires atom attributes "
	       "xorient, omega, torque");
}

/* ---------------------------------------------------------------------- */

int FixNVEGran::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixNVEGran::init()
{
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
  if (domain->dimension == 3) dtfrotate = dtf / INERTIA3D;
  else dtfrotate = dtf / INERTIA2D;
}

/* ---------------------------------------------------------------------- */

void FixNVEGran::initial_integrate()
{
  double dtfm;

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **xorient = atom->xorient;
  double **omega = atom->omega;
  double **torque = atom->torque;
  double *rmass = atom->rmass;
  double *radius = atom->radius;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      dtfm = dtf / rmass[i];
      v[i][0] += dtfm * f[i][0];
      v[i][1] += dtfm * f[i][1];
      v[i][2] += dtfm * f[i][2];
      x[i][0] += dtv * v[i][0];
      x[i][1] += dtv * v[i][1];
      x[i][2] += dtv * v[i][2];

      dtfm = dtfrotate / (radius[i]*radius[i]*rmass[i]);
      omega[i][0] += dtfm * torque[i][0];
      omega[i][1] += dtfm * torque[i][1];
      omega[i][2] += dtfm * torque[i][2];
      xorient[i][0] += dtv * omega[i][0];
      xorient[i][1] += dtv * omega[i][1];
      xorient[i][2] += dtv * omega[i][2];
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixNVEGran::final_integrate()
{
  double dtfm;

  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double **torque = atom->torque;
  double *rmass = atom->rmass;
  double *radius = atom->radius;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      dtfm = dtf / rmass[i];
      v[i][0] += dtfm * f[i][0];
      v[i][1] += dtfm * f[i][1];
      v[i][2] += dtfm * f[i][2];

      dtfm = dtfrotate / (radius[i]*radius[i]*rmass[i]);
      omega[i][0] += dtfm * torque[i][0];
      omega[i][1] += dtfm * torque[i][1];
      omega[i][2] += dtfm * torque[i][2];
    }
  }
}
