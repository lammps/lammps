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

#include "fix_nve_bpm_sphere.h"

#include "atom.h"
#include "domain.h"
#include "error.h"
#include "math_extra.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathExtra;

/* ---------------------------------------------------------------------- */

FixNVEBPMSphere::FixNVEBPMSphere(LAMMPS *_lmp, int narg, char **arg) : FixNVE(_lmp, narg, arg)
{
  if (narg < 3) error->all(FLERR, "Illegal fix nve/bpm/sphere command");

  time_integrate = 1;

  // process extra keywords
  // inertia = moment of inertia prefactor for sphere or disc

  inertia = 0.4;

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "disc") == 0) {
      inertia = 0.5;
      if (domain->dimension != 2)
        error->all(FLERR, "Fix nve/bpm/sphere disc requires 2d simulation");
      iarg++;
    } else
      error->all(FLERR, "Illegal fix nve/bpm/sphere command");
  }

  inv_inertia = 1.0 / inertia;

  // error checks

  if (!atom->quat_flag || !atom->sphere_flag)
    error->all(FLERR, "Fix nve/bpm/sphere requires atom style bpm/sphere");
}

/* ---------------------------------------------------------------------- */

void FixNVEBPMSphere::init()
{
  FixNVE::init();

  // check that all particles are finite-size spheres
  // no point particles allowed

  double *radius = atom->radius;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
      if (radius[i] == 0.0) error->one(FLERR, "Fix nve/bpm/sphere requires extended particles");
}

/* ---------------------------------------------------------------------- */

void FixNVEBPMSphere::initial_integrate(int /*vflag*/)
{
  double dtq, dtfm, dtirotate, particle_inertia;

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double **torque = atom->torque;
  double **quat = atom->quat;
  double *radius = atom->radius;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  // set timestep here since dt may have changed or come via rRESPA
  dtq = 0.5 * dtv;

  // update v,x,omega,quat for all particles
  // d_omega/dt = torque / inertia

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      dtfm = dtf / rmass[i];
      v[i][0] += dtfm * f[i][0];
      v[i][1] += dtfm * f[i][1];
      v[i][2] += dtfm * f[i][2];
      x[i][0] += dtv * v[i][0];
      x[i][1] += dtv * v[i][1];
      x[i][2] += dtv * v[i][2];

      particle_inertia = inertia * (radius[i] * radius[i] * rmass[i]);
      dtirotate = dtf / particle_inertia;
      omega[i][0] += dtirotate * torque[i][0];
      omega[i][1] += dtirotate * torque[i][1];
      omega[i][2] += dtirotate * torque[i][2];

      MathExtra::richardson_sphere(quat[i], omega[i], dtq);
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixNVEBPMSphere::final_integrate()
{
  double dtfm, dtirotate, particle_inertia;

  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double **torque = atom->torque;
  double *rmass = atom->rmass;
  double *radius = atom->radius;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  // update v,omega for all particles
  // d_omega/dt = torque / inertia

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      dtfm = dtf / rmass[i];
      v[i][0] += dtfm * f[i][0];
      v[i][1] += dtfm * f[i][1];
      v[i][2] += dtfm * f[i][2];

      particle_inertia = inertia * (radius[i] * radius[i] * rmass[i]);
      dtirotate = dtf / particle_inertia;
      omega[i][0] += dtirotate * torque[i][0];
      omega[i][1] += dtirotate * torque[i][1];
      omega[i][2] += dtirotate * torque[i][2];
    }
}
