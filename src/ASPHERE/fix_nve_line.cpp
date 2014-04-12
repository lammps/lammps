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

#include "math.h"
#include "stdio.h"
#include "string.h"
#include "fix_nve_line.h"
#include "atom.h"
#include "atom_vec_line.h"
#include "domain.h"
#include "math_const.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

#define INERTIA (1.0/12.0)     // moment of inertia prefactor for line segment

/* ---------------------------------------------------------------------- */

FixNVELine::FixNVELine(LAMMPS *lmp, int narg, char **arg) :
  FixNVE(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal fix nve/line command");

  time_integrate = 1;

  MINUSPI = -MY_PI;
  TWOPI = 2.0*MY_PI;
}

/* ---------------------------------------------------------------------- */

int FixNVELine::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  mask |= INITIAL_INTEGRATE_RESPA;
  mask |= FINAL_INTEGRATE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixNVELine::init()
{
  // error checks

  avec = (AtomVecLine *) atom->style_match("line");
  if (!avec) error->all(FLERR,"Fix nve/line requires atom style line");

  if (domain->dimension != 2)
    error->all(FLERR,"Fix nve/line can only be used for 2d simulations");

  // check that all particles are line segments
  // no point particles allowed

  int *line = atom->line;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      if (line[i] < 0) error->one(FLERR,"Fix nve/line requires line particles");
    }

  FixNVE::init();
}

/* ---------------------------------------------------------------------- */

void FixNVELine::initial_integrate(int vflag)
{
  double dtfm,dtirotate,length,theta;

  AtomVecLine::Bonus *bonus = avec->bonus;
  int *line = atom->line;
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double **torque = atom->torque;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  // set timestep here since dt may have changed or come via rRESPA

  double dtfrotate = dtf / INERTIA;

  // update v,x,omega,theta for all particles
  // d_omega/dt = torque / inertia
  // bound theta by -PI to PI

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      dtfm = dtf / rmass[i];
      v[i][0] += dtfm * f[i][0];
      v[i][1] += dtfm * f[i][1];
      x[i][0] += dtv * v[i][0];
      x[i][1] += dtv * v[i][1];

      length = bonus[line[i]].length;
      theta = bonus[line[i]].theta;
      dtirotate = dtfrotate / (length*length*rmass[i]);
      omega[i][2] += dtirotate * torque[i][2];
      theta += dtv * omega[i][2];
      while (theta <= MINUSPI) theta += TWOPI;
      while (theta > MY_PI) theta -= TWOPI;
      bonus[line[i]].theta = theta;
    }
}

/* ---------------------------------------------------------------------- */

void FixNVELine::final_integrate()
{
  double dtfm,dtirotate,length;

  AtomVecLine::Bonus *bonus = avec->bonus;
  int *line = atom->line;
  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double **torque = atom->torque;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  // set timestep here since dt may have changed or come via rRESPA

  double dtfrotate = dtf / INERTIA;

  // update v,omega for all particles
  // d_omega/dt = torque / inertia

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      dtfm = dtf / rmass[i];
      v[i][0] += dtfm * f[i][0];
      v[i][1] += dtfm * f[i][1];

      length = bonus[line[i]].length;
      dtirotate = dtfrotate / (length*length*rmass[i]);
      omega[i][2] += dtirotate * torque[i][2];
    }
}
