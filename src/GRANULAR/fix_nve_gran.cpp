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

  if (atom->check_style("granular") == 0)
    error->all("Must use fix nve/gran with atom style granular");
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
  if (force->dimension == 3)
    dtfphi = 0.5 * update->dt * force->ftm2v / INERTIA3D;
  else
    dtfphi = 0.5 * update->dt * force->ftm2v / INERTIA2D;
}

/* ---------------------------------------------------------------------- */

void FixNVEGran::initial_integrate()
{
  double dtfm;

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **phix = atom->phix;
  double **phiv = atom->phiv;
  double **phia = atom->phia;
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
      dtfm = dtfphi / (radius[i]*radius[i]*rmass[i]);
      phiv[i][0] += dtfm * phia[i][0];
      phiv[i][1] += dtfm * phia[i][1];
      phiv[i][2] += dtfm * phia[i][2];
      phix[i][0] += dtv * phiv[i][0];
      phix[i][1] += dtv * phiv[i][1];
      phix[i][2] += dtv * phiv[i][2];
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixNVEGran::final_integrate()
{
  double dtfm;

  double **v = atom->v;
  double **f = atom->f;
  double **phiv = atom->phiv;
  double **phia = atom->phia;
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
      dtfm = dtfphi / (radius[i]*radius[i]*rmass[i]);
      phiv[i][0] += dtfm * phia[i][0];
      phiv[i][1] += dtfm * phia[i][1];
      phiv[i][2] += dtfm * phia[i][2];
    }
  }
}
