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
   This is a time integrator for position, velocity and temperature (x,
   v and edpd_T) using the modified velocity-Verlet (MVV) algorithm.
   Setting verlet = 0.5 recovers the standard velocity-Verlet algorithm.

   Contributing author: Zhen Li (Brown University)
   Email: zhen_li@brown.edu

   Please cite the related publication:
   Z. Li, Y.-H. Tang, H. Lei, B. Caswell and G.E. Karniadakis. "Energy-
   conserving dissipative particle dynamics with temperature-dependent
   properties". Journal of Computational Physics, 2014, 265: 113-127.

   Z. Li, Y.-H. Tang , X. Li and G.E. Karniadakis. "Mesoscale modeling of
   phase transition dynamics of thermoresponsive polymers". Chemical
   Communications, 2015, 51: 11038-11040.
------------------------------------------------------------------------- */

#include "fix_mvv_edpd.h"
#include <cstring>
#include "atom.h"
#include "force.h"
#include "update.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixMvvEDPD::FixMvvEDPD(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (strcmp(style,"mvv/edpd") != 0 && narg < 3)
    error->all(FLERR,"Illegal fix mvv/edpd command");

  verlet = 0.5;
  if (narg > 3) verlet = utils::numeric(FLERR,arg[3],false,lmp);

  dynamic_group_allow = 1;
  time_integrate = 1;
}

/* ---------------------------------------------------------------------- */

int FixMvvEDPD::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixMvvEDPD::init()
{
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
}

/* ----------------------------------------------------------------------
   allow for both per-type and per-atom mass
------------------------------------------------------------------------- */

void FixMvvEDPD::initial_integrate(int /*vflag*/)
{
  double dtfm,dtT;
  // update v and x and cc of atoms in group

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double *edpd_temp = atom->edpd_temp;
  double *edpd_flux = atom->edpd_flux;
  double *edpd_cv = atom->edpd_cv;
  double **vest = atom->vest;
  double *vest_temp = atom->vest_temp;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  for (int i = 0; i < nlocal; i++)
  if (mask[i] & groupbit) {
     if (rmass) dtfm = dtf / rmass[i];
     else dtfm = dtf / mass[type[i]];

     dtT = 0.5 * dtv / edpd_cv[i];

     vest[i][0] = v[i][0] + dtfm * f[i][0];
     vest[i][1] = v[i][1] + dtfm * f[i][1];
     vest[i][2] = v[i][2] + dtfm * f[i][2];
     vest_temp[i] = edpd_temp[i] + dtT * edpd_flux[i];

     x[i][0] += dtv * vest[i][0];
     x[i][1] += dtv * vest[i][1];
     x[i][2] += dtv * vest[i][2];
     v[i][0] += 2.0 * verlet * dtfm * f[i][0];
     v[i][1] += 2.0 * verlet * dtfm * f[i][1];
     v[i][2] += 2.0 * verlet * dtfm * f[i][2];
     edpd_temp[i] += 2.0 * verlet * dtT * edpd_flux[i];
  }
}

/* ---------------------------------------------------------------------- */

void FixMvvEDPD::final_integrate()
{
  double dtfm, dtT;

  // update v and edpd_temp of atoms in group

  double **v = atom->v;
  double **f = atom->f;
  double *edpd_temp = atom->edpd_temp;
  double *edpd_flux = atom->edpd_flux;
  double *edpd_cv = atom->edpd_cv;
  double **vest = atom->vest;
  double *vest_temp = atom->vest_temp;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  for (int i = 0; i < nlocal; i++)
  if (mask[i] & groupbit) {
     if (rmass) dtfm = dtf / rmass[i];
     else dtfm = dtf / mass[type[i]];

     dtT = 0.5 * dtv / edpd_cv[i];

     v[i][0] = vest[i][0] + dtfm * f[i][0];
     v[i][1] = vest[i][1] + dtfm * f[i][1];
     v[i][2] = vest[i][2] + dtfm * f[i][2];
     edpd_temp[i] = vest_temp[i] + dtT * edpd_flux[i];
  }
}

/* ---------------------------------------------------------------------- */

void FixMvvEDPD::reset_dt()
{
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
}
