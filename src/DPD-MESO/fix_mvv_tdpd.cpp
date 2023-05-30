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
   This is a time integrator for position, velocity and concentration (x,
   v and cc) using the modified velocity-Verlet (MVV) algorithm.
   Setting verlet = 0.5 recovers the standard velocity-Verlet algorithm.

   Contributing author: Zhen Li (Clemson University)
   Email: zli7@clemson.edu

   Please cite the related publication:
   Z. Li, A. Yazdani, A. Tartakovsky and G.E. Karniadakis. "Transport
   dissipative particle dynamics model for mesoscopic advection-diffusion
   -reaction problems". The Journal of Chemical Physics, 2015, 143: 014101.
------------------------------------------------------------------------- */

#include "fix_mvv_tdpd.h"

#include "atom.h"
#include "error.h"
#include "force.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixMvvTDPD::FixMvvTDPD(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (strcmp(style,"tdpd/verlet") != 0 && narg < 3)
    error->all(FLERR,"Illegal fix mvv/tdpd command");

  verlet = 0.5;
  if (narg > 3) verlet = utils::numeric(FLERR,arg[3],false,lmp);

  cc_species = atom->cc_species;

  dynamic_group_allow = 1;
  time_integrate = 1;
}

/* ---------------------------------------------------------------------- */

int FixMvvTDPD::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixMvvTDPD::init()
{
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
  if (!force->pair_match("^tdpd",0))
    error->all(FLERR, "Must use pair style tdpd with fix mvv/tdpd");
}

/* ----------------------------------------------------------------------
   allow for both per-type and per-atom mass
------------------------------------------------------------------------- */

void FixMvvTDPD::initial_integrate(int /*vflag*/)
{
  double dtfm;
  // update v and x and cc of atoms in group

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **cc = atom->cc;
  double **cc_flux = atom->cc_flux;
  double **vest = atom->vest;
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

     vest[i][0] = v[i][0] + dtfm * f[i][0];
     vest[i][1] = v[i][1] + dtfm * f[i][1];
     vest[i][2] = v[i][2] + dtfm * f[i][2];

     x[i][0] += dtv * vest[i][0];
     x[i][1] += dtv * vest[i][1];
     x[i][2] += dtv * vest[i][2];
     v[i][0] += 2.0 * verlet * dtfm * f[i][0];
     v[i][1] += 2.0 * verlet * dtfm * f[i][1];
     v[i][2] += 2.0 * verlet * dtfm * f[i][2];
     for (int k = 0; k < cc_species; k++)
        cc[i][k] += 0.5 * dtv * cc_flux[i][k];
  }
}

/* ---------------------------------------------------------------------- */

void FixMvvTDPD::final_integrate()
{
  double dtfm;

  // update v of atoms in group

  double **v = atom->v;
  double **f = atom->f;
  double **cc = atom->cc;
  double **cc_flux = atom->cc_flux;
  double **vest = atom->vest;
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

     v[i][0] = vest[i][0] + dtfm * f[i][0];
     v[i][1] = vest[i][1] + dtfm * f[i][1];
     v[i][2] = vest[i][2] + dtfm * f[i][2];
     for (int k = 0; k < cc_species; k++)
        cc[i][k] += 0.5 * dtv * cc_flux[i][k];
  }
}

/* ---------------------------------------------------------------------- */

void FixMvvTDPD::reset_dt()
{
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
}
