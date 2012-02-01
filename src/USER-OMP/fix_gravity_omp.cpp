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
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "fix_gravity_omp.h"
#include "atom.h"
#include "update.h"
#include "domain.h"
#include "respa.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{CHUTE,SPHERICAL,GRADIENT,VECTOR};

/* ---------------------------------------------------------------------- */

FixGravityOMP::FixGravityOMP(LAMMPS *lmp, int narg, char **arg) :
  FixGravity(lmp, narg, arg) { }

/* ---------------------------------------------------------------------- */

void FixGravityOMP::post_force(int vflag)
{
  // update direction of gravity vector if gradient style

  if (style == GRADIENT) {
    if (domain->dimension == 3) {
      double phi_current = degree2rad * 
	(phi + (update->ntimestep - time_origin)*dt*phigrad*360.0);
      double theta_current = degree2rad * 
	(theta + (update->ntimestep - time_origin)*dt*thetagrad*360.0);
      xgrav = sin(theta_current) * cos(phi_current);
      ygrav = sin(theta_current) * sin(phi_current);
      zgrav = cos(theta_current);
    } else {
      double theta_current = degree2rad * 
	(theta + (update->ntimestep - time_origin)*dt*thetagrad*360.0);
      xgrav = sin(theta_current);
      ygrav = cos(theta_current);
    }
    xacc = magnitude*xgrav;
    yacc = magnitude*ygrav;
    zacc = magnitude*zgrav;
  }

  const double * const * const x = atom->x;
  double * const * const f = atom->f;
  double * const rmass = atom->rmass;
  double * const mass = atom->mass;
  int * const mask = atom->mask;
  int * const type = atom->type;
  const int nlocal = atom->nlocal;
  const double xacc_thr = xacc;
  const double yacc_thr = yacc;
  const double zacc_thr = zacc;
  double massone;
  
  int i;
  eflag = 0;
  double grav = 0.0;

  if (rmass) {
#if defined(_OPENMP)
#pragma omp parallel for private(i,massone) default(none) reduction(-:grav)
#endif
    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
	massone = rmass[i];
	f[i][0] += massone*xacc_thr;
	f[i][1] += massone*yacc_thr;
	f[i][2] += massone*zacc_thr;
	grav -= massone * (xacc_thr*x[i][0] + yacc_thr*x[i][1] + zacc_thr*x[i][2]);
      }
  } else {
#if defined(_OPENMP)
#pragma omp parallel for private(i,massone) default(none) reduction(-:grav)
#endif
    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
	massone = mass[type[i]];
	f[i][0] += massone*xacc_thr;
	f[i][1] += massone*yacc_thr;
	f[i][2] += massone*zacc_thr;
	grav -= massone * (xacc_thr*x[i][0] + yacc_thr*x[i][1] + zacc_thr*x[i][2]);
      }
  }
  egrav = grav;
}

/* ---------------------------------------------------------------------- */

void FixGravityOMP::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
}

