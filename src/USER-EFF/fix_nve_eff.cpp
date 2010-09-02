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

#include "stdio.h"
#include "string.h"
#include "fix_nve_eff.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "math.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixNVEEff::FixNVEEff(LAMMPS *lmp, int narg, char **arg) :
  FixNVE(lmp, narg, arg)
{
  // error check

  if (!atom->spin_flag || !atom->eradius_flag || 
      !atom->ervel_flag || !atom->erforce_flag) 
    error->all("Fix nve/eff requires atom attributes "
	       "spin, eradius, ervel, erforce");
}

/* ----------------------------------------------------------------------
   allow for both per-type and per-atom mass
------------------------------------------------------------------------- */

void FixNVEEff::initial_integrate(int vflag)
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
  double *rmass = atom->rmass;
  int *spin = atom->spin;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  // x + dt * [v + 0.5 * dt * (f / m)];
  
  if (rmass) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
	dtfm = dtf / rmass[i];
	v[i][0] += dtfm * f[i][0];
	v[i][1] += dtfm * f[i][1];
	v[i][2] += dtfm * f[i][2];
	x[i][0] += dtv * v[i][0];
	x[i][1] += dtv * v[i][1];
	x[i][2] += dtv * v[i][2];
	if (spin[i]) {
	  ervel[i] += dtfm * erforce[i] / 0.75;
	  eradius[i] += dtv * ervel[i];
	}
      }
    }
    
  } else {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
	dtfm = dtf / mass[type[i]];
	v[i][0] += dtfm * f[i][0];
	v[i][1] += dtfm * f[i][1];
	v[i][2] += dtfm * f[i][2];
	x[i][0] += dtv * v[i][0];
	x[i][1] += dtv * v[i][1];
	x[i][2] += dtv * v[i][2];
	if (spin[i]) {
	  ervel[i] += dtfm * erforce[i] / 0.75;
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
  double *rmass = atom->rmass;
  int *spin = atom->spin;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;
  
  // dyn_v[i] += m * dt * dyn_f[i];

  if (rmass) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
	dtfm = dtf / rmass[i];
	v[i][0] += dtfm * f[i][0];
	v[i][1] += dtfm * f[i][1];
	v[i][2] += dtfm * f[i][2];
	if (spin[i] != 0) 
	  ervel[i] += dtfm * erforce[i] / 0.75;
      }
    }
    
  } else {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
	dtfm = dtf / mass[type[i]];
	v[i][0] += dtfm * f[i][0];
	v[i][1] += dtfm * f[i][1];
	v[i][2] += dtfm * f[i][2];
	if (spin[i] != 0) 
	  ervel[i] += dtfm * erforce[i] / 0.75;
      }
    }
  }
}
