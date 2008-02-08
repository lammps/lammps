/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "math.h"
#include "string.h"
#include "fix_nve_dipole.h"
#include "atom.h"
#include "atom_vec.h"
#include "force.h"
#include "update.h"
#include "respa.h"
#include "error.h"

using namespace LAMMPS_NS;

// moment of inertia for a sphere

#define INERTIA 0.4

/* ---------------------------------------------------------------------- */

FixNVEDipole::FixNVEDipole(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 3) error->all("Illegal fix nve/dipole command");
  if (!atom->mu_flag || !atom->omega_flag || 
      !atom->torque_flag || !atom->avec->shape_type)
    error->all("Fix nve/dipole requires atom attributes "
	       "mu, omega, torque, shape");
  inertia = new double[atom->ntypes+1];
}

/* ---------------------------------------------------------------------- */

FixNVEDipole::~FixNVEDipole()
{
  delete [] inertia;
}

/* ---------------------------------------------------------------------- */

int FixNVEDipole::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  mask |= INITIAL_INTEGRATE_RESPA;
  mask |= FINAL_INTEGRATE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixNVEDipole::init()
{
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;

  if (strcmp(update->integrate_style,"respa") == 0)
    step_respa = ((Respa *) update->integrate)->step;

  // moment of inertia for each particle type

  double *mass = atom->mass;
  double **shape = atom->shape;

  for (int i = 1; i <= atom->ntypes; i++)
    inertia[i] = INERTIA * mass[i] * 0.25*shape[i][0]*shape[i][0];
}

/* ---------------------------------------------------------------------- */

void FixNVEDipole::initial_integrate(int vflag)
{
  double dtfm,msq,scale;

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **mu = atom->mu;
  double **omega = atom->omega;
  double **torque = atom->torque;
  double *mass = atom->mass;
  double *dipole = atom->dipole;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double g[3];

  // update v,x for all particles
  // update omega,mu for all dipoles
  // d_omega/dt = torque / inertia
  // d_mu/dt = omega cross mu
  // renormalize mu to dipole length

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      dtfm = dtf / mass[type[i]];
      v[i][0] += dtfm * f[i][0];
      v[i][1] += dtfm * f[i][1];
      v[i][2] += dtfm * f[i][2];
      x[i][0] += dtv * v[i][0];
      x[i][1] += dtv * v[i][1];
      x[i][2] += dtv * v[i][2];
      if (dipole[type[i]] > 0.0) {
	dtfm = dtf / inertia[type[i]];
	omega[i][0] += dtfm * torque[i][0];
	omega[i][1] += dtfm * torque[i][1];
	omega[i][2] += dtfm * torque[i][2];

	g[0] = mu[i][0] + dtv * (omega[i][1]*mu[i][2] - omega[i][2]*mu[i][1]);
	g[1] = mu[i][1] + dtv * (omega[i][2]*mu[i][0] - omega[i][0]*mu[i][2]);
	g[2] = mu[i][2] + dtv * (omega[i][0]*mu[i][1] - omega[i][1]*mu[i][0]);

	msq = g[0]*g[0] + g[1]*g[1] + g[2]*g[2];
	scale = dipole[type[i]]/sqrt(msq);
	mu[i][0] = g[0]*scale;
	mu[i][1] = g[1]*scale;
	mu[i][2] = g[2]*scale;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixNVEDipole::final_integrate()
{
  double dtfm;

  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double **torque = atom->torque;
  double *mass = atom->mass;
  double *dipole = atom->dipole;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  // update v for all particles
  // update omega for all dipoles

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      dtfm = dtf / mass[type[i]];
      v[i][0] += dtfm * f[i][0];
      v[i][1] += dtfm * f[i][1];
      v[i][2] += dtfm * f[i][2];
      if (dipole[type[i]] > 0.0) {
	dtfm = dtf / inertia[type[i]];
	omega[i][0] += dtfm * torque[i][0];
	omega[i][1] += dtfm * torque[i][1];
	omega[i][2] += dtfm * torque[i][2];
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixNVEDipole::initial_integrate_respa(int vflag, int ilevel, int flag)
{
  if (flag) return;             // only used by NPT,NPH

  dtv = step_respa[ilevel];
  dtf = 0.5 * step_respa[ilevel] * force->ftm2v;

  if (ilevel == 0) initial_integrate(vflag);
  else final_integrate();
}

/* ---------------------------------------------------------------------- */

void FixNVEDipole::final_integrate_respa(int ilevel)
{
  dtf = 0.5 * step_respa[ilevel] * force->ftm2v;
  final_integrate();
}

/* ---------------------------------------------------------------------- */

void FixNVEDipole::reset_dt()
{
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
}

