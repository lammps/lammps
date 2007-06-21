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

/* ----------------------------------------------------------------------
   Contributing author: Pieter in 't Veld (SNL)
------------------------------------------------------------------------- */

#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "fix_nvt_sllod.h"
#include "math_extra.h"
#include "atom.h"
#include "force.h"
#include "domain.h"
#include "comm.h"
#include "group.h"
#include "update.h"
#include "respa.h"
#include "modify.h"
#include "fix.h"
#include "fix_deform.h"
#include "compute.h"
#include "error.h"

using namespace LAMMPS_NS;

enum{NO_REMAP,X_REMAP,V_REMAP};                   // same as fix_deform.cpp

/* ---------------------------------------------------------------------- */

FixNVTSlodd::FixNVTSlodd(LAMMPS *lmp, int narg, char **arg) :
  FixNVT(lmp, narg, arg) {}

/* ---------------------------------------------------------------------- */

void FixNVTSlodd::init()
{
  FixNVT::init();

  // check fix deform remap settings

  int i;
  for (i = 0; i < modify->nfix; i++)
    if (strcmp(modify->fix[i]->style,"deform") == 0) {
      if (((FixDeform *) modify->fix[i])->remapflag != V_REMAP && 
	  comm->me == 0)
	error->warning("Using fix nvt/sllod with inconsistent fix deform remap option");
      break;
    }
  if (i == modify->nfix && comm->me == 0)
    error->warning("Using fix nvt/sllod with no fix deform defined");
}

/* ---------------------------------------------------------------------- */

void FixNVTSlodd::initial_integrate()
{
  double dtfm;

  double delta = update->ntimestep - update->firststep;
  delta /= update->nsteps;
  t_target = t_start + delta * (t_stop-t_start);

  // update eta_dot

  f_eta = t_freq*t_freq * (t_current/t_target - 1.0);
  eta_dot += f_eta*dthalf;
  eta_dot *= drag_factor;
  eta += dtv*eta_dot;
  factor = exp(-dthalf*eta_dot);

  // update vthermal and x of only atoms in NVT group
  // lamda = 0-1 triclinic lamda coords
  // vstream = streaming velocity = Hrate*lamda + Hratelo
  // vthermal = thermal velocity = v - vstream
  // vdelu = Hrate*Hinv*vthermal

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double *h_rate = domain->h_rate;
  double *h_ratelo = domain->h_ratelo;
  double h_two[6],lamda[3],vstream[3],vthermal[3],vdelu[3];
  MathExtra::multiply_shape_shape(h_rate,domain->h_inv,h_two);

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      dtfm = dtf / mass[type[i]];

      domain->x2lamda(x[i],lamda);
      vstream[0] = h_rate[0]*lamda[0] + h_rate[5]*lamda[1] + 
	h_rate[4]*lamda[2] + h_ratelo[0];
      vstream[1] = h_rate[1]*lamda[1] + h_rate[3]*lamda[2] + h_ratelo[1];
      vstream[2] = h_rate[2]*lamda[2] + h_ratelo[2];
      vthermal[0] = v[i][0] - vstream[0];
      vthermal[1] = v[i][1] - vstream[1];
      vthermal[2] = v[i][2] - vstream[2];
      vdelu[0] = h_two[0]*vthermal[0] + h_two[5]*vthermal[1] + 
	h_two[4]*vthermal[2];
      vdelu[1] = h_two[1]*vthermal[1] + h_two[3]*vthermal[2];
      vdelu[2] = h_two[2]*vthermal[2];

      v[i][0] = vstream[0] + 
	vthermal[0]*factor + dtfm*f[i][0] - dthalf*vdelu[0];
      v[i][1] = vstream[1] + 
	vthermal[1]*factor + dtfm*f[i][1] - dthalf*vdelu[1];
      v[i][2] = vstream[2] + 
	vthermal[2]*factor + dtfm*f[i][2] - dthalf*vdelu[2];
      x[i][0] += dtv * v[i][0];
      x[i][1] += dtv * v[i][1];
      x[i][2] += dtv * v[i][2];
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixNVTSlodd::final_integrate()
{
  double dtfm;

  // update vthermal of only atoms in NVT group
  // lamda = 0-1 triclinic lamda coords
  // vstream = streaming velocity = Hrate*lamda + Hratelo
  // vthermal = thermal velocity = v - vstream
  // vdelu = Hrate*Hinv*vthermal

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double *h_rate = domain->h_rate;
  double *h_ratelo = domain->h_ratelo;
  double h_two[6],lamda[3],vstream[3],vthermal[3],vdelu[3];
  MathExtra::multiply_shape_shape(h_rate,domain->h_inv,h_two);

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      dtfm = dtf / mass[type[i]];

      domain->x2lamda(x[i],lamda);
      vstream[0] = h_rate[0]*lamda[0] + h_rate[5]*lamda[1] + 
	h_rate[4]*lamda[2] + h_ratelo[0];
      vstream[1] = h_rate[1]*lamda[1] + h_rate[3]*lamda[2] + h_ratelo[1];
      vstream[2] = h_rate[2]*lamda[2] + h_ratelo[2];
      vthermal[0] = v[i][0] - vstream[0];
      vthermal[1] = v[i][1] - vstream[1];
      vthermal[2] = v[i][2] - vstream[2];
      vdelu[0] = h_two[0]*vthermal[0] + h_two[5]*vthermal[1] + 
	h_two[4]*vthermal[2];
      vdelu[1] = h_two[1]*vthermal[1] + h_two[3]*vthermal[2];
      vdelu[2] = h_two[2]*vthermal[2];

      v[i][0] = vstream[0] + 
	vthermal[0]*factor + dtfm*f[i][0] - dthalf*vdelu[0];
      v[i][1] = vstream[1] + 
	vthermal[1]*factor + dtfm*f[i][1] - dthalf*vdelu[1];
      v[i][2] = vstream[2] + 
	vthermal[2]*factor + dtfm*f[i][2] - dthalf*vdelu[2];
    }
  }

  // compute current T

  t_current = temperature->compute_scalar();

  // update eta_dot

  f_eta = t_freq*t_freq * (t_current/t_target - 1.0);
  eta_dot += f_eta*dthalf;
  eta_dot *= drag_factor;
}

/* ---------------------------------------------------------------------- */

void FixNVTSlodd::initial_integrate_respa(int ilevel, int flag)
{
  if (flag) return;             // only used by NPT,NPH

  // set timesteps by level

  double dtfm;
  dtv = step_respa[ilevel];
  dtf = 0.5 * step_respa[ilevel] * force->ftm2v;
  dthalf = 0.5 * step_respa[ilevel];

  // atom quantities

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  // outermost level - update eta_dot and apply to v
  // all other levels - NVE update of v

  if (ilevel == nlevels_respa-1) {
    double delta = update->ntimestep - update->firststep;
    delta /= update->nsteps;
    t_target = t_start + delta * (t_stop-t_start);
    
    // update eta_dot

    f_eta = t_freq*t_freq * (t_current/t_target - 1.0);
    eta_dot += f_eta*dthalf;
    eta_dot *= drag_factor;
    eta += dtv*eta_dot;
    factor = exp(-dthalf*eta_dot);
  } else factor = 1.0;

  // update v of only atoms in NVT group

  double *h_rate = domain->h_rate;
  double *h_ratelo = domain->h_ratelo;
  double h_two[6],lamda[3],vstream[3],vthermal[3],vdelu[3];
  MathExtra::multiply_shape_shape(h_rate,domain->h_inv,h_two);

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      dtfm = dtf / mass[type[i]];

      domain->x2lamda(x[i],lamda);
      vstream[0] = h_rate[0]*lamda[0] + h_rate[5]*lamda[1] + 
	h_rate[4]*lamda[2] + h_ratelo[0];
      vstream[1] = h_rate[1]*lamda[1] + h_rate[3]*lamda[2] + h_ratelo[1];
      vstream[2] = h_rate[2]*lamda[2] + h_ratelo[2];
      vthermal[0] = v[i][0] - vstream[0];
      vthermal[1] = v[i][1] - vstream[1];
      vthermal[2] = v[i][2] - vstream[2];
      vdelu[0] = h_two[0]*vthermal[0] + h_two[5]*vthermal[1] + 
	h_two[4]*vthermal[2];
      vdelu[1] = h_two[1]*vthermal[1] + h_two[3]*vthermal[2];
      vdelu[2] = h_two[2]*vthermal[2];

      v[i][0] = vstream[0] + 
	vthermal[0]*factor + dtfm*f[i][0] - dthalf*vdelu[0];
      v[i][1] = vstream[1] + 
	vthermal[1]*factor + dtfm*f[i][1] - dthalf*vdelu[1];
      v[i][2] = vstream[2] + 
	vthermal[2]*factor + dtfm*f[i][2] - dthalf*vdelu[2];
    }
  }

  // innermost level - also update x of only atoms in NVT group

  if (ilevel == 0) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
	x[i][0] += dtv * v[i][0];
	x[i][1] += dtv * v[i][1];
	x[i][2] += dtv * v[i][2];
      }
    }
  }
}
