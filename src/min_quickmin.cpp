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
#include "min_quickmin.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "output.h"
#include "timer.h"
#include "error.h"

using namespace LAMMPS_NS;

// same as in other min classes

enum{MAXITER,MAXEVAL,ETOL,FTOL,DOWNHILL,ZEROALPHA,ZEROFORCE,ZEROQUAD};

#define MIN(A,B) ((A) < (B)) ? (A) : (B)
#define MAX(A,B) ((A) > (B)) ? (A) : (B)

/* ---------------------------------------------------------------------- */

MinQuickmin::MinQuickmin(LAMMPS *lmp) : Min(lmp) {}

/* ---------------------------------------------------------------------- */

void MinQuickmin::init_style()
{
  dt = update->dt;
}

/* ---------------------------------------------------------------------- */

void MinQuickmin::setup_style()
{
  double **v = atom->v;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    v[i][0] = v[i][1] = v[i][2] = 0.0;
}

/* ----------------------------------------------------------------------
   set current vector lengths and pointers
   called after atoms have migrated
------------------------------------------------------------------------- */

void MinQuickmin::reset_vectors()
{
  // atomic dof

  nvec = 3 * atom->nlocal;
  if (nvec) xvec = atom->x[0];
  if (nvec) fvec = atom->f[0];
}

/* ---------------------------------------------------------------------- */

int MinQuickmin::iterate(int maxiter)
{
  int ntimestep;
  double vmax,vdotf,vdotfall,fdotf,fdotfall,scale;
  double dtvone,dtv,dtfm;

  alpha_final = 0.0;

  for (int iter = 0; iter < maxiter; iter++) {
    ntimestep = ++update->ntimestep;
    niter++;

    // zero velocity if anti-parallel to force
    // else project velocity in direction of force

    double **v = atom->v;
    double **f = atom->f;
    int nlocal = atom->nlocal;

    vdotf = 0.0;
    for (int i = 0; i < nlocal; i++)
      vdotf += v[i][0]*f[i][0] + v[i][1]*f[i][1] + v[i][2]*f[i][2];
    MPI_Allreduce(&vdotf,&vdotfall,1,MPI_DOUBLE,MPI_SUM,world);

    if (vdotfall < 0.0) {
      for (int i = 0; i < nlocal; i++)
	v[i][0] = v[i][1] = v[i][2] = 0.0;

    } else {
      fdotf = 0.0;
      for (int i = 0; i < nlocal; i++)
	fdotf += f[i][0]*f[i][0] + f[i][1]*f[i][1] + f[i][2]*f[i][2];
      MPI_Allreduce(&fdotf,&fdotfall,1,MPI_DOUBLE,MPI_SUM,world);
      if (fdotfall == 0.0) scale = 0.0;
      else scale = vdotfall/fdotfall;
      for (int i = 0; i < nlocal; i++) {
	v[i][0] = scale * f[i][0];
	v[i][1] = scale * f[i][1];
	v[i][2] = scale * f[i][2];
      }
    }

    // limit timestep so no particle moves further than dmax

    double *mass = atom->mass;
    int *type = atom->type;

    dtvone = dt;

    for (int i = 0; i < nlocal; i++) {
      vmax = MAX(fabs(v[i][0]),fabs(v[i][1]));
      vmax = MAX(vmax,fabs(v[i][2]));
      if (dtvone*vmax > dmax) dtvone = dmax/vmax;
    }
    MPI_Allreduce(&dtvone,&dtv,1,MPI_DOUBLE,MPI_MIN,world);

    // Euler integration step

    double **x = atom->x;
    
    for (int i = 0; i < nlocal; i++) {
      dtfm = dtv / mass[type[i]];
      x[i][0] += dtv * v[i][0];
      x[i][1] += dtv * v[i][1];
      x[i][2] += dtv * v[i][2];
      v[i][0] += dtfm * f[i][0];
      v[i][1] += dtfm * f[i][1];
      v[i][2] += dtfm * f[i][2];
    }
    
    eprevious = ecurrent;
    ecurrent = energy_force(0);
    neval++;

    // force tolerance criterion

    //fdotf = fnorm_sqr();
    //if (fdotf < update->ftol*update->ftol) return FTOL;

    // output for thermo, dump, restart files

    if (output->next == ntimestep) {
      timer->stamp();
      output->write(ntimestep);
      timer->stamp(TIME_OUTPUT);
    }
  }
  
  return MAXITER;
}
