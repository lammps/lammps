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
   Contributing authors: Julien Guénolé, CNRS and
                         Erik Bitzek, FAU Erlangen-Nuernberg
------------------------------------------------------------------------- */

#include <cmath>
#include "min_fire.h"
#include "universe.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "output.h"
#include "timer.h"
#include "error.h"
#include "variable.h"
#include "modify.h"
#include "compute.h"
#include "domain.h"
#include "neighbor.h"
#include "comm.h"

using namespace LAMMPS_NS;

// EPS_ENERGY = minimum normalization for energy tolerance

#define EPS_ENERGY 1.0e-8

/* ---------------------------------------------------------------------- */

MinFire::MinFire(LAMMPS *lmp) : Min(lmp) {}

/* ---------------------------------------------------------------------- */

void MinFire::init()
{
  Min::init();

  // simple parameters validation

  if (tmax < tmin) error->all(FLERR,"tmax has to be larger than tmin");
  if (dtgrow < 1.0) error->all(FLERR,"dtgrow has to be larger than 1.0");
  if (dtshrink > 1.0) error->all(FLERR,"dtshrink has to be smaller than 1.0");

  dt = update->dt;
  dtmax = tmax * dt;
  dtmin = tmin * dt;
  alpha = alpha0;
  last_negative = ntimestep_start = update->ntimestep;
  vdotf_negatif = 0;
}

/* ---------------------------------------------------------------------- */

void MinFire::setup_style()
{
  double **v = atom->v;
  int nlocal = atom->nlocal;

  // print the parameters used within fire into the log

  const char *s1[] = {"eulerimplicit","verlet","leapfrog","eulerexplicit"};
  const char *s2[] = {"no","yes"};

  if (comm->me == 0 && logfile) {
      fprintf(logfile,"  Parameters for fire: \n"
      "    dmax delaystep dtgrow dtshrink alpha0 alphashrink tmax tmin "
      "   integrator halfstepback relaxbox relaxbox_mod relaxbox_rate ptol \n"
      "    %4g %9i %6g %8g %6g %11g %4g %4g %13s %12s \n",
      dmax, delaystep, dtgrow, dtshrink, alpha0, alphashrink, tmax, tmin, 
      s1[integrator], s2[halfstepback_flag]);
  }

  // initialize the velocities

  for (int i = 0; i < nlocal; i++)
    v[i][0] = v[i][1] = v[i][2] = 0.0;
}

/* ----------------------------------------------------------------------
   set current vector lengths and pointers
   called after atoms have migrated
------------------------------------------------------------------------- */

void MinFire::reset_vectors()
{
  // atomic dof

  nvec = 3 * atom->nlocal;
  if (nvec) xvec = atom->x[0];
  if (nvec) fvec = atom->f[0];
}

/* ---------------------------------------------------------------------- */

int MinFire::iterate(int maxiter)
{
  bigint ntimestep;
  double vmax,vdotf,vdotfall,vdotv,vdotvall,fdotf,fdotfall;
  double scale1,scale2;
  double dtvone,dtv,dtf,dtfm;
  int flag,flagall;

  alpha_final = 0.0;

  // Leap Frog integration initialization

  if (integrator == 2) {

    double **f = atom->f;
    double **v = atom->v;
    double *rmass = atom->rmass;
    double *mass = atom->mass;
    int *type = atom->type;
    int nlocal = atom->nlocal;

    energy_force(0);
    neval++;

    dtf = -0.5 * dt * force->ftm2v;

    if (rmass) {
      for (int i = 0; i < nlocal; i++) {
        dtfm = dtf / rmass[i];
        v[i][0] = dtfm * f[i][0];
        v[i][1] = dtfm * f[i][1];
        v[i][2] = dtfm * f[i][2];
      }
    } else {
      for (int i = 0; i < nlocal; i++) {
        dtfm = dtf / mass[type[i]];
        v[i][0] = dtfm * f[i][0];
        v[i][1] = dtfm * f[i][1];
        v[i][2] = dtfm * f[i][2];
      }
    }
  }

  for (int iter = 0; iter < maxiter; iter++) {

    if (timer->check_timeout(niter))
      return TIMEOUT;

    ntimestep = ++update->ntimestep;
    niter++;

    // pointers

    int nlocal = atom->nlocal;
    double **v = atom->v;
    double **f = atom->f;
    double **x = atom->x;
    double *rmass = atom->rmass;
    double *mass = atom->mass;
    int *type = atom->type;

   // vdotfall = v dot f

    vdotf = 0.0;
    for (int i = 0; i < nlocal; i++)
      vdotf += v[i][0]*f[i][0] + v[i][1]*f[i][1] + v[i][2]*f[i][2];
    MPI_Allreduce(&vdotf,&vdotfall,1,MPI_DOUBLE,MPI_SUM,world);

    // sum vdotf over replicas, if necessary
    // this communicator would be invalid for multiprocess replicas

    if (update->multireplica == 1) {
      vdotf = vdotfall;
      MPI_Allreduce(&vdotf,&vdotfall,1,MPI_DOUBLE,MPI_SUM,universe->uworld);
    }

    // if (v dot f) > 0:
    // v = (1-alpha) v + alpha |v| Fhat
    // |v| = length of v, Fhat = unit f
    // Only: (1-alpha) and alpha |v| Fhat is calculated here
    // the modificatin of v is made wihtin the integration, after v update
    // if more than delaystep since v dot f was negative:
    // increase timestep, update global timestep and decrease alpha

    if (vdotfall > 0.0) {
      vdotv = 0.0;
      vdotf_negatif = 0;
      for (int i = 0; i < nlocal; i++)
        vdotv += v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2];
      MPI_Allreduce(&vdotv,&vdotvall,1,MPI_DOUBLE,MPI_SUM,world);

      // sum vdotv over replicas, if necessary
      // this communicator would be invalid for multiprocess replicas

      if (update->multireplica == 1) {
        vdotv = vdotvall;
        MPI_Allreduce(&vdotv,&vdotvall,1,MPI_DOUBLE,MPI_SUM,universe->uworld);
      }

      fdotf = 0.0;
      for (int i = 0; i < nlocal; i++)
        fdotf += f[i][0]*f[i][0] + f[i][1]*f[i][1] + f[i][2]*f[i][2];
      MPI_Allreduce(&fdotf,&fdotfall,1,MPI_DOUBLE,MPI_SUM,world);

      // sum fdotf over replicas, if necessary
      // this communicator would be invalid for multiprocess replicas

      if (update->multireplica == 1) {
        fdotf = fdotfall;
        MPI_Allreduce(&fdotf,&fdotfall,1,MPI_DOUBLE,MPI_SUM,universe->uworld);
      }

      scale1 = 1.0 - alpha;
      if (fdotfall <= 1e-20) scale2 = 0.0;
      else scale2 = alpha * sqrt(vdotvall/fdotfall);

      if (ntimestep - last_negative > delaystep) {
        dt = MIN(dt*dtgrow,dtmax);
        update->dt = dt;
        alpha *= alphashrink;
      }

    // else (v dot f) <= 0
    // if more than delaystep since starting the relaxation:
    // reset alpha
    //    if dt*dtshrink > dtmin:
    //    decrease timestep
    //    update global timestep (for thermo output)
    // half step back within the dynamics: x(t) = x(t-0.5*dt)
    // reset velocities: v = 0

    } else {
      last_negative = ntimestep;
      int delayflag = 1;
      if (ntimestep - ntimestep_start < delaystep && delaystep_start_flag)
        delayflag = 0;
      if (delayflag) {
        alpha = alpha0;
        if (dt*dtshrink >= dtmin) {
          dt *= dtshrink;
          update->dt = dt;
        }
      }

      // stopping criterion while stuck in a local bassin of the PES

      vdotf_negatif++;
      if (max_vdotf_negatif > 0 && vdotf_negatif > max_vdotf_negatif)
        return MAXVDOTF;

      // inertia correction
      
      if (halfstepback_flag) {
        for (int i = 0; i < nlocal; i++) {
          x[i][0] -= 0.5 * dtv * v[i][0];
          x[i][1] -= 0.5 * dtv * v[i][1];
          x[i][2] -= 0.5 * dtv * v[i][2];
        }
      }

      for (int i = 0; i < nlocal; i++)
        v[i][0] = v[i][1] = v[i][2] = 0.0;
    }

    // limit timestep so no particle moves further than dmax

    dtvone = dt;

    for (int i = 0; i < nlocal; i++) {
      vmax = MAX(fabs(v[i][0]),fabs(v[i][1]));
      vmax = MAX(vmax,fabs(v[i][2]));
      if (dtvone*vmax > dmax) dtvone = dmax/vmax;
    }

    MPI_Allreduce(&dtvone,&dtv,1,MPI_DOUBLE,MPI_MIN,world);

    // min dtv over replicas, if necessary
    // this communicator would be invalid for multiprocess replicas

    if (update->multireplica == 1) {
      dtvone = dtv;
      MPI_Allreduce(&dtvone,&dtv,1,MPI_DOUBLE,MPI_MIN,universe->uworld);
    }

    // Dynamic integration scheme:
    // 0: semi-implicit Euler 
    // 1: velocity Verlet
    // 2: leapfrog (initial half step before the iteration loop)
    // 3: explicit Euler 

    // Semi-implicit Euler OR Leap Frog integration

    if (integrator == 0 || integrator == 2) {

      dtf = dtv * force->ftm2v; 

      if (rmass) {
        for (int i = 0; i < nlocal; i++) {
          dtfm = dtf / rmass[i];
          v[i][0] += dtfm * f[i][0];
          v[i][1] += dtfm * f[i][1];
          v[i][2] += dtfm * f[i][2];
          if (vdotfall > 0.0) {
            v[i][0] = scale1*v[i][0] + scale2*f[i][0];
            v[i][1] = scale1*v[i][1] + scale2*f[i][1];
            v[i][2] = scale1*v[i][2] + scale2*f[i][2];
          }
          x[i][0] += dtv * v[i][0];
          x[i][1] += dtv * v[i][1];
          x[i][2] += dtv * v[i][2];
        }
      } else {
        for (int i = 0; i < nlocal; i++) {
          dtfm = dtf / mass[type[i]];
          v[i][0] += dtfm * f[i][0];
          v[i][1] += dtfm * f[i][1];
          v[i][2] += dtfm * f[i][2];
          if (vdotfall > 0.0) {
            v[i][0] = scale1*v[i][0] + scale2*f[i][0];
            v[i][1] = scale1*v[i][1] + scale2*f[i][1];
            v[i][2] = scale1*v[i][2] + scale2*f[i][2];
          }
          x[i][0] += dtv * v[i][0];
          x[i][1] += dtv * v[i][1];
          x[i][2] += dtv * v[i][2];
        }
      }

      eprevious = ecurrent;
      ecurrent = energy_force(0);
      neval++;

    // Velocity Verlet integration

    } else if (integrator == 1) {

      dtf = 0.5 * dtv * force->ftm2v;

      if (rmass) {
        for (int i = 0; i < nlocal; i++) {
          dtfm = dtf / rmass[i];
          v[i][0] += dtfm * f[i][0];
          v[i][1] += dtfm * f[i][1];
          v[i][2] += dtfm * f[i][2];
          if (vdotfall > 0.0) {
            v[i][0] = scale1*v[i][0] + scale2*f[i][0];
            v[i][1] = scale1*v[i][1] + scale2*f[i][1];
            v[i][2] = scale1*v[i][2] + scale2*f[i][2];
          }
          x[i][0] += dtv * v[i][0];
          x[i][1] += dtv * v[i][1];
          x[i][2] += dtv * v[i][2];
        }
      } else {
        for (int i = 0; i < nlocal; i++) {
          dtfm = dtf / mass[type[i]];
          v[i][0] += dtfm * f[i][0];
          v[i][1] += dtfm * f[i][1];
          v[i][2] += dtfm * f[i][2];
          if (vdotfall > 0.0) {
            v[i][0] = scale1*v[i][0] + scale2*f[i][0];
            v[i][1] = scale1*v[i][1] + scale2*f[i][1];
            v[i][2] = scale1*v[i][2] + scale2*f[i][2];
          }
          x[i][0] += dtv * v[i][0];
          x[i][1] += dtv * v[i][1];
          x[i][2] += dtv * v[i][2];
        }
      }
      
      eprevious = ecurrent;
      ecurrent = energy_force(0);
      neval++;

      if (rmass) {
        for (int i = 0; i < nlocal; i++) {
          dtfm = dtf / rmass[i];
          v[i][0] += dtfm * f[i][0];
          v[i][1] += dtfm * f[i][1];
          v[i][2] += dtfm * f[i][2];
          }
      } else {
        for (int i = 0; i < nlocal; i++) {
          dtfm = dtf / mass[type[i]];
          v[i][0] += dtfm * f[i][0];
          v[i][1] += dtfm * f[i][1];
          v[i][2] += dtfm * f[i][2];
        }
      }

    // Standard Euler integration

    } else if (integrator == 3) {

      dtf = dtv * force->ftm2v; 

      if (rmass) {
        for (int i = 0; i < nlocal; i++) {
          dtfm = dtf / rmass[i];
          if (vdotfall > 0.0) {
            v[i][0] = scale1*v[i][0] + scale2*f[i][0];
            v[i][1] = scale1*v[i][1] + scale2*f[i][1];
            v[i][2] = scale1*v[i][2] + scale2*f[i][2];
          }
          x[i][0] += dtv * v[i][0];
          x[i][1] += dtv * v[i][1];
          x[i][2] += dtv * v[i][2];
          v[i][0] += dtfm * f[i][0];
          v[i][1] += dtfm * f[i][1];
          v[i][2] += dtfm * f[i][2];
        }
      } else {
        for (int i = 0; i < nlocal; i++) {
          dtfm = dtf / mass[type[i]];
          if (vdotfall > 0.0) {
            v[i][0] = scale1*v[i][0] + scale2*f[i][0];
            v[i][1] = scale1*v[i][1] + scale2*f[i][1];
            v[i][2] = scale1*v[i][2] + scale2*f[i][2];
          }
          x[i][0] += dtv * v[i][0];
          x[i][1] += dtv * v[i][1];
          x[i][2] += dtv * v[i][2];
          v[i][0] += dtfm * f[i][0];
          v[i][1] += dtfm * f[i][1];
          v[i][2] += dtfm * f[i][2];
        }
      }

      eprevious = ecurrent;
      ecurrent = energy_force(0);
      neval++;
    }

    // energy tolerance criterion
    // only check after delaystep elapsed since velocties reset to 0
    // sync across replicas if running multi-replica minimization
    // reset the timestep to the initial value

    if (update->etol > 0.0 && ntimestep-last_negative > delaystep) {
      if (update->multireplica == 0) {
        if (fabs(ecurrent-eprevious) <
            update->etol * 0.5*(fabs(ecurrent) + fabs(eprevious) + EPS_ENERGY))
          return ETOL;
      } else {
        if (fabs(ecurrent-eprevious) <
            update->etol * 0.5*(fabs(ecurrent) + fabs(eprevious) + EPS_ENERGY))
          flag = 0;
        else flag = 1;
        MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,universe->uworld);
        if (flagall == 0)
          return ETOL;
      }
    }

    // force tolerance criterion
    // sync across replicas if running multi-replica minimization
    // reset the timestep to the initial value

    fdotf = 0.0;
    if (update->ftol > 0.0) {
      if (normstyle == MAX) fdotf = fnorm_max();        // max force norm
      else if (normstyle == INF) fdotf = fnorm_inf();   // inf force norm
      else if (normstyle == TWO) fdotf = fnorm_sqr();   // Euclidean force 2-norm
      else error->all(FLERR,"Illegal min_modify command");
      if (update->multireplica == 0) {
        if (fdotf < update->ftol*update->ftol) return FTOL;
      } else {
        if (fdotf < update->ftol*update->ftol) flag = 0;
        else flag = 1;
        MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,universe->uworld);
        if (flagall == 0) return FTOL;
      }
    }

    // output for thermo, dump, restart files

    if (output->next == ntimestep) {
      timer->stamp();
      output->write(ntimestep);
      timer->stamp(Timer::OUTPUT);
    }
  }

  return MAXITER;
}
