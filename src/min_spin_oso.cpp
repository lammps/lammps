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

/* ------------------------------------------------------------------------
   Contributing authors: Julien Tranchida (SNL)

   Please cite the related publication:
------------------------------------------------------------------------- */

#include <mpi.h>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include "min_spin_oso.h"
#include "universe.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "output.h"
#include "timer.h"
#include "error.h"
#include "modify.h"
#include "math_special.h"
#include "math_const.h"

using namespace LAMMPS_NS;
using namespace MathConst;

// EPS_ENERGY = minimum normalization for energy tolerance

#define EPS_ENERGY 1.0e-8

#define DELAYSTEP 5

void vm3(const double *m, const double *v, double *out);
void rodrigues_rotation(const double *upp_tr, double *out);

/* ---------------------------------------------------------------------- */

MinSpinOSO::MinSpinOSO(LAMMPS *lmp) : Min(lmp) {}

/* ---------------------------------------------------------------------- */

void MinSpinOSO::init()
{
  alpha_damp = 1.0;
  discrete_factor = 10.0;

  Min::init();

  dts = dt = update->dt;
  last_negative = update->ntimestep;
}

/* ---------------------------------------------------------------------- */

void MinSpinOSO::setup_style()
{
  double **v = atom->v;
  int nlocal = atom->nlocal;

  // check if the atom/spin style is defined

  if (!atom->sp_flag)
    error->all(FLERR,"min/spin_oso requires atom/spin style");

  for (int i = 0; i < nlocal; i++)
    v[i][0] = v[i][1] = v[i][2] = 0.0;
}

/* ---------------------------------------------------------------------- */

int MinSpinOSO::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0],"alpha_damp") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal fix_modify command");
    alpha_damp = force->numeric(FLERR,arg[1]);
    return 2;
  }
  if (strcmp(arg[0],"discrete_factor") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal fix_modify command");
    discrete_factor = force->numeric(FLERR,arg[1]);
    return 2;
  }
  return 0;
}

/* ----------------------------------------------------------------------
   set current vector lengths and pointers
   called after atoms have migrated
------------------------------------------------------------------------- */

void MinSpinOSO::reset_vectors()
{
  // atomic dof

  // size sp is 4N vector
  nvec = 4 * atom->nlocal;
  if (nvec) spvec = atom->sp[0];

  nvec = 3 * atom->nlocal;
  if (nvec) fmvec = atom->fm[0];

  if (nvec) xvec = atom->x[0];
  if (nvec) fvec = atom->f[0];
}

/* ----------------------------------------------------------------------
   minimization via damped spin dynamics
------------------------------------------------------------------------- */

int MinSpinOSO::iterate(int maxiter)
{
  bigint ntimestep;
  double fmdotfm;
  int flag, flagall;

  for (int iter = 0; iter < maxiter; iter++) {

    if (timer->check_timeout(niter))
      return TIMEOUT;

    ntimestep = ++update->ntimestep;
    niter++;

    // optimize timestep accross processes / replicas
    // need a force calculation for timestep optimization

    energy_force(0);
    dts = evaluate_dt();

    // apply damped precessional dynamics to the spins

    advance_spins(dts);

    eprevious = ecurrent;
    ecurrent = energy_force(0);
    neval++;

    //// energy tolerance criterion
    //// only check after DELAYSTEP elapsed since velocties reset to 0
    //// sync across replicas if running multi-replica minimization

    if (update->etol > 0.0 && ntimestep-last_negative > DELAYSTEP) {
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
        if (flagall == 0) return ETOL;
      }
    }

    // magnetic torque tolerance criterion
    // sync across replicas if running multi-replica minimization

    if (update->ftol > 0.0) {
      fmdotfm = fmnorm_sqr();
      if (update->multireplica == 0) {
        if (fmdotfm < update->ftol*update->ftol) return FTOL;
      } else {
        if (fmdotfm < update->ftol*update->ftol) flag = 0;
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

/* ----------------------------------------------------------------------
   evaluate max timestep
---------------------------------------------------------------------- */

double MinSpinOSO::evaluate_dt()
{
  double dtmax;
  double fmsq;
  double fmaxsqone,fmaxsqloc,fmaxsqall;
  int nlocal = atom->nlocal;
  double **fm = atom->fm;

  // finding max fm on this proc.

  fmsq = fmaxsqone = fmaxsqloc = fmaxsqall = 0.0;
  for (int i = 0; i < nlocal; i++) {
    fmsq = fm[i][0]*fm[i][0]+fm[i][1]*fm[i][1]+fm[i][2]*fm[i][2];
    fmaxsqone = MAX(fmaxsqone,fmsq);
  }

  // finding max fm on this replica

  fmaxsqloc = fmaxsqone;
  MPI_Allreduce(&fmaxsqone,&fmaxsqloc,1,MPI_DOUBLE,MPI_MAX,world);

  // finding max fm over all replicas, if necessary
  // this communicator would be invalid for multiprocess replicas

  fmaxsqall = fmaxsqloc;
  if (update->multireplica == 1) {
    fmaxsqall = fmaxsqloc;
    MPI_Allreduce(&fmaxsqloc,&fmaxsqall,1,MPI_DOUBLE,MPI_MAX,universe->uworld);
  }

  if (fmaxsqall == 0.0)
    error->all(FLERR,"Incorrect fmaxsqall calculation");

  // define max timestep by dividing by the
  // inverse of max frequency by discrete_factor

  dtmax = MY_2PI/(discrete_factor*sqrt(fmaxsqall));

  return dtmax;
}

/* ----------------------------------------------------------------------
   geometric damped advance of spins
---------------------------------------------------------------------- */

void MinSpinOSO::advance_spins(double dts)
{
  int nlocal = atom->nlocal;
  double **sp = atom->sp;
  double **fm = atom->fm;
  double tdampx, tdampy, tdampz;
  double f[3];  // upper triag. part of skew-symm. matr. to be exponented
  double rot_mat[9]; // exponential of a
  double s_new[3];

  // loop on all spins on proc.

  for (int i = 0; i < nlocal; i++) {

      // calc. damping torque
      tdampx = -alpha_damp*(fm[i][1]*sp[i][2] - fm[i][2]*sp[i][1]);
      tdampy = -alpha_damp*(fm[i][2]*sp[i][0] - fm[i][0]*sp[i][2]);
      tdampz = -alpha_damp*(fm[i][0]*sp[i][1] - fm[i][1]*sp[i][0]);

      // calculate rotation matrix
      f[0] = tdampz * dts;
      f[1] = -tdampy * dts;
      f[2] = tdampx * dts;
      rodrigues_rotation(f, rot_mat);

      // rotate spins
      vm3(rot_mat, sp[i], s_new);
      sp[i][0] = s_new[0];
      sp[i][1] = s_new[1];
      sp[i][2] = s_new[2];
  }

}

/* ----------------------------------------------------------------------
   compute and return ||mag. torque||_2^2
------------------------------------------------------------------------- */

double MinSpinOSO::fmnorm_sqr()
{
  int nlocal = atom->nlocal;
  double tx,ty,tz;
  double **sp = atom->sp;
  double **fm = atom->fm;

  // calc. magnetic torques

  double local_norm2_sqr = 0.0;
  for (int i = 0; i < nlocal; i++) {
    tx = (fm[i][1]*sp[i][2] - fm[i][2]*sp[i][1]);
    ty = (fm[i][2]*sp[i][0] - fm[i][0]*sp[i][2]);
    tz = (fm[i][0]*sp[i][1] - fm[i][1]*sp[i][0]);

    local_norm2_sqr += tx*tx + ty*ty + tz*tz;
  }

  // no extra atom calc. for spins

  if (nextra_atom)
    error->all(FLERR,"extra atom option not available yet");

  double norm2_sqr = 0.0;
  MPI_Allreduce(&local_norm2_sqr,&norm2_sqr,1,MPI_DOUBLE,MPI_SUM,world);

  return norm2_sqr;
}


void rodrigues_rotation(const double *upp_tr, double *out){

    /***
    * calculate 3x3 matrix exponential using Rodrigues' formula
    * (R. Murray, Z. Li, and S. Shankar Sastry,
    * A Mathematical Introduction to
    * Robotic Manipulation (1994), p. 28 and 30).
    *
    * upp_tr - vector x, y, z so that one calculate
    * U = exp(A) with A= [[0, x, y],
    *                     [-x, 0, z],
    *                     [-y, -z, 0]]
    ***/


    if (fabs(upp_tr[0]) < 1.0e-40 &&
        fabs(upp_tr[1]) < 1.0e-40 &&
        fabs(upp_tr[2]) < 1.0e-40){
        // if upp_tr is zero, return unity matrix
        int k;
        int m;
        for(k = 0; k < 3; k++){
            for(m = 0; m < 3; m++){
                if (m == k) out[3 * k + m] = 1.0;
                else out[3 * k + m] = 0.0;
            }
        }
        return;
    }

    double theta = sqrt(upp_tr[0] * upp_tr[0] +
                        upp_tr[1] * upp_tr[1] +
                        upp_tr[2] * upp_tr[2]);

    double A = cos(theta);
    double B = sin(theta);
    double D = 1 - A;
    double x = upp_tr[0]/theta;
    double y = upp_tr[1]/theta;
    double z = upp_tr[2]/theta;

    // diagonal elements of U
    out[0] = A + z * z * D;
    out[4] = A + y * y * D;
    out[8] = A + x * x * D;

    // off diagonal of U
    double s1 = -y * z *D;
    double s2 = x * z * D;
    double s3 = -x * y * D;

    double a1 = x * B;
    double a2 = y * B;
    double a3 = z * B;

    out[1] = s1 + a1;
    out[3] = s1 - a1;
    out[2] = s2 + a2;
    out[6] = s2 - a2;
    out[5] = s3 + a3;
    out[7] = s3 - a3;

}


void vm3(const double *m, const double *v, double *out){
    /***
     * out = vector^T x m,
     * m -- 3x3 matrix , v -- 3-d vector
     ***/

    int i;
    int j;

    for(i = 0; i < 3; i++){
        out[i] *= 0.0;
        for(j = 0; j < 3; j++){
            out[i] += *(m + 3 * j + i) * v[j];
        }
    }

}
