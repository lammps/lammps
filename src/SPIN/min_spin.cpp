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

/* ------------------------------------------------------------------------
   Contributing authors: Julien Tranchida (SNL)

   Please cite the related publication:
------------------------------------------------------------------------- */

#include "min_spin.h"

#include "atom.h"
#include "error.h"
#include "math_const.h"
#include "output.h"
#include "timer.h"
#include "universe.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace MathConst;

// EPS_ENERGY = minimum normalization for energy tolerance

#define EPS_ENERGY 1.0e-8

#define DELAYSTEP 5

/* ---------------------------------------------------------------------- */

MinSpin::MinSpin(LAMMPS *lmp) : Min(lmp) {
  alpha_damp = 1.0;
  discrete_factor = 10.0;
}

/* ---------------------------------------------------------------------- */

void MinSpin::init()
{
  Min::init();

  dts = dt = update->dt;
  last_negative = update->ntimestep;
}

/* ---------------------------------------------------------------------- */

void MinSpin::setup_style()
{
  double **v = atom->v;
  int nlocal = atom->nlocal;

  // check if the atom/spin style is defined

  if (!atom->sp_flag)
    error->all(FLERR,"min/spin requires atom/spin style");

  for (int i = 0; i < nlocal; i++)
    v[i][0] = v[i][1] = v[i][2] = 0.0;
}

/* ---------------------------------------------------------------------- */

int MinSpin::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0],"alpha_damp") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal min_modify command");
    alpha_damp = utils::numeric(FLERR,arg[1],false,lmp);
    return 2;
  }
  if (strcmp(arg[0],"discrete_factor") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal min_modify command");
    discrete_factor = utils::numeric(FLERR,arg[1],false,lmp);
    return 2;
  }
  return 0;
}

/* ----------------------------------------------------------------------
   set current vector lengths and pointers
   called after atoms have migrated
------------------------------------------------------------------------- */

void MinSpin::reset_vectors()
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

int MinSpin::iterate(int maxiter)
{
  bigint ntimestep;
  double fmdotfm,fmsq;
  int flag,flagall;

  for (int iter = 0; iter < maxiter; iter++) {

    if (timer->check_timeout(niter))
      return TIMEOUT;

    ntimestep = ++update->ntimestep;
    niter++;

    // optimize timestep across processes / replicas
    // need a force calculation for timestep optimization

    if (iter == 0) energy_force(0);
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

    fmdotfm = fmsq = 0.0;
    if (update->ftol > 0.0) {
      if (normstyle == MAX) fmsq = max_torque();        // max torque norm
      else if (normstyle == INF) fmsq = inf_torque();   // inf torque norm
      else if (normstyle == TWO) fmsq = total_torque(); // Euclidean torque 2-norm
      else error->all(FLERR,"Illegal min_modify command");
      fmdotfm = fmsq*fmsq;
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

double MinSpin::evaluate_dt()
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

void MinSpin::advance_spins(double dts)
{
  int nlocal = atom->nlocal;
  double **sp = atom->sp;
  double **fm = atom->fm;
  double tdampx,tdampy,tdampz;
  double fm2,energy,dts2;
  double cp[3],g[3];

  dts2 = dts*dts;

  // loop on all spins on proc.

  for (int i = 0; i < nlocal; i++) {

    // calc. damping torque

    tdampx = -alpha_damp*(fm[i][1]*sp[i][2] - fm[i][2]*sp[i][1]);
    tdampy = -alpha_damp*(fm[i][2]*sp[i][0] - fm[i][0]*sp[i][2]);
    tdampz = -alpha_damp*(fm[i][0]*sp[i][1] - fm[i][1]*sp[i][0]);

    // apply advance algorithm (geometric, norm preserving)

    fm2 = (tdampx*tdampx+tdampy*tdampy+tdampz*tdampz);
    energy = (sp[i][0]*tdampx)+(sp[i][1]*tdampy)+(sp[i][2]*tdampz);

    cp[0] = tdampy*sp[i][2]-tdampz*sp[i][1];
    cp[1] = tdampz*sp[i][0]-tdampx*sp[i][2];
    cp[2] = tdampx*sp[i][1]-tdampy*sp[i][0];

    g[0] = sp[i][0]+cp[0]*dts;
    g[1] = sp[i][1]+cp[1]*dts;
    g[2] = sp[i][2]+cp[2]*dts;

    g[0] += (tdampx*energy-0.5*sp[i][0]*fm2)*0.5*dts2;
    g[1] += (tdampy*energy-0.5*sp[i][1]*fm2)*0.5*dts2;
    g[2] += (tdampz*energy-0.5*sp[i][2]*fm2)*0.5*dts2;

    g[0] /= (1+0.25*fm2*dts2);
    g[1] /= (1+0.25*fm2*dts2);
    g[2] /= (1+0.25*fm2*dts2);

    sp[i][0] = g[0];
    sp[i][1] = g[1];
    sp[i][2] = g[2];

    // renormalization (check if necessary)

    // no comm. to atoms with same tag
    // because no need for simplecticity
  }
}
