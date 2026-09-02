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

#include "min_quickmin_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "error.h"
#include "force.h"
#include "output.h"
#include "timer.h"
#include "universe.h"
#include "update.h"

#include <cmath>

using namespace LAMMPS_NS;

// EPS_ENERGY = minimum normalization for energy tolerance

static constexpr double EPS_ENERGY = 1.0e-8;
static constexpr int DELAYSTEP = 5;

/* ---------------------------------------------------------------------- */

MinQuickMinKokkos::MinQuickMinKokkos(LAMMPS *lmp) : MinKokkos(lmp)
{
  atomKK = (AtomKokkos *) atom;
  kokkosable = 1;
}

/* ---------------------------------------------------------------------- */

void MinQuickMinKokkos::init()
{
  MinKokkos::init();

  dt = update->dt;
  last_negative = update->ntimestep;
}

/* ---------------------------------------------------------------------- */

void MinQuickMinKokkos::setup_style()
{
  atomKK->sync(Device,V_MASK);

  auto l_v = atomKK->k_v.view_device();

  Kokkos::parallel_for("min_quickmin/zero_v", atom->nlocal, LAMMPS_LAMBDA(const int i) {
    l_v(i,0) = l_v(i,1) = l_v(i,2) = 0.0;
  });

  atomKK->modified(Device,V_MASK);
}

/* ----------------------------------------------------------------------
   set current vector lengths and pointers
   called after atoms have migrated
------------------------------------------------------------------------- */

void MinQuickMinKokkos::reset_vectors()
{
  // atomic dof

  nvec = 3 * atom->nlocal;
  if (nvec) {
    auto d_x = atomKK->k_x.view_device();
    auto d_f = atomKK->k_f.view_device();
    xvec = DAT::t_kkfloat_1d(d_x.data(),nvec);
    fvec = DAT::t_kkacc_1d(d_f.data(),nvec);
  }
}

/* ----------------------------------------------------------------------
   minimization via QuickMin damped dynamics
------------------------------------------------------------------------- */

int MinQuickMinKokkos::iterate(int maxiter)
{
  bigint ntimestep;
  double vdotf,vdotfall,fdotf,fdotfall,scale;
  double dtvone,dtv,dtf;

  alpha_final = 0.0;

  for (int iter = 0; iter < maxiter; iter++) {

    if (timer->check_timeout(niter))
      return TIMEOUT;

    ntimestep = ++update->ntimestep;
    niter++;

    // re-acquire views and nlocal every iteration: energy_force() below can
    // trigger a reneighbor, which migrates atoms and may reallocate the
    // per-atom dual views

    atomKK->sync(Device,X_MASK|V_MASK|F_MASK|RMASS_MASK|TYPE_MASK);

    auto l_x = atomKK->k_x.view_device();
    auto l_v = atomKK->k_v.view_device();
    auto l_f = atomKK->k_f.view_device();
    auto l_rmass = atomKK->k_rmass.view_device();
    auto l_mass = atomKK->k_mass.view_device();
    auto l_type = atomKK->k_type.view_device();
    const int nlocal = atom->nlocal;

    // zero velocity if anti-parallel to force
    // else project velocity in direction of force

    vdotf = 0.0;
    Kokkos::parallel_reduce("min_quickmin/vdotf", nlocal, LAMMPS_LAMBDA(const int i, double &vdf) {
      vdf += static_cast<double>(static_cast<KK_ACC_FLOAT>(l_v(i,0))*l_f(i,0) +
                                 static_cast<KK_ACC_FLOAT>(l_v(i,1))*l_f(i,1) +
                                 static_cast<KK_ACC_FLOAT>(l_v(i,2))*l_f(i,2));
    },vdotf);
    MPI_Allreduce(&vdotf,&vdotfall,1,MPI_DOUBLE,MPI_SUM,world);

    // sum vdotf over replicas, if necessary
    // this communicator would be invalid for multiprocess replicas

    if (update->multireplica == 1) {
      vdotf = vdotfall;
      MPI_Allreduce(&vdotf,&vdotfall,1,MPI_DOUBLE,MPI_SUM,universe->uworld);
    }

    if (vdotfall < 0.0) {
      last_negative = ntimestep;

      Kokkos::parallel_for("min_quickmin/zero_v", nlocal, LAMMPS_LAMBDA(const int i) {
        l_v(i,0) = l_v(i,1) = l_v(i,2) = 0.0;
      });

    } else {
      fdotf = 0.0;
      Kokkos::parallel_reduce("min_quickmin/fdotf", nlocal, LAMMPS_LAMBDA(const int i, double &fdf) {
        fdf += static_cast<double>(l_f(i,0)*l_f(i,0) + l_f(i,1)*l_f(i,1) + l_f(i,2)*l_f(i,2));
      },fdotf);
      MPI_Allreduce(&fdotf,&fdotfall,1,MPI_DOUBLE,MPI_SUM,world);

      // sum fdotf over replicas, if necessary
      // this communicator would be invalid for multiprocess replicas

      if (update->multireplica == 1) {
        fdotf = fdotfall;
        MPI_Allreduce(&fdotf,&fdotfall,1,MPI_DOUBLE,MPI_SUM,universe->uworld);
      }

      if (fdotfall == 0.0) scale = 0.0;
      else scale = vdotfall/fdotfall;

      const KK_FLOAT l_scale = static_cast<KK_FLOAT>(scale);

      Kokkos::parallel_for("min_quickmin/project_v", nlocal, LAMMPS_LAMBDA(const int i) {
        l_v(i,0) = l_scale * static_cast<KK_FLOAT>(l_f(i,0));
        l_v(i,1) = l_scale * static_cast<KK_FLOAT>(l_f(i,1));
        l_v(i,2) = l_scale * static_cast<KK_FLOAT>(l_f(i,2));
      });
    }

    // limit timestep so no particle moves further than dmax

    dtvone = dt;
    const double l_dmax = dmax;

    Kokkos::parallel_reduce("min_quickmin/dtv_limit", nlocal, LAMMPS_LAMBDA(const int i, double &dtv_local) {
      const KK_FLOAT vmax = Kokkos::fmax(Kokkos::fabs(l_v(i,0)),
                             Kokkos::fmax(Kokkos::fabs(l_v(i,1)),Kokkos::fabs(l_v(i,2))));
      if (dtv_local * static_cast<double>(vmax) > l_dmax)
        dtv_local = l_dmax / static_cast<double>(vmax);
    },Kokkos::Min<double>(dtvone));
    dtvone = MIN(dtvone,dt);

    MPI_Allreduce(&dtvone,&dtv,1,MPI_DOUBLE,MPI_MIN,world);

    // min dtv over replicas, if necessary
    // this communicator would be invalid for multiprocess replicas

    if (update->multireplica == 1) {
      dtvone = dtv;
      MPI_Allreduce(&dtvone,&dtv,1,MPI_DOUBLE,MPI_MIN,universe->uworld);
    }

    dtf = dtv * force->ftm2v;

    // Euler integration step

    const KK_FLOAT l_dtv = static_cast<KK_FLOAT>(dtv);
    const KK_FLOAT l_dtf = static_cast<KK_FLOAT>(dtf);

    Kokkos::parallel_for("min_quickmin/integrate", nlocal, LAMMPS_LAMBDA(const int i) {
      const KK_FLOAT dtfm = l_dtf / (l_rmass.data() ? l_rmass(i) : l_mass(l_type(i)));
      l_x(i,0) += l_dtv * l_v(i,0);
      l_x(i,1) += l_dtv * l_v(i,1);
      l_x(i,2) += l_dtv * l_v(i,2);
      l_v(i,0) += dtfm * static_cast<KK_FLOAT>(l_f(i,0));
      l_v(i,1) += dtfm * static_cast<KK_FLOAT>(l_f(i,1));
      l_v(i,2) += dtfm * static_cast<KK_FLOAT>(l_f(i,2));
    });

    atomKK->modified(Device,X_MASK|V_MASK);

    eprevious = ecurrent;
    ecurrent = energy_force(0);
    neval++;

    // energy tolerance criterion
    // only check after DELAYSTEP elapsed since velocities reset to 0
    // sync across replicas if running multi-replica minimization

    if (update->etol > 0.0 && ntimestep-last_negative > DELAYSTEP) {
      const bool converged = (fabs(ecurrent-eprevious) <
          update->etol * 0.5*(fabs(ecurrent) + fabs(eprevious) + EPS_ENERGY));
      if (update->multireplica == 0) {
        if (converged) return ETOL;
      } else {
        int flag = converged ? 0 : 1;
        int flagall;
        MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,universe->uworld);
        if (flagall == 0) return ETOL;
      }
    }

    // force tolerance criterion
    // sync across replicas if running multi-replica minimization

    if (update->ftol > 0.0) {
      double fdotf_norm = 0.0;
      if (normstyle == MAX) fdotf_norm = fnorm_max();        // max force norm
      else if (normstyle == INF) fdotf_norm = fnorm_inf();   // inf force norm
      else if (normstyle == TWO) fdotf_norm = fnorm_sqr();   // Euclidean force 2-norm
      else error->all(FLERR,"Illegal min_modify command");

      const bool converged = (fdotf_norm < update->ftol*update->ftol);
      if (update->multireplica == 0) {
        if (converged) return FTOL;
      } else {
        int flag = converged ? 0 : 1;
        int flagall;
        MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,universe->uworld);
        if (flagall == 0) return FTOL;
      }
    }

    // output for thermo, dump, restart files

    if (output->next == ntimestep) {
      atomKK->sync(Host,ALL_MASK);

      timer->stamp();
      output->write(ntimestep);
      timer->stamp(Timer::OUTPUT);
    }
  }

  return MAXITER;
}
