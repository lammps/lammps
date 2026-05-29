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
   Contributing author: Mitch Murphy (alphataubio at gmail)
------------------------------------------------------------------------- */

#include "min_fire_kokkos.h"
#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "output.h"
#include "timer.h"
#include "universe.h"
#include "update.h"

using namespace LAMMPS_NS;

static constexpr KK_FLOAT EPS_ENERGY = 1.0e-8;

MinFireKokkos::MinFireKokkos(LAMMPS *lmp) : MinKokkos(lmp) {
  atomKK = (AtomKokkos *) atom;
  kokkosable = 1;
}

void MinFireKokkos::init() {
  MinKokkos::init();

  if (tmax < tmin) error->all(FLERR, "tmax has to be larger than tmin");
  if (dtgrow < 1.0) error->all(FLERR, "dtgrow has to be larger than 1.0");
  if (dtshrink > 1.0) error->all(FLERR, "dtshrink has to be smaller than 1.0");

  dt = update->dt;
  dtmax = tmax * dt;
  dtmin = tmin * dt;
  alpha = alpha0;
  last_negative = ntimestep_start = update->ntimestep;
  vdotf_negatif = 0;

  // bugfix for multiprocess replicas
  if (update->multireplica == 1) mpi_comm = universe->uworld;
  else mpi_comm = world;
}

void MinFireKokkos::setup_style() {
  atomKK->sync(Device, V_MASK);
  auto l_v = atomKK->k_v.view_device();
  int nlocal = atom->nlocal;

  // print the parameters used within fire/abcfire into the log

  const char *integrator_names[] = {"eulerimplicit", "verlet", "leapfrog", "eulerexplicit"};
  const char *yesno[] = {"no", "yes"};

  if (comm->me == 0)
    utils::logmesg(lmp,
                   "  Parameters for {}:\n"
                   "    {:^5} {:^9} {:^6} {:^8} {:^6} {:^11} {:^4} {:^4} {:^14} {:^12} {:^11}\n"
                   "    {:^5} {:^9} {:^6} {:^8} {:^6} {:^11} {:^4} {:^4} {:^14} {:^12} {:^11}\n",
                   update->minimize_style, "dmax", "delaystep", "dtgrow", "dtshrink", "alpha0",
                   "alphashrink", "tmax", "tmin", "integrator", "halfstepback", "abcfire", dmax,
                   delaystep, dtgrow, dtshrink, alpha0, alphashrink, tmax, tmin,
                   integrator_names[integrator], yesno[halfstepback_flag], yesno[abcflag]);

  // initialize the velocities

  Kokkos::parallel_for("min_fire/zero_v", nlocal, LAMMPS_LAMBDA(const int i) {
    l_v(i,0) = l_v(i,1) = l_v(i,2) = 0.0;
  });
  atomKK->modified(Device, V_MASK);
}

void MinFireKokkos::reset_vectors() {
  nvec = 3 * atom->nlocal;
  if (nvec) {
    auto d_x = atomKK->k_x.view_device();
    auto d_f = atomKK->k_f.view_device();
    xvec = DAT::t_kkfloat_1d(d_x.data(), nvec);
    fvec = DAT::t_kkacc_1d(d_f.data(), nvec);
  }
}

int MinFireKokkos::iterate(int maxiter) {
  if (integrator == EULERIMPLICIT) {
    return abcflag ? run_iterate<EULERIMPLICIT, true>(maxiter) : run_iterate<EULERIMPLICIT, false>(maxiter);
  } else if (integrator == VERLET) {
    return abcflag ? run_iterate<VERLET, true>(maxiter) : run_iterate<VERLET, false>(maxiter);
  } else if (integrator == LEAPFROG) {
    return abcflag ? run_iterate<LEAPFROG, true>(maxiter) : run_iterate<LEAPFROG, false>(maxiter);
  } else if (integrator == EULEREXPLICIT) {
    return abcflag ? run_iterate<EULEREXPLICIT, true>(maxiter) : run_iterate<EULEREXPLICIT, false>(maxiter);
  }
  return MAXITER;
}

template <int INTEGRATOR, bool ABCFLAG>
int MinFireKokkos::run_iterate(int maxiter) {
  KK_FLOAT scale1 = 0.0, scale2 = 0.0; // Initialize to zero
  alpha_final = 0.0;
  int flagv0 = 1;

  atomKK->sync(Device, X_MASK | V_MASK | F_MASK | RMASS_MASK | TYPE_MASK);
  auto l_x = atomKK->k_x.view_device();
  auto l_v = atomKK->k_v.view_device();
  auto l_f = atomKK->k_f.view_device();
  auto l_rmass = atomKK->k_rmass.view_device();
  auto l_mass = atomKK->k_mass.view_device();
  auto l_type = atomKK->k_type.view_device();
  int nlocal = atom->nlocal;

  if constexpr (INTEGRATOR == LEAPFROG) {
    // ghost position/vel might change on host during legacy comm
    // but already just synced top of run_iterate()
    energy_force(0);
    neval++;
    KK_FLOAT dtf = -0.5 * dt * static_cast<KK_FLOAT>(force->ftm2v);
    Kokkos::parallel_for("min_fire/leapfrog_init", atom->nlocal, LAMMPS_LAMBDA(const int i) {
      const KK_FLOAT dtfm = dtf / (l_rmass.data() ? l_rmass(i) : l_mass(l_type(i)));
      l_v(i,0) = dtfm * l_f(i,0);
      l_v(i,1) = dtfm * l_f(i,1);
      l_v(i,2) = dtfm * l_f(i,2);
    });
  }

  for (int iter = 0; iter < maxiter; iter++) {
    if (timer->check_timeout(niter)) return TIMEOUT;
    bigint ntimestep = ++update->ntimestep;
    niter++;

    KK_ACC_FLOAT vdotf_local, vdotf_all;
    Kokkos::parallel_reduce("min_fire/vdotf", nlocal, LAMMPS_LAMBDA(const int i, KK_ACC_FLOAT &l_vdf) {
      // vdf += v0*f0 + v1*f1 + v2*f2
      l_vdf = Kokkos::fma(l_v(i,0), l_f(i,0),
              Kokkos::fma(l_v(i,1), l_f(i,1),
              Kokkos::fma(l_v(i,2), l_f(i,2), l_vdf)));
    }, vdotf_local);
    // bugfix for multiprocess replicas
    MPI_Allreduce(&vdotf_local, &vdotf_all, 1, MPI_KK_ACC_FLOAT, MPI_SUM, mpi_comm);

    if (vdotf_all > 0.0) {
      KK_ACC_FLOAT vdotv_local, fdotf_local, vdotv_all, fdotf_all;
      Kokkos::parallel_reduce("min_fire/norms", nlocal, LAMMPS_LAMBDA(const int i, KK_ACC_FLOAT &l_vv) {
        // vv += v0*v0 + v1*v1 + v2*v2
        l_vv = Kokkos::fma(l_v(i,0), l_v(i,0),
               Kokkos::fma(l_v(i,1), l_v(i,1),
               Kokkos::fma(l_v(i,2), l_v(i,2), l_vv)));
      }, vdotv_local);

      Kokkos::parallel_reduce("min_fire/fnorms", nlocal, LAMMPS_LAMBDA(const int i, KK_ACC_FLOAT &l_ff) {
        // ff += f0*f0 + f1*f1 + f2*f2
        l_ff = Kokkos::fma(l_f(i,0), l_f(i,0),
               Kokkos::fma(l_f(i,1), l_f(i,1),
               Kokkos::fma(l_f(i,2), l_f(i,2), l_ff)));
      }, fdotf_local);

      // bugfix for multiprocess replicas
      MPI_Allreduce(&vdotv_local, &vdotv_all, 1, MPI_KK_ACC_FLOAT, MPI_SUM, mpi_comm);
      MPI_Allreduce(&fdotf_local, &fdotf_all, 1, MPI_KK_ACC_FLOAT, MPI_SUM, mpi_comm);

      if constexpr (ABCFLAG) {
        if (alpha < 1e-10) alpha = 1e-10;
        const KK_FLOAT abc = (1.0 - pow(1.0 - alpha, (KK_FLOAT)(ntimestep - last_negative)));
        scale1 = (1.0 - alpha) / abc;
        scale2 = (fdotf_all <= 1e-20) ? 0.0 : (alpha * sqrt(vdotv_all / fdotf_all)) / abc;
      } else {
        scale1 = 1.0 - alpha;
        scale2 = (fdotf_all <= 1e-20) ? 0.0 : alpha * sqrt(vdotv_all / fdotf_all);
      }

      if (ntimestep - last_negative > delaystep) {
        dt = fmin(dt * dtgrow, dtmax);
        update->dt = dt;
        alpha *= alphashrink;
      }
      vdotf_negatif = 0;
    } else {
      last_negative = ntimestep;
      if (ntimestep - ntimestep_start >= delaystep || !delaystep_start_flag) {
        alpha = alpha0;
        if (dt * dtshrink >= dtmin) {
          dt *= dtshrink;
          update->dt = dt;
        }
      }

      vdotf_negatif++;
      if (max_vdotf_negatif > 0 && vdotf_negatif > max_vdotf_negatif) return MAXVDOTF;

      auto l_dt = dt;
      auto l_halfstepback_flag = halfstepback_flag;

      Kokkos::parallel_for("min_fire/inertia_reset", nlocal, LAMMPS_LAMBDA(const int i) {
        if (l_halfstepback_flag) {
          l_x(i,0) -= 0.5 * l_dt * l_v(i,0);
          l_x(i,1) -= 0.5 * l_dt * l_v(i,1);
          l_x(i,2) -= 0.5 * l_dt * l_v(i,2);
        }
        l_v(i,0) = l_v(i,1) = l_v(i,2) = 0.0;
      });
      flagv0 = 1;
    }

    if (!ABCFLAG && flagv0) {
      atomKK->modified(Device, X_MASK | V_MASK);
      atomKK->sync(Host, X_MASK | V_MASK);
      energy_force(0); // ghost position/vel might change on host during legacy comm
      atomKK->sync(Device, X_MASK | V_MASK);
      neval++;
      const KK_FLOAT dtf_init = dt * force->ftm2v;
      Kokkos::parallel_for("min_fire/v_init", nlocal, LAMMPS_LAMBDA(const int i) {
        const KK_FLOAT dtfm = dtf_init / (l_rmass.data() ? l_rmass(i) : l_mass(l_type(i)));
        l_v(i,0) = dtfm * l_f(i,0);
        l_v(i,1) = dtfm * l_f(i,1);
        l_v(i,2) = dtfm * l_f(i,2);
      });
    }

    KK_ACC_FLOAT dtvone = static_cast<KK_ACC_FLOAT>(dt);
    auto l_dmax = dmax;
    if constexpr (!ABCFLAG) {
      Kokkos::parallel_reduce("min_fire/dtv_limit", nlocal, LAMMPS_LAMBDA(const int i, KK_ACC_FLOAT &dtmin_local) {
         const KK_FLOAT vmax = fmax(fabs(l_v(i,0)), fmax(fabs(l_v(i,1)), fabs(l_v(i,2))));
        if (dtmin_local * vmax > l_dmax) dtmin_local = l_dmax / vmax;
      }, Kokkos::Min<KK_ACC_FLOAT>(dtvone));
      dtvone = Kokkos::min(dtvone, static_cast<KK_ACC_FLOAT>(dt));
    }
    KK_ACC_FLOAT dtv;
    // bugfix for multiprocess replicas
    MPI_Allreduce(&dtvone, &dtv, 1, MPI_KK_ACC_FLOAT, MPI_MIN, mpi_comm);

    if (flagv0) {
      Kokkos::parallel_for("min_fire/final_v_zero", nlocal, LAMMPS_LAMBDA(const int i) {
        l_v(i,0) = l_v(i,1) = l_v(i,2) = 0.0;
      });
    }

    const KK_FLOAT dtf_final = dtv * force->ftm2v;
    const KK_FLOAT dtf_half = 0.5 * dtf_final;
    Kokkos::parallel_for("min_fire/integrate", nlocal, LAMMPS_LAMBDA(const int i) {
      const KK_FLOAT mass_val = (l_rmass.data() ? l_rmass(i) : l_mass(l_type(i)));
      const KK_FLOAT dtfm = dtf_final / mass_val;
      if (INTEGRATOR == EULERIMPLICIT || INTEGRATOR == LEAPFROG) {
        l_v(i,0) += dtfm * l_f(i,0);
        l_v(i,1) += dtfm * l_f(i,1);
        l_v(i,2) += dtfm * l_f(i,2);
        if (vdotf_all > 0.0) {
          l_v(i,0) = scale1 * l_v(i,0) + scale2 * l_f(i,0);
          l_v(i,1) = scale1 * l_v(i,1) + scale2 * l_f(i,1);
          l_v(i,2) = scale1 * l_v(i,2) + scale2 * l_f(i,2);
          if (ABCFLAG) {
            // make sure that the displacement is not larger than dmax
            if (fabs(l_v(i,0)*dtv) > l_dmax) l_v(i,0) = l_dmax/dtv * l_v(i,0)/fabs(l_v(i,0));
            if (fabs(l_v(i,1)*dtv) > l_dmax) l_v(i,1) = l_dmax/dtv * l_v(i,1)/fabs(l_v(i,1));
            if (fabs(l_v(i,2)*dtv) > l_dmax) l_v(i,2) = l_dmax/dtv * l_v(i,2)/fabs(l_v(i,2));
          }
        }
        l_x(i,0) += dtv * l_v(i,0);
        l_x(i,1) += dtv * l_v(i,1);
        l_x(i,2) += dtv * l_v(i,2);
      } else if (INTEGRATOR == VERLET) {
        const KK_FLOAT dtfm_half = dtf_half / mass_val;
        l_v(i,0) += dtfm_half * l_f(i,0);
        l_v(i,1) += dtfm_half * l_f(i,1);
        l_v(i,2) += dtfm_half * l_f(i,2);
        if (vdotf_all > 0.0) {
          l_v(i,0) = scale1 * l_v(i,0) + scale2 * l_f(i,0);
          l_v(i,1) = scale1 * l_v(i,1) + scale2 * l_f(i,1);
          l_v(i,2) = scale1 * l_v(i,2) + scale2 * l_f(i,2);
          if (ABCFLAG) {
            // make sure that the displacement is not larger than dmax
            if (fabs(l_v(i,0)*dtv) > l_dmax) l_v(i,0) = l_dmax/dtv * l_v(i,0)/fabs(l_v(i,0));
            if (fabs(l_v(i,1)*dtv) > l_dmax) l_v(i,1) = l_dmax/dtv * l_v(i,1)/fabs(l_v(i,1));
            if (fabs(l_v(i,2)*dtv) > l_dmax) l_v(i,2) = l_dmax/dtv * l_v(i,2)/fabs(l_v(i,2));
          }
        }
        l_x(i,0) += dtv * l_v(i,0);
        l_x(i,1) += dtv * l_v(i,1);
        l_x(i,2) += dtv * l_v(i,2);
      } else if (INTEGRATOR == EULEREXPLICIT) {
        if (vdotf_all > 0.0) {
          l_v(i,0) = scale1 * l_v(i,0) + scale2 * l_f(i,0);
          l_v(i,1) = scale1 * l_v(i,1) + scale2 * l_f(i,1);
          l_v(i,2) = scale1 * l_v(i,2) + scale2 * l_f(i,2);
          if (ABCFLAG) {
            // make sure that the displacement is not larger than dmax
            if (fabs(l_v(i,0)*dtv) > l_dmax) l_v(i,0) = l_dmax/dtv * l_v(i,0)/fabs(l_v(i,0));
            if (fabs(l_v(i,1)*dtv) > l_dmax) l_v(i,1) = l_dmax/dtv * l_v(i,1)/fabs(l_v(i,1));
            if (fabs(l_v(i,2)*dtv) > l_dmax) l_v(i,2) = l_dmax/dtv * l_v(i,2)/fabs(l_v(i,2));
          }
        }
        l_x(i,0) += dtv * l_v(i,0);
        l_x(i,1) += dtv * l_v(i,1);
        l_x(i,2) += dtv * l_v(i,2);
        l_v(i,0) += dtfm * l_f(i,0);
        l_v(i,1) += dtfm * l_f(i,1);
        l_v(i,2) += dtfm * l_f(i,2);
      }
    });

    eprevious = ecurrent;
    atomKK->modified(Device, X_MASK | V_MASK);
    atomKK->sync(Host, X_MASK | V_MASK);
    ecurrent = energy_force(0); // ghost position/vel might change on host during legacy comm

    // energy_force() may migrate atoms / rebuild atom storage, so any cached
    // loop bounds and device views can become stale and must be refreshed.
    atomKK->sync(Device, X_MASK | V_MASK | F_MASK | TYPE_MASK | RMASS_MASK);
    nlocal = atom->nlocal;
    l_x = atomKK->k_x.view<DeviceType>();
    l_v = atomKK->k_v.view<DeviceType>();
    l_f = atomKK->k_f.view<DeviceType>();
    l_type = atomKK->k_type.view<DeviceType>();
    l_rmass = atomKK->k_rmass.view<DeviceType>();
    neval++;

    if constexpr (INTEGRATOR == VERLET) {
      Kokkos::parallel_for("min_fire/verlet_v_final", nlocal, LAMMPS_LAMBDA(const int i) {
        const KK_FLOAT dtfm_half = dtf_half / (l_rmass.data() ? l_rmass(i) : l_mass(l_type(i)));
        l_v(i,0) += dtfm_half * l_f(i,0);
        l_v(i,1) += dtfm_half * l_f(i,1);
        l_v(i,2) += dtfm_half * l_f(i,2);
      });
    }
    flagv0 = 0;

    // -------------------------------------------------
    // Corrected ETOL Check
    // -------------------------------------------------
    if (update->etol > 0.0 && ntimestep - last_negative > delaystep) {
      bool local_converged = (fabs(ecurrent - eprevious) <
          update->etol * 0.5 * (fabs(ecurrent) + fabs(eprevious) + EPS_ENERGY));

      if (update->multireplica == 0) {
        if (local_converged) return ETOL;
      } else {
        // MUST communicate regardless of local state
        int local_flag = local_converged ? 0 : 1;
        int global_flag;
        MPI_Allreduce(&local_flag, &global_flag, 1, MPI_INT, MPI_SUM, universe->uworld);
        // Only return if EVERYONE (sum=0) is converged
        if (global_flag == 0) return ETOL;
      }
    }

    // -------------------------------------------------
    // Corrected FTOL Check
    // -------------------------------------------------
    if (update->ftol > 0.0) {
      const KK_FLOAT fdotf = (normstyle == MAX) ? fnorm_max() : (normstyle == INF ? fnorm_inf() : fnorm_sqr());
      bool local_converged = (fdotf < update->ftol * update->ftol);

      if (update->multireplica == 0) {
        if (local_converged) return FTOL;
      } else {
        // MUST communicate regardless of local state
        int local_flag = local_converged ? 0 : 1;
        int global_flag;
        MPI_Allreduce(&local_flag, &global_flag, 1, MPI_INT, MPI_SUM, universe->uworld);
        if (global_flag == 0) return FTOL;
      }
    }

    if (output->next == ntimestep) {
      timer->stamp();
      output->write(ntimestep);
      timer->stamp(Timer::OUTPUT);
    }
  }
  atomKK->modified(Device, X_MASK | V_MASK);
  return MAXITER;
}
