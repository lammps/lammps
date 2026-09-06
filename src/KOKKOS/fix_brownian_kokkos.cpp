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

#include "fix_brownian_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "random_mars.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixBrownianKokkos<DeviceType>::FixBrownianKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixBrownian(lmp, narg, arg),
#ifndef LMP_KOKKOS_DEBUG_RNG
  rand_pool(seed + comm->me)
#else
  rand_pool(seed + comm->me, lmp)
#endif
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | V_MASK | F_MASK | MASK_MASK;
  datamask_modify = X_MASK | V_MASK;

#ifdef LMP_KOKKOS_DEBUG_RNG
  rand_pool.init(rng,seed + comm->me);
#endif
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixBrownianKokkos<DeviceType>::~FixBrownianKokkos()
{
  if (copymode) return;

#ifdef LMP_KOKKOS_DEBUG_RNG
  rand_pool.destroy();
#endif
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixBrownianKokkos<DeviceType>::init()
{
  FixBrownian::init();

  if (utils::strmatch(update->integrate_style,"^respa"))
    error->all(FLERR,"Cannot (yet) use respa with fix brownian/kk");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixBrownianKokkos<DeviceType>::initial_integrate(int /*vflag*/)
{
  if (domain->dimension == 2) {
    if (!noise_flag) {
      initial_integrate_kokkos<0,0,1>();
    } else if (gaussian_noise_flag) {
      initial_integrate_kokkos<0,1,1>();
    } else {
      initial_integrate_kokkos<1,0,1>();
    }
  } else {
    if (!noise_flag) {
      initial_integrate_kokkos<0,0,0>();
    } else if (gaussian_noise_flag) {
      initial_integrate_kokkos<0,1,0>();
    } else {
      initial_integrate_kokkos<1,0,0>();
    }
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int Tp_UNIFORM, int Tp_GAUSS, int Tp_2D>
void FixBrownianKokkos<DeviceType>::initial_integrate_kokkos()
{
  atomKK->sync(execution_space,datamask_read);

  d_x = atomKK->k_x.view<DeviceType>();
  d_v = atomKK->k_v.view<DeviceType>();
  d_f = atomKK->k_f.view<DeviceType>();
  d_mask = atomKK->k_mask.view<DeviceType>();

  l_dt = static_cast<KK_FLOAT>(dt);
  l_g1 = static_cast<KK_FLOAT>(g1);
  l_g2 = static_cast<KK_FLOAT>(g2);

  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  copymode = 1;
  Kokkos::parallel_for(
    Kokkos::RangePolicy<DeviceType,TagFixBrownian<Tp_UNIFORM,Tp_GAUSS,Tp_2D>>(0,nlocal),*this);
  copymode = 0;

  atomKK->modified(execution_space,datamask_modify);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int Tp_UNIFORM, int Tp_GAUSS, int Tp_2D>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixBrownianKokkos<DeviceType>::operator()(TagFixBrownian<Tp_UNIFORM,Tp_GAUSS,Tp_2D>,
                                               const int &i) const
{
  if (d_mask[i] & groupbit) {
    KK_FLOAT dx, dy, dz;

    rand_type rand_gen = rand_pool.get_state();

    if (Tp_2D) {
      dz = 0.0;
      if (Tp_UNIFORM) {
        dx = l_dt * (l_g1 * static_cast<KK_FLOAT>(d_f(i,0))
                     + l_g2 * static_cast<KK_FLOAT>(rand_gen.drand() - 0.5));
        dy = l_dt * (l_g1 * static_cast<KK_FLOAT>(d_f(i,1))
                     + l_g2 * static_cast<KK_FLOAT>(rand_gen.drand() - 0.5));
      } else if (Tp_GAUSS) {
        dx = l_dt * (l_g1 * static_cast<KK_FLOAT>(d_f(i,0))
                     + l_g2 * static_cast<KK_FLOAT>(rand_gen.normal()));
        dy = l_dt * (l_g1 * static_cast<KK_FLOAT>(d_f(i,1))
                     + l_g2 * static_cast<KK_FLOAT>(rand_gen.normal()));
      } else {
        dx = l_dt * l_g1 * static_cast<KK_FLOAT>(d_f(i,0));
        dy = l_dt * l_g1 * static_cast<KK_FLOAT>(d_f(i,1));
      }
    } else {
      if (Tp_UNIFORM) {
        dx = l_dt * (l_g1 * static_cast<KK_FLOAT>(d_f(i,0))
                     + l_g2 * static_cast<KK_FLOAT>(rand_gen.drand() - 0.5));
        dy = l_dt * (l_g1 * static_cast<KK_FLOAT>(d_f(i,1))
                     + l_g2 * static_cast<KK_FLOAT>(rand_gen.drand() - 0.5));
        dz = l_dt * (l_g1 * static_cast<KK_FLOAT>(d_f(i,2))
                     + l_g2 * static_cast<KK_FLOAT>(rand_gen.drand() - 0.5));
      } else if (Tp_GAUSS) {
        dx = l_dt * (l_g1 * static_cast<KK_FLOAT>(d_f(i,0))
                     + l_g2 * static_cast<KK_FLOAT>(rand_gen.normal()));
        dy = l_dt * (l_g1 * static_cast<KK_FLOAT>(d_f(i,1))
                     + l_g2 * static_cast<KK_FLOAT>(rand_gen.normal()));
        dz = l_dt * (l_g1 * static_cast<KK_FLOAT>(d_f(i,2))
                     + l_g2 * static_cast<KK_FLOAT>(rand_gen.normal()));
      } else {
        dx = l_dt * l_g1 * static_cast<KK_FLOAT>(d_f(i,0));
        dy = l_dt * l_g1 * static_cast<KK_FLOAT>(d_f(i,1));
        dz = l_dt * l_g1 * static_cast<KK_FLOAT>(d_f(i,2));
      }
    }

    rand_pool.free_state(rand_gen);

    d_x(i,0) += dx;
    d_v(i,0) = dx / l_dt;

    d_x(i,1) += dy;
    d_v(i,1) = dy / l_dt;

    d_x(i,2) += dz;
    d_v(i,2) = dz / l_dt;
  }
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixBrownianKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixBrownianKokkos<LMPHostType>;
#endif
}
