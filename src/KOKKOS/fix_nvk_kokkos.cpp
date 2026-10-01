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

#include "fix_nvk_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "error.h"
#include "force.h"
#include "update.h"

#include <cmath>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixNVKKokkos<DeviceType>::FixNVKKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixNVK(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | V_MASK | F_MASK | MASK_MASK | TYPE_MASK | RMASS_MASK;
  datamask_modify = X_MASK | V_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixNVKKokkos<DeviceType>::init()
{
  FixNVK::init();

  if (utils::strmatch(update->integrate_style,"^respa"))
    error->all(FLERR,"Cannot (yet) use respa with fix nvk/kk");
}

/* ----------------------------------------------------------------------
   one half of the isokinetic step: reduce a and b (Minary 2003 eqs 4.12,
   4.13), form s and sdot on the host, then rescale v (and advance x in the
   first half) per eqs 4.15-4.17
------------------------------------------------------------------------- */

template<class DeviceType>
void FixNVKKokkos<DeviceType>::integrate(int xupdate)
{
  atomKK->sync(execution_space,datamask_read);

  d_x = atomKK->k_x.view<DeviceType>();
  d_v = atomKK->k_v.view<DeviceType>();
  d_f = atomKK->k_f.view<DeviceType>();
  d_rmass = atomKK->k_rmass.view<DeviceType>();
  d_mass = atomKK->k_mass.view<DeviceType>();
  d_type = atomKK->k_type.view<DeviceType>();
  d_mask = atomKK->k_mask.view<DeviceType>();
  atomKK->k_mass.template sync<DeviceType>();

  const int rmass_flag = (atomKK->rmass) ? 1 : 0;

  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  s_FixNVK_ab ab;

  copymode = 1;
  if (rmass_flag)
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType,TagFixNVKReduce<1>>(0,nlocal),*this,ab);
  else
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType,TagFixNVKReduce<0>>(0,nlocal),*this,ab);
  copymode = 0;

  double a, b;
  MPI_Allreduce(&ab.a,&a,1,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(&ab.b,&b,1,MPI_DOUBLE,MPI_SUM,world);
  a /= (2.0*K_target);                    // units of inverse time
  b /= (2.0*K_target * force->mvv2e);     // units of inverse time squared
  const double sqtb = sqrt(b);
  const double s = a/b * (cosh(dtf*sqtb) - 1.0) + sinh(dtf*sqtb) / sqtb;
  const double sdot = a/b * sqtb * sinh(dtf*sqtb) + cosh(dtf*sqtb);

  l_s = static_cast<KK_FLOAT>(s);
  l_sdot = static_cast<KK_FLOAT>(sdot);
  l_dtv = static_cast<KK_FLOAT>(dtv);
  l_ftm2v = static_cast<KK_FLOAT>(force->ftm2v);

  copymode = 1;
  if (rmass_flag) {
    if (xupdate)
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagFixNVKUpdate<1,1>>(0,nlocal),*this);
    else
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagFixNVKUpdate<1,0>>(0,nlocal),*this);
  } else {
    if (xupdate)
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagFixNVKUpdate<0,1>>(0,nlocal),*this);
    else
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagFixNVKUpdate<0,0>>(0,nlocal),*this);
  }
  copymode = 0;

  atomKK->modified(execution_space, xupdate ? (X_MASK | V_MASK) : V_MASK);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixNVKKokkos<DeviceType>::initial_integrate(int /*vflag*/)
{
  integrate(1);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixNVKKokkos<DeviceType>::final_integrate()
{
  integrate(0);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int RMASS>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixNVKKokkos<DeviceType>::operator()(TagFixNVKReduce<RMASS>, const int &i,
                                          s_FixNVK_ab &ab) const
{
  if (d_mask[i] & groupbit) {
    const double fx = static_cast<double>(d_f(i,0));
    const double fy = static_cast<double>(d_f(i,1));
    const double fz = static_cast<double>(d_f(i,2));
    const double massone = RMASS ? static_cast<double>(d_rmass(i))
                                 : static_cast<double>(d_mass(d_type(i)));
    ab.a += fx*static_cast<double>(d_v(i,0)) + fy*static_cast<double>(d_v(i,1)) +
            fz*static_cast<double>(d_v(i,2));
    ab.b += (fx*fx + fy*fy + fz*fz) / massone;
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int RMASS, int XUPDATE>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixNVKKokkos<DeviceType>::operator()(TagFixNVKUpdate<RMASS,XUPDATE>, const int &i) const
{
  if (d_mask[i] & groupbit) {
    const KK_FLOAT massone = RMASS ? d_rmass(i) : d_mass(d_type(i));
    const KK_FLOAT sm = l_s / massone;
    d_v(i,0) = (d_v(i,0) + static_cast<KK_FLOAT>(d_f(i,0)) * sm * l_ftm2v) / l_sdot;
    d_v(i,1) = (d_v(i,1) + static_cast<KK_FLOAT>(d_f(i,1)) * sm * l_ftm2v) / l_sdot;
    d_v(i,2) = (d_v(i,2) + static_cast<KK_FLOAT>(d_f(i,2)) * sm * l_ftm2v) / l_sdot;
    if (XUPDATE) {
      d_x(i,0) += l_dtv * d_v(i,0);
      d_x(i,1) += l_dtv * d_v(i,1);
      d_x(i,2) += l_dtv * d_v(i,2);
    }
  }
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixNVKKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixNVKKokkos<LMPHostType>;
#endif
}
