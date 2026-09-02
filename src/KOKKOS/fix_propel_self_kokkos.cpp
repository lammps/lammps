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

#include "fix_propel_self_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "domain_kokkos.h"
#include "error.h"
#include "memory_kokkos.h"
#include "update.h"

using namespace LAMMPS_NS;

enum { DIPOLE, VELOCITY, QUAT };

static constexpr double TOL = 1.0e-14;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixPropelSelfKokkos<DeviceType>::FixPropelSelfKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixPropelSelf(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | V_MASK | F_MASK | MU_MASK | IMAGE_MASK | MASK_MASK;
  datamask_modify = F_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixPropelSelfKokkos<DeviceType>::~FixPropelSelfKokkos()
{
  if (copymode) return;
  memoryKK->destroy_kokkos(k_vatom,vatom);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixPropelSelfKokkos<DeviceType>::init()
{
  FixPropelSelf::init();

  // the quat mode needs the ellipsoid bonus arrays, which
  // AtomVecEllipsoidKokkos does not expose to device kernels

  if (mode == QUAT)
    error->all(FLERR,"Cannot (yet) use fix propel/self/kk with option quat");

  if (utils::strmatch(update->integrate_style,"^respa"))
    error->all(FLERR,"Cannot (yet) use respa with fix propel/self/kk");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixPropelSelfKokkos<DeviceType>::post_force(int vflag)
{
  atomKK->sync(execution_space,datamask_read);

  // energy and virial setup

  if (vflag) v_init(vflag);
  else evflag = 0;

  if (vflag_atom) {
    memoryKK->destroy_kokkos(k_vatom,vatom);
    memoryKK->create_kokkos(k_vatom,vatom,maxvatom,"propel/self:vatom");
    d_vatom = k_vatom.template view<DeviceType>();
  }

  d_x = atomKK->k_x.view<DeviceType>();
  d_v = atomKK->k_v.view<DeviceType>();
  d_f = atomKK->k_f.view<DeviceType>();
  d_image = atomKK->k_image.view<DeviceType>();
  d_mask = atomKK->k_mask.view<DeviceType>();
  if (mode == DIPOLE) d_mu = atomKK->k_mu.view<DeviceType>();

  // domain data has to be copied by value: a kernel cannot dereference domain

  prd = Few<double,3>(domain->prd);
  h = Few<double,6>(domain->h);
  triclinic = domain->triclinic;

  const int nlocal = atom->nlocal;

  double result[6] = {0.0,0.0,0.0,0.0,0.0,0.0};

  copymode = 1;
  if (mode == DIPOLE)
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType,TagFixPropelSelfDipole>(0,nlocal),
                            *this,result);
  else
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType,TagFixPropelSelfVelocity>(0,nlocal),
                            *this,result);
  copymode = 0;

  atomKK->modified(execution_space,datamask_modify);

  if (vflag_global)
    for (int k = 0; k < 6; k++) virial[k] += result[k];

  if (vflag_atom) {
    k_vatom.template modify<DeviceType>();
    k_vatom.sync_host();
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixPropelSelfKokkos<DeviceType>::operator()(TagFixPropelSelfDipole, const int &i,
                                                 value_type result) const
{
  if (d_mask[i] & groupbit) {
    const KK_FLOAT fx = static_cast<KK_FLOAT>(magnitude) * d_mu(i,0);
    const KK_FLOAT fy = static_cast<KK_FLOAT>(magnitude) * d_mu(i,1);
    const KK_FLOAT fz = static_cast<KK_FLOAT>(magnitude) * d_mu(i,2);

    d_f(i,0) += static_cast<KK_ACC_FLOAT>(fx);
    d_f(i,1) += static_cast<KK_ACC_FLOAT>(fy);
    d_f(i,2) += static_cast<KK_ACC_FLOAT>(fz);

    if (evflag) tally(result,i,fx,fy,fz);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixPropelSelfKokkos<DeviceType>::operator()(TagFixPropelSelfVelocity, const int &i,
                                                 value_type result) const
{
  if (d_mask[i] & groupbit) {
    const KK_FLOAT nv2 = d_v(i,0)*d_v(i,0) + d_v(i,1)*d_v(i,1) + d_v(i,2)*d_v(i,2);
    KK_FLOAT fnorm = 0.0;

    // without this check fnorm blows up for a nearly stationary particle

    if (nv2 > static_cast<KK_FLOAT>(TOL))
      fnorm = static_cast<KK_FLOAT>(magnitude) / Kokkos::sqrt(nv2);

    const KK_FLOAT fx = fnorm * d_v(i,0);
    const KK_FLOAT fy = fnorm * d_v(i,1);
    const KK_FLOAT fz = fnorm * d_v(i,2);

    d_f(i,0) += static_cast<KK_ACC_FLOAT>(fx);
    d_f(i,1) += static_cast<KK_ACC_FLOAT>(fy);
    d_f(i,2) += static_cast<KK_ACC_FLOAT>(fz);

    if (evflag) tally(result,i,fx,fy,fz);
  }
}

/* ----------------------------------------------------------------------
   virial contribution of the active force, using unwrapped coordinates
------------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixPropelSelfKokkos<DeviceType>::tally(value_type result, const int &i,
                                            KK_FLOAT fx, KK_FLOAT fy, KK_FLOAT fz) const
{
  Few<double,3> x_i;
  x_i[0] = static_cast<double>(d_x(i,0));
  x_i[1] = static_cast<double>(d_x(i,1));
  x_i[2] = static_cast<double>(d_x(i,2));
  auto unwrap = DomainKokkos::unmap(prd,h,triclinic,x_i,d_image(i));

  KK_FLOAT vi[6];
  vi[0] = fx * static_cast<KK_FLOAT>(unwrap[0]);
  vi[1] = fy * static_cast<KK_FLOAT>(unwrap[1]);
  vi[2] = fz * static_cast<KK_FLOAT>(unwrap[2]);
  vi[3] = fx * static_cast<KK_FLOAT>(unwrap[1]);
  vi[4] = fx * static_cast<KK_FLOAT>(unwrap[2]);
  vi[5] = fy * static_cast<KK_FLOAT>(unwrap[2]);

  if (vflag_global)
    for (int k = 0; k < 6; k++) result[k] += static_cast<double>(vi[k]);

  if (vflag_atom)
    for (int k = 0; k < 6; k++)
      Kokkos::atomic_add(&(d_vatom(i,k)),static_cast<KK_ACC_FLOAT>(vi[k]));
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixPropelSelfKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixPropelSelfKokkos<LMPHostType>;
#endif
}
