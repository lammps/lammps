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

#include "fix_flow_gauss_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "error.h"
#include "group_kokkos.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixFlowGaussKokkos<DeviceType>::FixFlowGaussKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixFlowGauss(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = F_MASK | V_MASK | MASK_MASK | TYPE_MASK | RMASS_MASK;
  datamask_modify = F_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixFlowGaussKokkos<DeviceType>::init()
{
  FixFlowGauss::init();

  if (utils::strmatch(update->integrate_style,"^respa"))
    error->all(FLERR,"Cannot (yet) use respa with fix flow/gauss/kk");
}

/* ----------------------------------------------------------------------
   same as the base setup(), but the group mass is reduced on the device
   and the respa branch is unreachable (rejected in init())
------------------------------------------------------------------------- */

template<class DeviceType>
void FixFlowGaussKokkos<DeviceType>::setup(int vflag)
{
  // need to compute work done if fix_modify energy yes is set

  if (thermo_energy) workflag = true;

  // get total mass of group

  auto *groupKK = (GroupKokkos *) group;
  mTot = groupKK->mass_kk<DeviceType>(igroup);
  if (mTot <= 0.0)
    error->all(FLERR,"Invalid group mass in fix flow/gauss");

  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixFlowGaussKokkos<DeviceType>::post_force(int /*vflag*/)
{
  atomKK->sync(execution_space,datamask_read);

  d_f = atomKK->k_f.view<DeviceType>();
  d_v = atomKK->k_v.view<DeviceType>();
  d_mask = atomKK->k_mask.view<DeviceType>();
  d_type = atomKK->k_type.view<DeviceType>();
  d_rmass = atomKK->k_rmass.view<DeviceType>();
  d_mass = atomKK->k_mass.view<DeviceType>();
  atomKK->k_mass.template sync<DeviceType>();

  const int nlocal = atom->nlocal;

  // find the total force on all atoms in the group

  s_KK_double3 fsum;

  copymode = 1;
  Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType,TagFixFlowGaussReduce>(0,nlocal),
                          *this,fsum);
  copymode = 0;

  // add the processor sums together

  double f_thisProc[3] = {fsum.d0,fsum.d1,fsum.d2};
  MPI_Allreduce(f_thisProc,f_tot,3,MPI_DOUBLE,MPI_SUM,world);

  // a conserved direction contributes no applied acceleration

  for (int ii = 0; ii < 3; ii++)
    if (!flow[ii]) f_tot[ii] = 0.0;

  // compute applied acceleration

  for (int ii = 0; ii < 3; ii++) {
    a_app[ii] = -f_tot[ii] / mTot;
    m_a_app[ii] = static_cast<KK_FLOAT>(a_app[ii]);
  }

  // apply added acceleration to each atom
  // the added energy is more costly, so only accumulate it if requested

  copymode = 1;
  if (workflag) {
    double peAdded = 0.0;
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType,TagFixFlowGaussApplyWork>(0,nlocal),
                            *this,peAdded);
    double pe_tmp = 0.0;
    MPI_Allreduce(&peAdded,&pe_tmp,1,MPI_DOUBLE,MPI_SUM,world);
    pe_tot += pe_tmp;
  } else {
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagFixFlowGaussApply>(0,nlocal),*this);
  }
  copymode = 0;

  atomKK->modified(execution_space,datamask_modify);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixFlowGaussKokkos<DeviceType>::operator()(TagFixFlowGaussReduce, const int &i,
                                                s_KK_double3 &fsum) const
{
  if (d_mask[i] & groupbit) {
    fsum.d0 += static_cast<double>(d_f(i,0));
    fsum.d1 += static_cast<double>(d_f(i,1));
    fsum.d2 += static_cast<double>(d_f(i,2));
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixFlowGaussKokkos<DeviceType>::operator()(TagFixFlowGaussApply, const int &i) const
{
  if (d_mask[i] & groupbit) {
    KK_FLOAT f_app[3];
    applied_force(i,f_app);

    // f_app[jj] is 0 if flow[jj] is false

    d_f(i,0) += static_cast<KK_ACC_FLOAT>(f_app[0]);
    d_f(i,1) += static_cast<KK_ACC_FLOAT>(f_app[1]);
    d_f(i,2) += static_cast<KK_ACC_FLOAT>(f_app[2]);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixFlowGaussKokkos<DeviceType>::operator()(TagFixFlowGaussApplyWork, const int &i,
                                                double &peAdded) const
{
  if (d_mask[i] & groupbit) {
    KK_FLOAT f_app[3];
    applied_force(i,f_app);

    // f_app[jj] is 0 if flow[jj] is false

    d_f(i,0) += static_cast<KK_ACC_FLOAT>(f_app[0]);
    d_f(i,1) += static_cast<KK_ACC_FLOAT>(f_app[1]);
    d_f(i,2) += static_cast<KK_ACC_FLOAT>(f_app[2]);

    peAdded += static_cast<double>(f_app[0]*d_v(i,0) + f_app[1]*d_v(i,1) + f_app[2]*d_v(i,2));
  }
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixFlowGaussKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixFlowGaussKokkos<LMPHostType>;
#endif
}
