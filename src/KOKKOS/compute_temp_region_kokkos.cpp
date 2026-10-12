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

#include "compute_temp_region_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "kokkos_base.h"
#include "region.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
ComputeTempRegionKokkos<DeviceType>::ComputeTempRegionKokkos(LAMMPS *lmp, int narg, char **arg) :
  ComputeTempRegion(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  datamask_read = X_MASK | V_MASK | MASK_MASK | RMASS_MASK | TYPE_MASK;
  datamask_modify = EMPTY_MASK;

  maxbias = 0;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ComputeTempRegionKokkos<DeviceType>::init()
{
  ComputeTempRegion::init();

  if (!dynamic_cast<KokkosBase*>(region))
    error->all(FLERR, "Cannot use compute temp/region/kk with region style {} that has no KOKKOS support",
               region->style);
}

/* ----------------------------------------------------------------------
   evaluate the region for all atoms in the group into d_match
------------------------------------------------------------------------- */

template<class DeviceType>
void ComputeTempRegionKokkos<DeviceType>::region_match_all()
{
  region->prematch();

  const int nlocal = atom->nlocal;
  if ((int) k_match.extent(0) < nlocal)
    k_match = DAT::tdual_int_1d("temp/region:k_match", atom->nmax);

  KokkosBase *regionKKBase = dynamic_cast<KokkosBase*>(region);
  regionKKBase->match_all_kokkos(groupbit, k_match);
  k_match.template sync<DeviceType>();
  d_match = k_match.template view<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
double ComputeTempRegionKokkos<DeviceType>::compute_scalar()
{
  invoked_scalar = update->ntimestep;

  region_match_all();

  atomKK->sync(execution_space,datamask_read);
  atomKK->k_mass.sync<DeviceType>();

  v = atomKK->k_v.view<DeviceType>();
  if (atomKK->rmass)
    rmass = atomKK->k_rmass.view<DeviceType>();
  else
    mass = atomKK->k_mass.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();
  int nlocal = atom->nlocal;

  CTEMP t_kk;

  copymode = 1;
  if (atomKK->rmass)
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagComputeTempRegionScalar<1> >(0,nlocal),*this,t_kk);
  else
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagComputeTempRegionScalar<0> >(0,nlocal),*this,t_kk);
  copymode = 0;

  double tarray[2], tarray_all[2];
  tarray[0] = t_kk.t0;
  tarray[1] = t_kk.t1;
  MPI_Allreduce(tarray, tarray_all, 2, MPI_DOUBLE, MPI_SUM, world);
  dof = domain->dimension * tarray_all[0] - extra_dof;
  if (dof < 0.0 && tarray_all[0] > 0.0)
    error->all(FLERR, "Temperature compute degrees of freedom < 0");
  if (dof > 0)
    scalar = force->mvv2e * tarray_all[1] / (dof * force->boltz);
  else
    scalar = 0.0;
  return scalar;
}

template<class DeviceType>
template<int RMASS>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void ComputeTempRegionKokkos<DeviceType>::operator()(TagComputeTempRegionScalar<RMASS>, const int &i, CTEMP& t_kk) const {
  if ((mask[i] & groupbit) && d_match[i]) {
    KK_FLOAT massone = 0.0;
    if (RMASS) massone = rmass[i];
    else massone = mass[type[i]];
    t_kk.t0 += 1.0;
    t_kk.t1 += static_cast<double>((v(i,0)*v(i,0) + v(i,1)*v(i,1) + v(i,2)*v(i,2)) * massone);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ComputeTempRegionKokkos<DeviceType>::compute_vector()
{
  int i;

  invoked_vector = update->ntimestep;

  region_match_all();

  atomKK->sync(execution_space,datamask_read);
  atomKK->k_mass.sync<DeviceType>();

  v = atomKK->k_v.view<DeviceType>();
  if (atomKK->rmass)
    rmass = atomKK->k_rmass.view<DeviceType>();
  else
    mass = atomKK->k_mass.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();
  int nlocal = atom->nlocal;

  double t[6];
  for (i = 0; i < 6; i++) t[i] = 0.0;
  CTEMP t_kk;

  copymode = 1;
  if (atomKK->rmass)
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagComputeTempRegionVector<1> >(0,nlocal),*this,t_kk);
  else
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagComputeTempRegionVector<0> >(0,nlocal),*this,t_kk);
  copymode = 0;

  t[0] = t_kk.t0;
  t[1] = t_kk.t1;
  t[2] = t_kk.t2;
  t[3] = t_kk.t3;
  t[4] = t_kk.t4;
  t[5] = t_kk.t5;

  MPI_Allreduce(t,vector,6,MPI_DOUBLE,MPI_SUM,world);
  for (i = 0; i < 6; i++) vector[i] *= force->mvv2e;
}

template<class DeviceType>
template<int RMASS>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void ComputeTempRegionKokkos<DeviceType>::operator()(TagComputeTempRegionVector<RMASS>, const int &i, CTEMP& t_kk) const {
  if ((mask[i] & groupbit) && d_match[i]) {
    KK_FLOAT massone = 0.0;
    if (RMASS) massone = rmass[i];
    else massone = mass[type[i]];
    t_kk.t0 += static_cast<double>(massone * v(i,0)*v(i,0));
    t_kk.t1 += static_cast<double>(massone * v(i,1)*v(i,1));
    t_kk.t2 += static_cast<double>(massone * v(i,2)*v(i,2));
    t_kk.t3 += static_cast<double>(massone * v(i,0)*v(i,1));
    t_kk.t4 += static_cast<double>(massone * v(i,0)*v(i,2));
    t_kk.t5 += static_cast<double>(massone * v(i,1)*v(i,2));
  }
}

/* ----------------------------------------------------------------------
   remove velocity bias from all atoms to leave thermal velocity
   the bias is the full velocity of atoms outside the region
------------------------------------------------------------------------- */

template<class DeviceType>
void ComputeTempRegionKokkos<DeviceType>::remove_bias_all()
{
  remove_bias_all_kk();
  atomKK->sync(Host,V_MASK);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ComputeTempRegionKokkos<DeviceType>::remove_bias_all_kk()
{
  region_match_all();

  atomKK->sync(execution_space,V_MASK|MASK_MASK);
  v = atomKK->k_v.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();
  int nlocal = atom->nlocal;

  if (atom->nmax > maxbias) {
    maxbias = atom->nmax;
    vbiasall = typename AT::t_kkfloat_1d_3("temp/region/kk:vbiasall", maxbias);
  }

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeTempRegionRemoveBias >(0,nlocal),*this);
  copymode = 0;

  atomKK->modified(execution_space,V_MASK);
}

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void ComputeTempRegionKokkos<DeviceType>::operator()(TagComputeTempRegionRemoveBias, const int &i) const {
  if (mask[i] & groupbit) {
    if (d_match[i]) {
      vbiasall(i,0) = vbiasall(i,1) = vbiasall(i,2) = 0.0;
    } else {
      vbiasall(i,0) = v(i,0);
      vbiasall(i,1) = v(i,1);
      vbiasall(i,2) = v(i,2);
      v(i,0) = v(i,1) = v(i,2) = 0.0;
    }
  }
}

/* ----------------------------------------------------------------------
   add back in velocity bias to all atoms removed by remove_bias_all()
   assume remove_bias_all() was previously called
------------------------------------------------------------------------- */

template<class DeviceType>
void ComputeTempRegionKokkos<DeviceType>::restore_bias_all()
{
  restore_bias_all_kk();
  atomKK->sync(Host,V_MASK);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ComputeTempRegionKokkos<DeviceType>::restore_bias_all_kk()
{
  atomKK->sync(execution_space,V_MASK|MASK_MASK);
  v = atomKK->k_v.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();
  int nlocal = atom->nlocal;

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeTempRegionRestoreBias >(0,nlocal),*this);
  copymode = 0;

  atomKK->modified(execution_space,V_MASK);
}

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void ComputeTempRegionKokkos<DeviceType>::operator()(TagComputeTempRegionRestoreBias, const int &i) const {
  if (mask[i] & groupbit) {
    v(i,0) += vbiasall(i,0);
    v(i,1) += vbiasall(i,1);
    v(i,2) += vbiasall(i,2);
  }
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class ComputeTempRegionKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class ComputeTempRegionKokkos<LMPHostType>;
#endif
}
