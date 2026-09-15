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

#include "compute_entropy_atom_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "domain.h"
#include "force.h"
#include "math_const.h"
#include "memory_kokkos.h"
#include "neigh_list_kokkos.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "pair.h"
#include "update.h"

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
ComputeEntropyAtomKokkos<DeviceType>::ComputeEntropyAtomKokkos(LAMMPS *lmp, int narg, char **arg) :
  ComputeEntropyAtom(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | MASK_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
ComputeEntropyAtomKokkos<DeviceType>::~ComputeEntropyAtomKokkos()
{
  if (copymode) return;

  memoryKK->destroy_kokkos(k_pair_entropy,pair_entropy);
  memoryKK->destroy_kokkos(k_pair_entropy_avg,pair_entropy_avg);
  pair_entropy = nullptr;
  pair_entropy_avg = nullptr;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ComputeEntropyAtomKokkos<DeviceType>::init()
{
  ComputeEntropyAtom::init();

  // adjust neighbor list request for KOKKOS

  auto request = neighbor->find_request(this);
  request->set_kokkos_host(std::is_same_v<DeviceType,LMPHostType> &&
                           !std::is_same_v<DeviceType,LMPDeviceType>);
  request->set_kokkos_device(std::is_same_v<DeviceType,LMPDeviceType>);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ComputeEntropyAtomKokkos<DeviceType>::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow pair_entropy and pair_entropy_avg arrays if necessary

  if (atom->nmax > nmax) {
    memoryKK->destroy_kokkos(k_pair_entropy,pair_entropy);
    memoryKK->destroy_kokkos(k_pair_entropy_avg,pair_entropy_avg);
    nmax = atom->nmax;
    memoryKK->create_kokkos(k_pair_entropy,pair_entropy,nmax,"entropy/atom:pair_entropy");
    d_pair_entropy = k_pair_entropy.template view<DeviceType>();
    if (!avg_flag) {
      vector_atom = pair_entropy;
    } else {
      memoryKK->create_kokkos(k_pair_entropy_avg,pair_entropy_avg,nmax,
                              "entropy/atom:pair_entropy_avg");
      d_pair_entropy_avg = k_pair_entropy_avg.template view<DeviceType>();
      vector_atom = pair_entropy_avg;
    }
  }

  // invoke occasional neighbor list build (if not perpetual)

  if (!avg_flag) neighbor->build_one(list);

  inum = list->inum + list->gnum;
  auto *k_list = static_cast<NeighListKokkos<DeviceType>*>(list);
  d_numneigh = k_list->d_numneigh;
  d_neighbors = k_list->d_neighbors;
  d_ilist = k_list->d_ilist;

  // the bin edges are the same for every atom

  if ((int)d_rbin.extent(0) != nbin) {
    d_rbin = typename AT::t_kkfloat_1d("entropy/atom:rbin",nbin);
    d_rbinsq = typename AT::t_kkfloat_1d("entropy/atom:rbinsq",nbin);
  }
  auto h_rbin = Kokkos::create_mirror_view(d_rbin);
  auto h_rbinsq = Kokkos::create_mirror_view(d_rbinsq);
  for (int i = 0; i < nbin; i++) {
    h_rbin(i) = static_cast<KK_FLOAT>(i*deltar);
    h_rbinsq(i) = h_rbin(i)*h_rbin(i);
  }
  Kokkos::deep_copy(d_rbin,h_rbin);
  Kokkos::deep_copy(d_rbinsq,h_rbinsq);

  // one g(r) histogram per atom, so that the kernel needs no scratch

  if ((int)d_gofr.extent(0) < inum || (int)d_gofr.extent(1) != nbin)
    d_gofr = typename AT::t_kkfloat_2d("entropy/atom:gofr",inum,nbin);

  const double volume = domain->xprd * domain->yprd * domain->zprd;

  sigmasq2_kk = static_cast<KK_FLOAT>(2.0*sigma*sigma);
  density_kk = static_cast<KK_FLOAT>(atom->natoms / volume);
  deltar_kk = static_cast<KK_FLOAT>(deltar);
  cutsq_kk = static_cast<KK_FLOAT>(cutsq);
  cutsq2_kk = static_cast<KK_FLOAT>(cutsq2);
  groupbit_kk = groupbit;

  const double neigh_cutoff = force->pair->cutforce + neighbor->skin;
  local_volume_kk = static_cast<KK_FLOAT>((4.0/3.0)*MY_PI*
                      neigh_cutoff*neigh_cutoff*neigh_cutoff);

  atomKK->sync(execution_space,datamask_read);
  x = atomKK->k_x.template view<DeviceType>();
  mask = atomKK->k_mask.template view<DeviceType>();

  copymode = 1;

  if (local_flag)
    Kokkos::parallel_for("ComputeEntropyAtom",
      Kokkos::RangePolicy<DeviceType, TagComputeEntropyAtom<1> >(0,inum),*this);
  else
    Kokkos::parallel_for("ComputeEntropyAtom",
      Kokkos::RangePolicy<DeviceType, TagComputeEntropyAtom<0> >(0,inum),*this);

  if (avg_flag)
    Kokkos::parallel_for("ComputeEntropyAtomAvg",
      Kokkos::RangePolicy<DeviceType, TagComputeEntropyAtomAvg>(0,inum),*this);

  copymode = 0;

  k_pair_entropy.template modify<DeviceType>();
  k_pair_entropy.sync_host();
  if (avg_flag) {
    k_pair_entropy_avg.template modify<DeviceType>();
    k_pair_entropy_avg.sync_host();
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int LOCAL>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void ComputeEntropyAtomKokkos<DeviceType>::operator()(TagComputeEntropyAtom<LOCAL>, const int &ii) const
{
  const int i = d_ilist[ii];
  if (!(mask[i] & groupbit_kk)) return;

  const KK_FLOAT xtmp = x(i,0);
  const KK_FLOAT ytmp = x(i,1);
  const KK_FLOAT ztmp = x(i,2);
  const int jnum = d_numneigh[i];

  KK_FLOAT density = density_kk;
  if (LOCAL) density = jnum / local_volume_kk;

  // kernel normalization: g(r) times the Gaussian

  const KK_FLOAT invNormConstantBase = static_cast<KK_FLOAT>(1.0) /
    (static_cast<KK_FLOAT>(4.0)*static_cast<KK_FLOAT>(MY_PI)*density *
     Kokkos::sqrt(static_cast<KK_FLOAT>(2.0)*static_cast<KK_FLOAT>(MY_PI))*static_cast<KK_FLOAT>(sigma));

  for (int k = 0; k < nbin; k++) d_gofr(ii,k) = static_cast<KK_FLOAT>(0.0);

  for (int jj = 0; jj < jnum; jj++) {
    int j = d_neighbors(i,jj);
    j &= NEIGHMASK;

    const KK_FLOAT delx = xtmp - x(j,0);
    const KK_FLOAT dely = ytmp - x(j,1);
    const KK_FLOAT delz = ztmp - x(j,2);
    const KK_FLOAT rsq = delx*delx + dely*dely + delz*delz;

    if (rsq < cutsq_kk) {
      const KK_FLOAT r = Kokkos::sqrt(rsq);
      const int bin = static_cast<int>(Kokkos::floor(r/deltar_kk));
      int minbin = bin - deltabin;
      if (minbin < 0) minbin = 0;
      if (minbin > nbin-1) minbin = nbin-1;
      int maxbin = bin + deltabin;
      if (maxbin > nbin-1) maxbin = nbin-1;
      for (int k = minbin; k < maxbin+1; k++) {
        const KK_FLOAT invNormKernel = invNormConstantBase/d_rbinsq[k];
        const KK_FLOAT distance = r - d_rbin[k];
        d_gofr(ii,k) += invNormKernel*Kokkos::exp(-distance*distance/sigmasq2_kk);
      }
    }
  }

  // integrand, then the trapezoid rule over the bins

  KK_FLOAT value = static_cast<KK_FLOAT>(0.0);
  for (int k = 0; k < nbin; k++) {
    const KK_FLOAT g = d_gofr(ii,k);
    KK_FLOAT integrand;
    if (g < static_cast<KK_FLOAT>(1.e-10))
      integrand = d_rbinsq[k];
    else
      integrand = (g*Kokkos::log(g) - g + static_cast<KK_FLOAT>(1.0))*d_rbinsq[k];

    // the CPU style sums the interior bins and then adds 0.5*integrand[0] and
    // 0.5*integrand[nbin-1] as two separate terms, so with a single bin that
    // bin is counted twice at half weight

    KK_FLOAT weight = static_cast<KK_FLOAT>(0.0);
    if ((k > 0) && (k < nbin-1)) weight += static_cast<KK_FLOAT>(1.0);
    if (k == 0) weight += static_cast<KK_FLOAT>(0.5);
    if (k == nbin-1) weight += static_cast<KK_FLOAT>(0.5);
    value += weight*integrand;
  }
  value *= deltar_kk;

  d_pair_entropy[i] = static_cast<KK_FLOAT>(-2.0)*static_cast<KK_FLOAT>(MY_PI)*density*value;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void ComputeEntropyAtomKokkos<DeviceType>::operator()(TagComputeEntropyAtomAvg, const int &ii) const
{
  const int i = d_ilist[ii];
  if (!(mask[i] & groupbit_kk)) return;

  const KK_FLOAT xtmp = x(i,0);
  const KK_FLOAT ytmp = x(i,1);
  const KK_FLOAT ztmp = x(i,2);
  const int jnum = d_numneigh[i];

  KK_FLOAT sum = d_pair_entropy[i];
  KK_FLOAT counter = static_cast<KK_FLOAT>(1.0);

  for (int jj = 0; jj < jnum; jj++) {
    int j = d_neighbors(i,jj);
    j &= NEIGHMASK;

    const KK_FLOAT delx = xtmp - x(j,0);
    const KK_FLOAT dely = ytmp - x(j,1);
    const KK_FLOAT delz = ztmp - x(j,2);
    const KK_FLOAT rsq = delx*delx + dely*dely + delz*delz;

    if (rsq < cutsq2_kk) {
      sum += d_pair_entropy[j];
      counter += static_cast<KK_FLOAT>(1.0);
    }
  }

  d_pair_entropy_avg[i] = sum/counter;
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class ComputeEntropyAtomKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class ComputeEntropyAtomKokkos<LMPHostType>;
#endif
}
