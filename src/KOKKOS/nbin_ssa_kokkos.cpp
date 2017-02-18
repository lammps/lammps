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

/* ----------------------------------------------------------------------
   Contributing authors:
   James Larentzos (ARL) and Timothy I. Mattox (Engility Corporation)
------------------------------------------------------------------------- */

#include "nbin_ssa_kokkos.h"
#include "neighbor.h"
#include "atom_kokkos.h"
#include "group.h"
#include "domain.h"
#include "comm.h"
#include "update.h"
#include "error.h"
#include "atom_masks.h"

// #include "memory.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
NBinSSAKokkos<DeviceType>::NBinSSAKokkos(LAMMPS *lmp) : NBinStandard(lmp)
{
  atoms_per_bin = ghosts_per_gbin = 16;

  d_resize = typename AT::t_int_scalar("NBinSSAKokkos::d_resize");
#ifndef KOKKOS_USE_CUDA_UVM
  h_resize = Kokkos::create_mirror_view(d_resize);
#else
  h_resize = d_resize;
#endif
  h_resize() = 1;

  k_gbincount = DAT::tdual_int_1d("NBinSSAKokkos::gbincount",8);
  gbincount = k_gbincount.view<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void NBinSSAKokkos<DeviceType>::bin_atoms_setup(int nall)
{
  if (mbins > (int) k_bins.d_view.dimension_0()) {
    k_bins = DAT::tdual_int_2d("NBinSSAKokkos::bins",mbins,atoms_per_bin);
    bins = k_bins.view<DeviceType>();

    k_bincount = DAT::tdual_int_1d("NBinSSAKokkos::bincount",mbins);
    bincount = k_bincount.view<DeviceType>();
  }

  ghosts_per_gbin = atom->nghost / 7; // estimate needed size

  if (ghosts_per_gbin > (int) k_gbins.d_view.dimension_1()) {
    k_gbins = DAT::tdual_int_2d("NBinSSAKokkos::gbins",8,ghosts_per_gbin);
    gbins = k_gbins.view<DeviceType>();
  }

  // Clear the local bin extent bounding box.
  h_lbinxlo() = mbinx - 1; // Safe to = stencil->sx + 1
  h_lbinylo() = mbiny - 1; // Safe to = stencil->sy + 1
  h_lbinzlo() = mbinz - 1; // Safe to = stencil->sz + 1
  h_lbinxhi() = 0; // Safe to = mbinx - stencil->sx - 1
  h_lbinyhi() = 0; // Safe to = mbiny - stencil->sy - 1
  h_lbinzhi() = 0; // Safe to = mbinz - stencil->sz - 1
  deep_copy(d_lbinxlo, h_lbinxlo);
  deep_copy(d_lbinylo, h_lbinylo);
  deep_copy(d_lbinzlo, h_lbinzlo);
  deep_copy(d_lbinxhi, h_lbinxhi);
  deep_copy(d_lbinyhi, h_lbinyhi);
  deep_copy(d_lbinzhi, h_lbinzhi);
}

/* ----------------------------------------------------------------------
   bin owned and ghost atoms for the Shardlow Splitting Algorithm (SSA)
   local atoms are in distinct bins (binhead[]) from the ghosts
   ghost atoms are "binned" in gairhead_ssa[] instead
     ghosts which are not in an Active Interaction Region (AIR) are skipped
------------------------------------------------------------------------- */

template<class DeviceType>
void NBinSSAKokkos<DeviceType>::bin_atoms()
{
  last_bin = update->ntimestep;

  int i;

  // bin the ghost atoms
  h_resize() = 1;
  while(h_resize() > 0) {
    h_resize() = 0;
    deep_copy(d_resize, h_resize);

    for (int i = 0; i < 8; i++) {
      k_gbincount.h_view(i) = 0;
    }
    k_gbincount.modify<LMPHostType>();
    k_gbincount.sync<DeviceType>();
    DeviceType::fence(); // FIXME?

    atomKK->sync(ExecutionSpaceFromDevice<DeviceType>::space,X_MASK);
    x = atomKK->k_x.view<DeviceType>();

    // I don't think these two lines need to be repeated here... - TIM 20170216
    sublo_[0] = domain->sublo[0];
    sublo_[1] = domain->sublo[1];
    sublo_[2] = domain->sublo[2];
    subhi_[0] = domain->subhi[0];
    subhi_[1] = domain->subhi[1];
    subhi_[2] = domain->subhi[2];

    NPairSSAKokkosBinGhostsFunctor<DeviceType> f(*this);

    Kokkos::parallel_for(atom->nghost, f);
    DeviceType::fence();

    deep_copy(h_resize, d_resize);
    if(h_resize()) {
      k_gbincount.modify<DeviceType>();
      k_gbincount.sync<DeviceType>();
      for (i = 1; i < 8; i++) {
        if (k_gbincount.h_view(i) > ghosts_per_gbin) {
          ghosts_per_gbin = k_gbincount.h_view(i);
        }
      }
      k_gbins = DAT::tdual_int_2d("gbins", 8, ghosts_per_gbin);
      gbins = k_gbins.view<DeviceType>();
    }
  }
  c_gbins = gbins; // gbins won't change until the next bin_atoms

  // bin the local atoms
  h_resize() = 1;
  while(h_resize() > 0) {
    h_resize() = 0;
    deep_copy(d_resize, h_resize);

    MemsetZeroFunctor<DeviceType> f_zero;
    f_zero.ptr = (void*) k_bincount.view<DeviceType>().ptr_on_device();
    Kokkos::parallel_for(mbins, f_zero);
    DeviceType::fence();

    atomKK->sync(ExecutionSpaceFromDevice<DeviceType>::space,X_MASK);
    x = atomKK->k_x.view<DeviceType>();

    // I don't think these two lines need to be repeated here... - TIM 20170216
    bboxlo_[0] = bboxlo[0]; bboxlo_[1] = bboxlo[1]; bboxlo_[2] = bboxlo[2];
    bboxhi_[0] = bboxhi[0]; bboxhi_[1] = bboxhi[1]; bboxhi_[2] = bboxhi[2];

    NPairSSAKokkosBinAtomsFunctor<DeviceType> f(*this);

    Kokkos::parallel_for(atom->nlocal, f);
    DeviceType::fence();

    deep_copy(h_resize, d_resize);
    if(h_resize()) {

      atoms_per_bin += 16;
      k_bins = DAT::tdual_int_2d("bins", mbins, atoms_per_bin);
      bins = k_bins.view<DeviceType>();
    }
  }
  deep_copy(h_lbinxlo, d_lbinxlo);
  deep_copy(h_lbinylo, d_lbinylo);
  deep_copy(h_lbinzlo, d_lbinzlo);
  deep_copy(h_lbinxhi, d_lbinxhi);
  deep_copy(h_lbinyhi, d_lbinyhi);
  deep_copy(h_lbinzhi, d_lbinzhi);
  c_bins = bins; // bins won't change until the next bin_atoms
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void NBinSSAKokkos<DeviceType>::binGhostsItem(const int &i_) const
{
  const int i = i_ + atom->nlocal;
  const int iAIR = coord2ssaAIR(x(i, 0), x(i, 1), x(i, 2));
  if (iAIR > 0) { // include only ghost atoms in an AIR
    const int ac = Kokkos::atomic_fetch_add(&gbincount[iAIR], (int)1);
    if(ac < (int) gbins.dimension_1()) {
      gbins(iAIR, ac) = i;
    } else {
      d_resize() = 1;
    }
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void NBinSSAKokkos<DeviceType>::binAtomsItem(const int &i) const
{
  int loc[3];
  const int ibin = coord2bin(x(i, 0), x(i, 1), x(i, 2), &(loc[0]));

  // Find the bounding box of the local atoms in the bins
  if (loc[0] < d_lbinxlo()) Kokkos::atomic_fetch_min(&d_lbinxlo(),loc[0]);
  if (loc[0] >= d_lbinxhi()) Kokkos::atomic_fetch_max(&d_lbinxhi(),loc[0] + 1);
  if (loc[1] < d_lbinylo()) Kokkos::atomic_fetch_min(&d_lbinylo(),loc[1]);
  if (loc[1] >= d_lbinyhi()) Kokkos::atomic_fetch_max(&d_lbinyhi(),loc[1] + 1);
  if (loc[2] < d_lbinzlo()) Kokkos::atomic_fetch_min(&d_lbinzlo(),loc[2]);
  if (loc[2] >= d_lbinzhi()) Kokkos::atomic_fetch_max(&d_lbinzhi(),loc[2] + 1);

  const int ac = Kokkos::atomic_fetch_add(&(bincount[ibin]), (int)1);
  if(ac < (int) bins.dimension_1()) {
    bins(ibin, ac) = i;
  } else {
    d_resize() = 1;
  }
}

namespace LAMMPS_NS {
template class NBinSSAKokkos<LMPDeviceType>;
#ifdef KOKKOS_HAVE_CUDA
template class NBinSSAKokkos<LMPHostType>;
#endif
}
