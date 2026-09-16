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

#include "fix_oxdna_npair_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "kokkos.h"
#include "memory_kokkos.h"
#include "neighbor.h"
#include "neigh_list_kokkos.h"
#include "neigh_request.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixOxdnaNpairKokkos<DeviceType>::FixOxdnaNpairKokkos(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  datamask_read = X_MASK | TYPE_MASK;
  datamask_modify = EMPTY_MASK;

  k_screened_pair_count = DAT::tdual_int_scalar("FixOxdnaNpair:screened_pair_count");
  screened_max_atoms = 0;
  screened_max_neigh = 0;
  screened_pair_count = 0;
  screen_cut_max = 0.0;
  screen_cutsq = static_cast<KK_FLOAT>(4.0);
  force_screening_all_backends = false;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixOxdnaNpairKokkos<DeviceType>::~FixOxdnaNpairKokkos() = default;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixOxdnaNpairKokkos<DeviceType>::init()
{
  // adjust neighbor list request for KOKKOS
  neighbor->add_request(this);
  neighflag = lmp->kokkos->neighflag;
  auto request = neighbor->find_request(this);
  request->set_kokkos_host(std::is_same_v<DeviceType,LMPHostType> &&
                           !std::is_same_v<DeviceType,LMPDeviceType>);
  request->set_kokkos_device(std::is_same_v<DeviceType,LMPDeviceType>);
  if (neighflag == FULL) request->enable_full();

  last_allocate = -1;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int FixOxdnaNpairKokkos<DeviceType>::setmask()
{
  int mask = 0;
  mask |= MIN_PRE_FORCE;
  mask |= PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixOxdnaNpairKokkos<DeviceType>::init_list(int, class NeighList* ptr)
{
  this->list = ptr;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixOxdnaNpairKokkos<DeviceType>::min_setup_pre_force(int vflag)
{
  update_screen_cutsq();
  min_pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixOxdnaNpairKokkos<DeviceType>::min_pre_force(int /*vflag*/)
{
  if ((force_screening_all_backends || execution_space != HostKK) &&
      last_allocate != neighbor->lastcall) {
     compute_neigh_screen_to_npair();
     last_allocate = neighbor->lastcall;
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixOxdnaNpairKokkos<DeviceType>::setup_pre_force(int vflag)
{
  update_screen_cutsq();
  pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixOxdnaNpairKokkos<DeviceType>::pre_force(int /*vflag*/)
{
  if ((force_screening_all_backends || execution_space != HostKK) &&
      last_allocate != neighbor->lastcall) {
     compute_neigh_screen_to_npair();
     last_allocate = neighbor->lastcall;
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixOxdnaNpairKokkos<DeviceType>::update_screen_cutsq()
{
  // Derive the COM screen cutoff from the cutoffs registered by the consuming
  // pair styles (hbond / xstk / coaxstk) in their init_one, then add the
  // neighbor skin. Since this screened list is rebuilt only when the neighbor
  // list rebuilds, a skin margin is required to keep the filtered pair list
  // valid between rebuilds (same Verlet-list principle as the base neighbor
  // list itself).
  // However, the pair_styles already add in an extra 0.4*nx margin
  // to their cutoffs to account for the base-site offset, so we can be
  // a little cheeky and half the neighbor skin margin too.
  const KK_FLOAT base_screen_cut = (screen_cut_max > 0.0) ? screen_cut_max : 2.0;
  const KK_FLOAT screen_cut_with_skin = base_screen_cut + (0.5 * neighbor->skin);
  screen_cutsq = screen_cut_with_skin * screen_cut_with_skin;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixOxdnaNpairKokkos<DeviceType>::compute_neigh_screen_to_npair()
{
  // get the neighbor list and neighbors used in operator()
  NeighListKokkos<DeviceType>* k_list = static_cast<NeighListKokkos<DeviceType>*>(this->list);
  d_neighbors = k_list->d_neighbors;
  anum = this->list->inum;
  d_alist = k_list->d_ilist;
  d_numneigh = k_list->d_numneigh;

  // reallocate screened neighbor arrays if necessary
  screened_pair_count = 0;
  const int max_atoms = atom->nmax;
  const int max_neigh = d_neighbors.extent(1);
  if (max_atoms > screened_max_atoms || max_neigh > screened_max_neigh) {
    screened_max_atoms = max_atoms;
    screened_max_neigh = max_neigh;
    MemKK::realloc_kokkos(k_numneigh_screened, "FixOxdnaNpair:numneigh_screened",
                          screened_max_atoms);
    MemKK::realloc_kokkos(k_screened_offsets, "FixOxdnaNpair:screened_offsets",
                          screened_max_atoms + 1);
    MemKK::realloc_kokkos(k_pairs_screened, "FixOxdnaNpair:pairs_screened",
              screened_max_atoms * screened_max_neigh);
    d_numneigh_screened = k_numneigh_screened.template view<DeviceType>();
    d_screened_offsets = k_screened_offsets.template view<DeviceType>();
    d_pairs_screened = k_pairs_screened.template view<DeviceType>();
  }

  atomKK->sync(execution_space, datamask_read);
  x = atomKK->k_x.view<DeviceType>();

  // Pass 1 (count): "TagFixOxdnaNpairNeighScreen" loops over each atom a and its
  // raw neighbours, runs 'screen_pair_fast' (a cheap CoM distance bool) for each,
  // and records only the surviving count per atom in d_numneigh_screened. No
  // per-atom survivor list is stored - the fill pass below re-screens instead,
  // which avoids an nmax x max_neigh scratch matrix.
  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixOxdnaNpairNeighScreen>(0, anum), *this);
  copymode = 0;

  // Perhaps "_local" suffixes are a little deceiving - these are shallow copies and point
  // to the same data as the non "_local" views. They're just for use in the lambda below
  // to avoid "this->" captures which the compiler would not like.
  const auto d_alist_local = d_alist;
  const auto d_numneigh_screened_local = d_numneigh_screened;
  const auto d_screened_offsets_local = d_screened_offsets;
  const int anum_local = anum;

  // Pass 2 (scan): prefix-sum the per-atom screened counts (in neighbor-list
  // order) into d_screened_offsets, giving the starting flat index for each atom.
  // E.g. counts [2,0,3] -> offsets [0,2,2,5]; offsets(anum) is the total. The
  // Kokkos docs explain parallel_scan / prefix sum / "update" / "final".
  Kokkos::parallel_scan(
    Kokkos::RangePolicy<DeviceType>(0, anum + 1),
    KOKKOS_LAMBDA(const int i, int &update, const bool final) {
      if (i < anum_local) {
        if (final) d_screened_offsets_local(i) = update;
        update += d_numneigh_screened_local(d_alist_local(i));
      } else if (final) {
        d_screened_offsets_local(anum_local) = update;
      }
    });

  // After the parallel_scan, the subview just gives us the value at the last offset,
  // i.e. the total screened pair count.
  // screened_pair_count is a host-side int, so we can't just directly read the value
  // from the device-side d_screened_offsets view (would trigger a seg-fault).
  Kokkos::deep_copy(
    k_screened_pair_count.view_host(), Kokkos::subview(d_screened_offsets_local, anum_local));
  screened_pair_count = k_screened_pair_count.view_host()();

  // Pass 3 (fill): re-screen each atom's neighbours and write its survivors as
  // packed (a,b) uint64 keys directly at d_screened_offsets(i)..+count. The
  // ComputeGPUPair functors then run one thread per flat pair index, unpacking
  // a (upper 32 bits) and b (lower 32 bits, special-bond bits preserved) with a
  // single global load.
  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixOxdnaNpairFill>(0, anum), *this);
  copymode = 0;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
bool FixOxdnaNpairKokkos<DeviceType>::screen_pair_fast(const int &braw,
                                                       const KK_FLOAT &a_com0,
                                                       const KK_FLOAT &a_com1,
                                                       const KK_FLOAT &a_com2) const
{
  const int b = braw & NEIGHMASK;

  const KK_FLOAT b_com0 = x(b,0);
  const KK_FLOAT b_com1 = x(b,1);
  const KK_FLOAT b_com2 = x(b,2);

  KK_FLOAT delr_com[3];
  delr_com[0] = a_com0 - b_com0;
  delr_com[1] = a_com1 - b_com1;
  delr_com[2] = a_com2 - b_com2;

  // fma is fused-multipy-add op
  const KK_FLOAT rsq_com = fma(delr_com[2], delr_com[2],
                           fma(delr_com[1], delr_com[1], delr_com[0] * delr_com[0]));

  // Boolean screen against the derived COM cutoff (set in
  // compute_neigh_screen_to_npair from the consuming styles' registered cutoffs).
  return (rsq_com < screen_cutsq);
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixOxdnaNpairKokkos<DeviceType>::operator()(TagFixOxdnaNpairNeighScreen, const int &ia) const
{
  const int a = d_alist(ia);
  const int bnum = d_numneigh(a);
  const KK_FLOAT a_com0 = x(a,0);
  const KK_FLOAT a_com1 = x(a,1);
  const KK_FLOAT a_com2 = x(a,2);

  int nscreen = 0;
  for (int ib = 0; ib < bnum; ib++) {
    const int braw = d_neighbors(a,ib);
    if (screen_pair_fast(braw, a_com0, a_com1, a_com2)) nscreen++;
  }
  d_numneigh_screened(a) = nscreen;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixOxdnaNpairKokkos<DeviceType>::operator()(TagFixOxdnaNpairFill, const int &ia) const
{
  const int a = d_alist(ia);
  const int bnum = d_numneigh(a);
  const KK_FLOAT a_com0 = x(a,0);
  const KK_FLOAT a_com1 = x(a,1);
  const KK_FLOAT a_com2 = x(a,2);

  // Re-screen with the same predicate used in pass 1; write survivors as packed
  // (a, braw) keys at the scanned base offset for this atom. braw keeps the
  // special-bond bits so ComputeGPUPair can apply special_lj/sbmask on unpack.
  int nscreen = d_screened_offsets(ia);
  for (int ib = 0; ib < bnum; ib++) {
    const int braw = d_neighbors(a,ib);
    if (screen_pair_fast(braw, a_com0, a_com1, a_com2)) {
      d_pairs_screened(nscreen++) =
        (static_cast<uint64_t>(static_cast<uint32_t>(a)) << 32) |
        static_cast<uint64_t>(static_cast<uint32_t>(braw));
    }
  }
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixOxdnaNpairKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixOxdnaNpairKokkos<LMPHostType>;
#endif
}    // namespace LAMMPS_NS
