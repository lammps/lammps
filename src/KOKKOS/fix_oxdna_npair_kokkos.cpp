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
#include "force.h"
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

  k_pair_counts = DAT::tdual_int_1d("FixOxdnaNpair:pair_counts", 2);
  screened_max_atoms = 0;
  screened_max_neigh = 0;
  screened_pair_count = 0;
  screen_cut_max = 0.0;
  screen_cutsq = static_cast<KK_FLOAT>(4.0);
  special_skip[0] = special_skip[1] = special_skip[2] = special_skip[3] = 0;
  coax_list_requested = false;
  coax_active = 0;
  coax_pair_count = 0;
  coax_max_atoms = 0;
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

  // the screen keeps only pairs closer than screen_cut_max + skin, so ask for
  // a list with that cutoff: it is then trimmed from a longer list (e.g. that
  // of oxdna*/dh) and the screening passes loop over fewer neighbors.
  // screen_cut_max was registered by the pair styles in init_one(), which
  // runs before the fixes are initialized.
  if (screen_cut_max > 0.0) request->set_cutoff_fixed(screen_cut_max);

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
      last_allocate != neighbor->nbuild) {
     compute_neigh_screen_to_npair();
     last_allocate = neighbor->nbuild;
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
      last_allocate != neighbor->nbuild) {
     compute_neigh_screen_to_npair();
     last_allocate = neighbor->nbuild;
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixOxdnaNpairKokkos<DeviceType>::update_screen_cutsq()
{
  // Derive the COM screen cutoff from the cutoffs registered by the consuming
  // pair styles (hbond / xstk / coaxstk) in their init_one, then add the
  // neighbor skin. Since this screened list is rebuilt only when the neighbor
  // list rebuilds, the full skin is required to keep the filtered pair list
  // valid between rebuilds: two atoms can approach each other by up to one
  // skin distance before the next rebuild (same Verlet-list principle as the
  // base neighbor list itself).
  const KK_FLOAT base_screen_cut = static_cast<KK_FLOAT>((screen_cut_max > 0.0) ? screen_cut_max : 2.0);
  const KK_FLOAT screen_cut_with_skin = base_screen_cut + static_cast<KK_FLOAT>(neighbor->skin);
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
    d_numneigh_screened = k_numneigh_screened.template view<DeviceType>();
    d_screened_offsets = k_screened_offsets.template view<DeviceType>();
  }

  atomKK->sync(execution_space, datamask_read);
  x = atomKK->k_x.view<DeviceType>();

  // the list of strand-end pairs for coaxial stacking (see request_coax_list())
  // is built in the same passes as the screened list, in the same pair order

  coax_active = coax_list_requested ? 1 : 0;
  coax_pair_count = 0;
  if (coax_active) {
    if (atom->nmax > coax_max_atoms) {
      coax_max_atoms = atom->nmax;
      MemKK::realloc_kokkos(k_numneigh_coax, "FixOxdnaNpair:numneigh_coax", coax_max_atoms);
      MemKK::realloc_kokkos(k_coax_offsets, "FixOxdnaNpair:coax_offsets", coax_max_atoms + 1);
      d_numneigh_coax = k_numneigh_coax.template view<DeviceType>();
      d_coax_offsets = k_coax_offsets.template view<DeviceType>();
    }
    atomKK->sync(execution_space, CG_DNA_MASK);
    id3p = atomKK->k_id3p.template view<DeviceType>();
    id5p = atomKK->k_id5p.template view<DeviceType>();
  }

  // Pairs whose special-bond weight is zero (1-2 bonded partners, which stay in
  // the neighbor list for the bonded excluded volume) are skipped by every
  // consumer of the screened list, so leave them out of it.
  for (int m = 0; m < 4; m++) special_skip[m] = (force->special_lj[m] == 0.0) ? 1 : 0;

  // Pass 1 (count): "TagFixOxdnaNpairNeighScreen" loops over each atom a and its
  // raw neighbours, runs 'screen_pair_fast' (a cheap CoM distance bool) for each,
  // and records only the surviving counts per atom in d_numneigh_screened (and
  // d_numneigh_coax). No per-atom survivor list is stored - the fill pass below
  // re-screens instead, which avoids an nmax x max_neigh scratch matrix.
  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixOxdnaNpairNeighScreen>(0, anum), *this);
  copymode = 0;

  // Perhaps "_local" suffixes are a little deceiving - these are shallow copies and point
  // to the same data as the non "_local" views. They're just for use in the lambda below
  // to avoid "this->" captures which the compiler would not like.
  const auto d_alist_local = d_alist;
  const auto d_numneigh_screened_local = d_numneigh_screened;
  const auto d_screened_offsets_local = d_screened_offsets;
  const auto d_numneigh_coax_local = d_numneigh_coax;
  const auto d_coax_offsets_local = d_coax_offsets;
  const auto d_pair_counts_local = k_pair_counts.template view<DeviceType>();
  const int anum_local = anum;
  const int coax_local = coax_active;

  // Pass 2 (scan): prefix-sum the per-atom counts of both lists (in neighbor-list
  // order) into their offsets, giving the starting flat index for each atom.
  // E.g. counts [2,0,3] -> offsets [0,2,2,5]; offsets(anum) is the total, which
  // is also stored in d_pair_counts for a single readback of both totals.
  // The Kokkos docs explain parallel_scan / prefix sum / "update" / "final".
  Kokkos::parallel_scan(
    Kokkos::RangePolicy<DeviceType>(0, anum + 1),
    KOKKOS_LAMBDA(const int i, FixOxdnaNpairScanCounts &update, const bool final) {
      if (i < anum_local) {
        const int a = d_alist_local(i);
        if (final) {
          d_screened_offsets_local(i) = update.screened;
          if (coax_local) d_coax_offsets_local(i) = update.coax;
        }
        update.screened += d_numneigh_screened_local(a);
        if (coax_local) update.coax += d_numneigh_coax_local(a);
      } else if (final) {
        d_screened_offsets_local(anum_local) = update.screened;
        if (coax_local) d_coax_offsets_local(anum_local) = update.coax;
        d_pair_counts_local(0) = update.screened;
        d_pair_counts_local(1) = update.coax;
      }
    });

  // one readback of both totals, needed on the host to size the lists and the
  // kernels that run over them
  k_pair_counts.template modify<DeviceType>();
  k_pair_counts.sync_host();
  screened_pair_count = k_pair_counts.view_host()(0);
  if (coax_active) coax_pair_count = k_pair_counts.view_host()(1);

  // size the packed pair lists by the number of pairs that survived screening,
  // with some headroom to avoid reallocating at every rebuild

  if ((bigint) screened_pair_count > (bigint) k_pairs_screened.extent(0)) {
    const bigint newsize = (bigint) screened_pair_count + screened_pair_count / 5 + 1;
    MemKK::realloc_kokkos(k_pairs_screened, "FixOxdnaNpair:pairs_screened", (size_t) newsize);
  }
  d_pairs_screened = k_pairs_screened.template view<DeviceType>();
  if (coax_active) {
    if ((bigint) coax_pair_count > (bigint) k_pairs_coax.extent(0)) {
      const bigint newsize = (bigint) coax_pair_count + coax_pair_count / 5 + 1;
      MemKK::realloc_kokkos(k_pairs_coax, "FixOxdnaNpair:pairs_coax", (size_t) newsize);
    }
    d_pairs_coax = k_pairs_coax.template view<DeviceType>();
  }

  // Pass 3 (fill): re-screen each atom's neighbours and write its survivors as
  // packed (a,b) uint64 keys directly at d_screened_offsets(i)..+count (and the
  // strand-end pairs among them at d_coax_offsets(i)..). The ComputeGPUPair
  // functors then run one thread per flat pair index, unpacking a (upper 32
  // bits) and b (lower 32 bits, special-bond bits preserved) with a single
  // global load.
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
  if (special_skip[braw >> SBBITS & 3]) return false;
  const int b = braw & NEIGHMASK;

  const KK_FLOAT b_com0 = x(b,0);
  const KK_FLOAT b_com1 = x(b,1);
  const KK_FLOAT b_com2 = x(b,2);

  KK_FLOAT delr_com[3];
  delr_com[0] = a_com0 - b_com0;
  delr_com[1] = a_com1 - b_com1;
  delr_com[2] = a_com2 - b_com2;

  // fma is fused-multipy-add op
  const KK_FLOAT rsq_com = Kokkos::fma(delr_com[2], delr_com[2],
                           Kokkos::fma(delr_com[1], delr_com[1], delr_com[0] * delr_com[0]));

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

  // coaxial stacking acts only between two strand ends
  const bool a_end = coax_active && is_strand_end(a);

  int nscreen = 0;
  int ncoax = 0;
  for (int ib = 0; ib < bnum; ib++) {
    const int braw = d_neighbors(a,ib);
    if (screen_pair_fast(braw, a_com0, a_com1, a_com2)) {
      nscreen++;
      if (a_end && is_strand_end(braw & NEIGHMASK)) ncoax++;
    }
  }
  d_numneigh_screened(a) = nscreen;
  if (coax_active) d_numneigh_coax(a) = ncoax;
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
  const bool a_end = coax_active && is_strand_end(a);
  int nscreen = d_screened_offsets(ia);
  int ncoax = a_end ? d_coax_offsets(ia) : 0;
  for (int ib = 0; ib < bnum; ib++) {
    const int braw = d_neighbors(a,ib);
    if (screen_pair_fast(braw, a_com0, a_com1, a_com2)) {
      const uint64_t key =
        (static_cast<uint64_t>(static_cast<uint32_t>(a)) << 32) |
        static_cast<uint64_t>(static_cast<uint32_t>(braw));
      d_pairs_screened(nscreen++) = key;
      if (a_end && is_strand_end(braw & NEIGHMASK)) d_pairs_coax(ncoax++) = key;
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
