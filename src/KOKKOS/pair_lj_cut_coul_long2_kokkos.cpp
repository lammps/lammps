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

#include "pair_lj_cut_coul_long2_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "ewald_const.h"
#include "force.h"
#include "kokkos.h"
#include "neigh_list_kokkos.h"

using namespace LAMMPS_NS;
using namespace EwaldConst;

/* ----------------------------------------------------------------------
   Hand-written force-only kernel for lj/cut/coul/long, HALFTHREAD + newton.

   Warp-per-atom over a HALF list: one warp's vector lanes split atom i's
   neighbor loop, the i-force is reduced across the lanes, and each pair
   scatters -f to j with an atomic.  This is the configuration the stock
   neigh/thread path cannot express (it requires newton off + a full list,
   doubling the pair count); a half list keeps the pair count halved while
   still getting the memory-level parallelism that hides the x[j] gather
   latency.  Coulomb uses the analytical Ewald correction (no table lookup).
   Parameters travel in the functor (constant memory on CUDA) via the stack
   arrays, as in the STACKPARAMS path of the generic kernel.
------------------------------------------------------------------------- */

namespace {

template<class DeviceType>
struct LJCL2Force {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef Kokkos::TeamPolicy<DeviceType, Kokkos::IndexType<int>> policy_type;
  typedef typename policy_type::member_type member_type;

  typename AT::t_kkfloat_1d_3_lr_randomread x;
  typename AT::t_kkacc_1d_3 f;
  typename AT::t_kkfloat_1d_randomread q;
  typename AT::t_int_1d_randomread type;
  typename AT::t_int_1d d_ilist;
  typename AT::t_int_1d d_numneigh;
  typename AT::t_neighbors_2d d_neighbors;        // LayoutLeft  (atom-fast)
  typename AT::t_neighbors_2d_lr d_neighbors_t;   // LayoutRight (jj-fast), transpose

  int nlocal, newton_pair, inum, atoms_per_team;
  bool full;   // full list (newton off): sum i-force only, no atomic j-scatter
  bool use_transpose;  // read neighbor indices from the LayoutRight transpose list.
                       // In warp-per-atom, the vector lanes read d_neighbors(i,jj)
                       // for consecutive jj; LayoutLeft strides those by nmax
                       // (fully uncoalesced), LayoutRight makes them contiguous.
  KK_FLOAT g_ewald, qqrd2e;
  KK_FLOAT special_lj[4], special_coul[4];
  params_lj_coul m_params[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
  KK_FLOAT m_cutsq[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
  KK_FLOAT m_cut_ljsq[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
  KK_FLOAT m_cut_coulsq[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];

  // Force contribution of one (i,jj) neighbor pair: always folds the i-force into
  // (fx,fy,fz).  On a HALF list (full==false) it also atomic-scatters -f to j.
  // On a FULL list (full==true) j supplies its own force from its own row, so the
  // scatter is skipped -- that removes every per-pair atomic, the throughput wall
  // on A100, at the cost of computing each pair twice.
  KOKKOS_INLINE_FUNCTION
  void pair_contrib(const int i, const KK_FLOAT xtmp, const KK_FLOAT ytmp,
                    const KK_FLOAT ztmp, const KK_FLOAT qtmp, const int itype,
                    const int jraw, KK_ACC_FLOAT& fx, KK_ACC_FLOAT& fy,
                    KK_ACC_FLOAT& fz) const {
    const int sb = (jraw >> SBBITS) & 3;
    const KK_FLOAT factor_lj = special_lj[sb];
    const KK_FLOAT factor_coul = special_coul[sb];
    const int j = jraw & NEIGHMASK;

    const KK_FLOAT delx = xtmp - x(j,0);
    const KK_FLOAT dely = ytmp - x(j,1);
    const KK_FLOAT delz = ztmp - x(j,2);
    const int jtype = type(j);
    const KK_FLOAT rsq = delx*delx + dely*dely + delz*delz;

    if (rsq < m_cutsq[itype][jtype]) {
      // one reciprocal-sqrt supplies rinv, r2inv, and r (was 2 divides + 1 sqrt)
      const KK_FLOAT rinv = Kokkos::rsqrt(rsq);
      const KK_FLOAT r2inv = rinv*rinv;
      KK_FLOAT fpair = static_cast<KK_FLOAT>(0.0);

      if (rsq < m_cut_ljsq[itype][jtype]) {
        const KK_FLOAT r6inv = r2inv*r2inv*r2inv;
        fpair += factor_lj *
          r6inv*(m_params[itype][jtype].lj1*r6inv - m_params[itype][jtype].lj2) * r2inv;
      }

      if (rsq < m_cut_coulsq[itype][jtype]) {
        // analytical Ewald real-space correction (Abramowitz-Stegun 7.1.26)
        const KK_FLOAT r = rsq*rinv;
        const KK_FLOAT grij = g_ewald * r;
        const KK_FLOAT expm2 = Kokkos::exp(-grij*grij);
        const KK_FLOAT t = static_cast<KK_FLOAT>(1.0) /
          (static_cast<KK_FLOAT>(1.0) + static_cast<KK_FLOAT>(EWALD_P)*grij);
        const KK_FLOAT erfc = t * (static_cast<KK_FLOAT>(A1)+t*(static_cast<KK_FLOAT>(A2)+
                              t * (static_cast<KK_FLOAT>(A3)+t*(static_cast<KK_FLOAT>(A4)+
                              t * static_cast<KK_FLOAT>(A5))))) * expm2;
        const KK_FLOAT prefactor = qqrd2e * qtmp*q(j) * rinv;
        KK_FLOAT forcecoul = prefactor * (erfc + static_cast<KK_FLOAT>(EWALD_F)*grij*expm2);
        if (factor_coul < static_cast<KK_FLOAT>(1.0))
          forcecoul -= (static_cast<KK_FLOAT>(1.0)-factor_coul)*prefactor;
        fpair += forcecoul * r2inv;
      }

      fx += static_cast<KK_ACC_FLOAT>(delx*fpair);
      fy += static_cast<KK_ACC_FLOAT>(dely*fpair);
      fz += static_cast<KK_ACC_FLOAT>(delz*fpair);

      if (!full && (newton_pair || j < nlocal)) {
        Kokkos::atomic_add(&f(j,0), -static_cast<KK_ACC_FLOAT>(delx*fpair));
        Kokkos::atomic_add(&f(j,1), -static_cast<KK_ACC_FLOAT>(dely*fpair));
        Kokkos::atomic_add(&f(j,2), -static_cast<KK_ACC_FLOAT>(delz*fpair));
      }
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const member_type& team) const {
    const int atom_in_team = team.league_rank()*atoms_per_team + team.team_rank();
    if (atom_in_team >= inum) return;

    const int i = d_ilist(atom_in_team);
    const KK_FLOAT xtmp = x(i,0);
    const KK_FLOAT ytmp = x(i,1);
    const KK_FLOAT ztmp = x(i,2);
    const KK_FLOAT qtmp = q(i);
    const int itype = type(i);
    const int jnum = d_numneigh(i);

    KK_ACC_FLOAT fxtmp = 0.0, fytmp = 0.0, fztmp = 0.0;

    // vector lanes of this atom's warp split the neighbor loop; the i-force
    // contributions are reduced across lanes into (fxtmp,fytmp,fztmp).  (A 2-way
    // ILP unroll here regressed the A100 ~30%: it is throughput/atomic-bound, not
    // latency-bound, so extra in-flight gathers only deepen the memory queues.)
    Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(team, jnum),
      [&] (const int jj, KK_ACC_FLOAT& fx, KK_ACC_FLOAT& fy, KK_ACC_FLOAT& fz) {
      // LayoutRight transpose makes these consecutive-jj reads coalesce across the
      // warp's vector lanes; LayoutLeft strides them by nmax (uncoalesced).
      const int jraw = use_transpose ? d_neighbors_t(i,jj) : d_neighbors(i,jj);
      pair_contrib(i, xtmp, ytmp, ztmp, qtmp, itype, jraw, fx, fy, fz);
    }, fxtmp, fytmp, fztmp);

    // one lane per atom writes the reduced i-force to the global array.  On a full
    // list this atom is the sole writer of f(i) AND the framework does not zero
    // f(i) on this path (it relies on the pair kernel to own it, cf. the
    // "NEIGHFLAG == FULL && ZEROFLAG" store in pair_kokkos.h), so assign rather
    // than add -- pair runs first among force styles, so nothing is clobbered.
    // On a half list other atoms' j-scatter also hit f(i), so it must be atomic.
    Kokkos::single(Kokkos::PerThread(team), [&] () {
      if (full) {
        f(i,0) = fxtmp;
        f(i,1) = fytmp;
        f(i,2) = fztmp;
      } else {
        Kokkos::atomic_add(&f(i,0), fxtmp);
        Kokkos::atomic_add(&f(i,1), fytmp);
        Kokkos::atomic_add(&f(i,2), fztmp);
      }
    });
  }
};

}  // anonymous namespace

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairLJCutCoulLong2Kokkos<DeviceType>::PairLJCutCoulLong2Kokkos(LAMMPS *lmp) :
  PairLJCutCoulLongKokkos<DeviceType>(lmp) {}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairLJCutCoulLong2Kokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
  // Optimized force-only paths (the melt hot path), selected by neighbor style:
  //   HALFTHREAD + newton on  -> half list, warp-per-atom, atomic j-scatter
  //   FULL                    -> full list, warp-per-atom, no atomics (2x pairs)
  // Everything else -- energy/virial steps for thermo, plain-half lists, many atom
  // types -- uses the proven base-class kernel unchanged.
  const bool use_half = (this->neighflag == HALFTHREAD) && this->force->newton_pair;
  const bool use_full = (this->neighflag == FULL);
  if (eflag_in || vflag_in ||
      (!use_half && !use_full) ||
      this->atom->ntypes > MAX_TYPES_STACKPARAMS) {
    PairLJCutCoulLongKokkos<DeviceType>::compute(eflag_in, vflag_in);
    return;
  }

  this->eflag = eflag_in;
  this->vflag = vflag_in;
  this->ev_init(eflag_in, vflag_in, 0);

  this->atomKK->sync(this->execution_space, X_MASK | F_MASK | TYPE_MASK | Q_MASK);
  this->atomKK->modified(this->execution_space, F_MASK);

  LJCL2Force<DeviceType> ff;
  ff.x    = this->atomKK->k_x.template view<DeviceType>();
  ff.f    = this->atomKK->k_f.template view<DeviceType>();
  ff.q    = this->atomKK->k_q.template view<DeviceType>();
  ff.type = this->atomKK->k_type.template view<DeviceType>();

  auto* k_list = static_cast<NeighListKokkos<DeviceType>*>(this->list);
  ff.d_ilist    = k_list->d_ilist;
  ff.d_numneigh = k_list->d_numneigh;
  ff.d_neighbors= k_list->d_neighbors;
  ff.d_neighbors_t = k_list->d_neighbors_transpose;  // empty view unless neigh/transpose on
  ff.use_transpose = this->lmp->kokkos->neigh_transpose;
  ff.inum       = k_list->inum;

  ff.full        = use_full;
  ff.nlocal      = this->atom->nlocal;
  ff.newton_pair = this->force->newton_pair;
  ff.g_ewald     = static_cast<KK_FLOAT>(this->g_ewald);
  ff.qqrd2e      = static_cast<KK_FLOAT>(this->force->qqrd2e);
  for (int k = 0; k < 4; k++) {
    ff.special_lj[k]   = static_cast<KK_FLOAT>(this->force->special_lj[k]);
    ff.special_coul[k] = static_cast<KK_FLOAT>(this->force->special_coul[k]);
  }

  const int nt = this->atom->ntypes;
  for (int i = 1; i <= nt; i++)
    for (int j = 1; j <= nt; j++) {
      ff.m_params[i][j]      = this->m_params[i][j];
      ff.m_cutsq[i][j]       = this->m_cutsq[i][j];
      ff.m_cut_ljsq[i][j]    = this->m_cut_ljsq[i][j];
      ff.m_cut_coulsq[i][j]  = this->m_cut_coulsq[i][j];
    }

  // warp-per-atom launch shape, from the standard KOKKOS package options -- the
  // same machinery as the stock neigh/thread path:
  //   threads/per/atom  -> vector_length (vector lanes splitting each atom's
  //                        neighbor loop; power of 2, <= 32)
  //   pair/team/size    -> team threads; atoms_per_team = team / vector_length
  // Tune at runtime with no rebuild, e.g.:
  //   -pk kokkos neigh half newton on binsize 8.0 neigh/transpose on \
  //      threads/per/atom 8 pair/team/size 64
  //
  // MEASURED OPTIMUM on A100 (melt, 505 neighs/atom half list, single precision):
  // threads/per/atom 8, pair/team/size 64 (= vector_length 8, atoms_per_team 8),
  // used as the DEFAULT when the user does not set these options -- the stock
  // package defaults (tpa 1, team 128) would pick vector_length 1, the worst shape
  // here (one-thread-per-atom, ~13.8 vs 26 ns/day).  Cost is monotonic in
  // vector_length (vl 8 -> 26, 16 -> 23.7, 32 -> 20.8): the kernel is limited by
  // the cross-lane reduction + atomic overhead, not gather latency or occupancy,
  // so fewer lanes per atom wins.  See kokkos_neigh.md sections 16.7, 16.11.  The
  // full-list path has ~2x the neighbors per atom, so its optimum may differ.
  auto *kk = this->lmp->kokkos;
  const int vector_length = kk->threads_per_atom_set ? kk->threads_per_atom : 8;
  const int team_size     = kk->pair_team_size_set   ? kk->pair_team_size   : 64;
  const int atoms_per_team = (team_size >= vector_length) ? team_size / vector_length : 1;
  ff.atoms_per_team = atoms_per_team;
  const int nteams = (k_list->inum + atoms_per_team - 1) / atoms_per_team;

  using policy_t = Kokkos::TeamPolicy<DeviceType, Kokkos::IndexType<int>>;
  policy_t policy(nteams, atoms_per_team, vector_length);
  Kokkos::parallel_for("PairLJCutCoulLong2::force", policy, ff);
}

namespace LAMMPS_NS {
template class PairLJCutCoulLong2Kokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairLJCutCoulLong2Kokkos<LMPHostType>;
#endif
}
