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
#include "comm.h"
#include "error.h"
#include "ewald_const.h"
#include "force.h"
#include "kokkos.h"
#include "neigh_list_kokkos.h"
#include "neighbor.h"
#include "update.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>

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

// ---------------------------------------------------------------------------
// i-batched union neighbor list: packing of one entry into one int.
//   bits  0..21  j                (requires nall < 2^22 = 4194304)
//   bits 22..29  member mask      (bit k = the k-th i atom of the group)
//   bits 30..31  special-bond code (only ever set on a single-member entry,
//                                   since the factors are per pair)
// ---------------------------------------------------------------------------
enum { UJ_BITS = 22, UMASK_SHIFT = 22, USB_SHIFT = 30, UPROBE = 8 };
static constexpr int UJ_MASK    = (1 << UJ_BITS) - 1;
static constexpr int UMASK_ALL  = 0xFF << UMASK_SHIFT;

KOKKOS_INLINE_FUNCTION int upopcount(int v) {
  int c = 0;
  for (int b = 0; b < 8; b++) c += (v >> b) & 1;
  return c;
}

// ---------------------------------------------------------------------------
// Packed j-side gather record (LMP_LJCL2_PACK=1).
//
// The flat kernel's per-entry cost is set by the NUMBER of memory instructions
// it issues, not by the bytes it moves.  Per entry it currently issues, from
// three unrelated arrays: x(j,0), x(j,1), x(j,2) -- three separate 32-bit loads,
// because a LayoutRight float*[3] has a 12-byte stride and nvcc cannot vectorize
// a 4-byte-aligned access -- plus type(j), plus q(j) on the ~24.5% of entries
// that pass the cutoff.  That is ~4.25 scattered loads x 210.5M entries/step.
//
// kokkos_neigh.md 18.2 measured this kernel at 47.5% issue with 68% of the warp
// slots RESIDENT and DRAM at only 14.7%: warps that are present but stalled on
// memory instructions whose data L2 is already absorbing.  That is a request-rate
// signature, not a bandwidth or an occupancy one.
//
// Packing (x,y,z,q) into one 16-byte-aligned record collapses four loads into a
// single LDG.128 that always lands inside one 32-byte sector, and demotes the
// type to a byte array whose sectors carry 32 atoms instead of 8.  Unconditional
// loads per entry: 4.25 -> 2.  Sectors touched: ~1.65 -> ~1.05.  The packed array
// is also one contiguous 16 B/atom region rather than three separate ones, so its
// L2 working set is strictly smaller than the three it replaces.
//
// LJCL2XQ itself lives in the header (the pair style owns the array).
// ---------------------------------------------------------------------------

template<class DeviceType>
struct LJCL2Pack {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  typename AT::t_kkfloat_1d_3_lr_randomread x;
  typename AT::t_kkfloat_1d_randomread q;
  typename AT::t_int_1d_randomread type;
  Kokkos::View<LJCL2XQ*, DeviceType> d_xq;
  Kokkos::View<unsigned char*, DeviceType> d_type8;

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i) const {
    LJCL2XQ r;
    r.x = x(i,0); r.y = x(i,1); r.z = x(i,2); r.q = q(i);
    d_xq(i) = r;
    d_type8(i) = static_cast<unsigned char>(type(i));
  }
};

// ---------------------------------------------------------------------------
// Kernel variants (LMP_LJCL2_VARIANT).  1-3 are DIAGNOSTIC: they delete part of
// the pair body, so their forces are wrong and their thermo is meaningless.
// They exist because every cost attribution in kokkos_neigh.md so far is a fit
// to end-to-end timings, and the packed-gather result (section 21) contradicted
// it.  Running the same kernel with one component removed measures that
// component directly.  Variant 4 is a real candidate, not a diagnostic.
// ---------------------------------------------------------------------------
enum {
  LJCL2_VAR_NORMAL   = 0,
  LJCL2_VAR_NOATOMIC = 1,  // drop the j-side atomic scatter (i-force still summed)
  LJCL2_VAR_NOCOUL   = 2,  // drop the Ewald block (exp, divide, poly); keep atomics
  LJCL2_VAR_WALKONLY = 3,  // list walk + gather + distance test only
  LJCL2_VAR_FASTMATH = 4   // normal, with __expf / __fdividef for the Ewald block
};

// The build has no -use_fast_math (CMAKE_CXX_FLAGS is empty), so `Kokkos::exp`
// is accurate expf and `1.0f/x` is an IEEE-correct division -- roughly a dozen
// instructions each, both inside the in-cutoff branch.  Over the argument
// ranges here (-grij^2 in [-12,0]; 1 + EWALD_P*grij in [1,4]) the intrinsics are
// accurate to ~2 ulp, well below the ~1e-7 accuracy of the Abramowitz-Stegun
// erfc approximation they feed, so variant 4 is a free-lunch candidate rather
// than a precision trade.  Guarded on sizeof so a double build folds it away.
template<bool FAST>
KOKKOS_INLINE_FUNCTION KK_FLOAT ljcl2_exp(const KK_FLOAT v) {
#if defined(__CUDA_ARCH__)
  if (FAST && sizeof(KK_FLOAT) == sizeof(float))
    return static_cast<KK_FLOAT>(__expf(static_cast<float>(v)));
#endif
  return Kokkos::exp(v);
}

template<bool FAST>
KOKKOS_INLINE_FUNCTION KK_FLOAT ljcl2_recip(const KK_FLOAT v) {
#if defined(__CUDA_ARCH__)
  if (FAST && sizeof(KK_FLOAT) == sizeof(float))
    return static_cast<KK_FLOAT>(__fdividef(1.0f, static_cast<float>(v)));
#endif
  return static_cast<KK_FLOAT>(1.0)/v;
}

// Per-pair physics, shared by the flat and the union kernel so the two cannot
// drift apart.  Both evaluate in the same order, so they agree bit for bit.
template<class DeviceType>
struct LJCL2Common {
  typedef ArrayTypes<DeviceType> AT;

  typename AT::t_kkfloat_1d_3_lr_randomread x;
  typename AT::t_kkacc_1d_3 f;
  typename AT::t_kkfloat_1d_randomread q;
  typename AT::t_int_1d_randomread type;
  typename AT::t_int_1d d_ilist;

  int nlocal, newton_pair, inum;
  KK_FLOAT g_ewald, qqrd2e;
  KK_FLOAT special_lj[4], special_coul[4];
  params_lj_coul m_params[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
  KK_FLOAT m_cutsq[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
  KK_FLOAT m_cut_ljsq[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
  KK_FLOAT m_cut_coulsq[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];

  // fpair for one pair already known to satisfy rsq < m_cutsq[itype][jtype].
  // FAST selects the intrinsic exp/reciprocal; NOCOUL drops the Ewald block
  // (diagnostic variant 2 only -- the result is not physical).
  template<bool FAST = false, bool NOCOUL = false>
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT pair_fpair(const KK_FLOAT rsq, const int itype, const int jtype,
                      const KK_FLOAT qiqj, const int sb) const {
    // one reciprocal-sqrt supplies rinv, r2inv, and r (was 2 divides + 1 sqrt)
    const KK_FLOAT rinv = Kokkos::rsqrt(rsq);
    const KK_FLOAT r2inv = rinv*rinv;
    KK_FLOAT fpair = static_cast<KK_FLOAT>(0.0);

    if (rsq < m_cut_ljsq[itype][jtype]) {
      const KK_FLOAT r6inv = r2inv*r2inv*r2inv;
      fpair += special_lj[sb] *
        r6inv*(m_params[itype][jtype].lj1*r6inv - m_params[itype][jtype].lj2) * r2inv;
    }

    if (!NOCOUL && rsq < m_cut_coulsq[itype][jtype]) {
      // analytical Ewald real-space correction (Abramowitz-Stegun 7.1.26)
      const KK_FLOAT r = rsq*rinv;
      const KK_FLOAT grij = g_ewald * r;
      const KK_FLOAT expm2 = ljcl2_exp<FAST>(-grij*grij);
      const KK_FLOAT t = ljcl2_recip<FAST>(
        static_cast<KK_FLOAT>(1.0) + static_cast<KK_FLOAT>(EWALD_P)*grij);
      const KK_FLOAT erfc = t * (static_cast<KK_FLOAT>(A1)+t*(static_cast<KK_FLOAT>(A2)+
                            t * (static_cast<KK_FLOAT>(A3)+t*(static_cast<KK_FLOAT>(A4)+
                            t * static_cast<KK_FLOAT>(A5))))) * expm2;
      const KK_FLOAT prefactor = qqrd2e * qiqj * rinv;
      KK_FLOAT forcecoul = prefactor * (erfc + static_cast<KK_FLOAT>(EWALD_F)*grij*expm2);
      const KK_FLOAT factor_coul = special_coul[sb];
      if (factor_coul < static_cast<KK_FLOAT>(1.0))
        forcecoul -= (static_cast<KK_FLOAT>(1.0)-factor_coul)*prefactor;
      fpair += forcecoul * r2inv;
    }
    return fpair;
  }
};

template<class DeviceType, bool FULL, bool DUAL, bool PACK, int VAR = LJCL2_VAR_NORMAL>
struct LJCL2Force : public LJCL2Common<DeviceType> {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef Kokkos::TeamPolicy<DeviceType, Kokkos::IndexType<int>> policy_type;
  typedef typename policy_type::member_type member_type;

  typename AT::t_int_1d d_numneigh;
  typename AT::t_neighbors_2d d_neighbors;        // LayoutLeft  (atom-fast)
  typename AT::t_neighbors_2d_lr d_neighbors_t;   // LayoutRight (jj-fast), transpose

  // packed j-side gather (PACK only); empty views otherwise
  Kokkos::View<const LJCL2XQ*, DeviceType> d_xq;
  Kokkos::View<const unsigned char*, DeviceType> d_type8;

  // fused dual-cutoff inner list (see the header); dual==0 leaves every branch
  // below compiled out of the hot path by the uniform-across-warp flags
  typename AT::t_neighbors_2d_lr d_inbr;
  typename AT::t_int_1d d_innum;
  typename AT::t_int_1d d_ictl;
  KK_FLOAT inner_cutsq;

  int atoms_per_team;
  bool use_transpose;  // read neighbor indices from the LayoutRight transpose list.
                       // In warp-per-atom, the vector lanes read d_neighbors(i,jj)
                       // for consecutive jj; LayoutLeft strides those by nmax
                       // (fully uncoalesced), LayoutRight makes them contiguous.

  // Force contribution of one (i,jj) neighbor pair: always folds the i-force into
  // (fx,fy,fz).  On a HALF list (full==false) it also atomic-scatters -f to j.
  // On a FULL list (full==true) j supplies its own force from its own row, so the
  // scatter is skipped -- that removes every per-pair atomic, the throughput wall
  // on A100, at the cost of computing each pair twice.
  // With store!=0 the pair is also appended to the inner list when it is within
  // the inner cutoff -- free apart from the store itself, since rsq is already in
  // a register.  The inner row has the same capacity as the master row and the
  // inner list is a subset of it, so the append cannot overflow.
  //
  // The append is a contended per-atom atomic.  An atomic-free variant (lane L
  // owns jj = L, L+vl, ... and writes keepers at L + vl*t) was built and
  // measured 10% SLOWER: splitting the loop into a refresh branch and a tight
  // branch inlines this body twice (REG 43 -> 55) and gives the tight loop a
  // per-lane trip count that nvcc can no longer unroll, which costs more than
  // the atomics saved.  See kokkos_neigh.md 20.8 before trying it again.
  KOKKOS_INLINE_FUNCTION
  void pair_contrib(const int i, const KK_FLOAT xtmp, const KK_FLOAT ytmp,
                    const KK_FLOAT ztmp, const KK_FLOAT qtmp, const int itype,
                    const int jraw, KK_ACC_FLOAT& fx, KK_ACC_FLOAT& fy,
                    KK_ACC_FLOAT& fz, const int store) const {
    const int sb = (jraw >> SBBITS) & 3;
    const int j = jraw & NEIGHMASK;

    // One 16-byte load for (x,y,z,q) plus one byte for the type, instead of
    // three 32-bit loads from x, one from type and one from q.  qj is carried
    // in a register either way; in the unpacked path q(j) is deliberately read
    // only inside the cutoff branch, in the packed path it rides along free.
    KK_FLOAT xj, yj, zj, qj;
    int jtype;
    if (PACK) {
      const LJCL2XQ j4 = d_xq(j);
      xj = j4.x; yj = j4.y; zj = j4.z; qj = j4.q;
      jtype = static_cast<int>(d_type8(j));
    } else {
      xj = this->x(j,0); yj = this->x(j,1); zj = this->x(j,2);
      qj = static_cast<KK_FLOAT>(0.0);
      jtype = this->type(j);
    }

    const KK_FLOAT delx = xtmp - xj;
    const KK_FLOAT dely = ytmp - yj;
    const KK_FLOAT delz = ztmp - zj;
    const KK_FLOAT rsq = delx*delx + dely*dely + delz*delz;

    if (DUAL && store && rsq < inner_cutsq)
      d_inbr(i, Kokkos::atomic_fetch_add(&d_innum(i), 1)) = jraw;

    // Diagnostic variant 3: stop after the gather and the distance test.  The
    // in-cutoff accumulation of the raw separation keeps the loads and the
    // compare from being dead-code-eliminated without costing a divide, an exp
    // or an atomic, so the arm measures exactly "walk + gather + test".
    if (VAR == LJCL2_VAR_WALKONLY) {
      if (rsq < this->m_cutsq[itype][jtype]) {
        fx += static_cast<KK_ACC_FLOAT>(delx);
        fy += static_cast<KK_ACC_FLOAT>(dely);
        fz += static_cast<KK_ACC_FLOAT>(delz);
      }
      return;
    }

    if (rsq < this->m_cutsq[itype][jtype]) {
      const KK_FLOAT fpair =
        this->template pair_fpair<VAR == LJCL2_VAR_FASTMATH, VAR == LJCL2_VAR_NOCOUL>(
          rsq, itype, jtype, qtmp*(PACK ? qj : this->q(j)), sb);

      fx += static_cast<KK_ACC_FLOAT>(delx*fpair);
      fy += static_cast<KK_ACC_FLOAT>(dely*fpair);
      fz += static_cast<KK_ACC_FLOAT>(delz*fpair);

      if (VAR != LJCL2_VAR_NOATOMIC &&
          !FULL && (this->newton_pair || j < this->nlocal)) {
        Kokkos::atomic_add(&this->f(j,0), -static_cast<KK_ACC_FLOAT>(delx*fpair));
        Kokkos::atomic_add(&this->f(j,1), -static_cast<KK_ACC_FLOAT>(dely*fpair));
        Kokkos::atomic_add(&this->f(j,2), -static_cast<KK_ACC_FLOAT>(delz*fpair));
      }
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const member_type& team) const {
    const int atom_in_team = team.league_rank()*atoms_per_team + team.team_rank();
    if (atom_in_team >= this->inum) return;

    const int i = this->d_ilist(atom_in_team);
    KK_FLOAT xtmp, ytmp, ztmp, qtmp;
    int itype;
    if (PACK) {
      const LJCL2XQ i4 = d_xq(i);
      xtmp = i4.x; ytmp = i4.y; ztmp = i4.z; qtmp = i4.q;
      itype = static_cast<int>(d_type8(i));
    } else {
      xtmp = this->x(i,0); ytmp = this->x(i,1); ztmp = this->x(i,2);
      qtmp = this->q(i);
      itype = this->type(i);
    }

    // dual-cutoff mode: walk (and refill) the master list on a refresh step,
    // walk the inner list otherwise.  d_ictl(0) is uniform across the launch,
    // so the branch is warp-uniform and free.
    const int refresh = DUAL ? d_ictl(0) : 1;
    const int inner = DUAL && !refresh;
    const int jnum = inner ? d_innum(i) : d_numneigh(i);
    const int store = DUAL && refresh;

    KK_ACC_FLOAT fxtmp = 0.0, fytmp = 0.0, fztmp = 0.0;

    // vector lanes of this atom's warp split the neighbor loop; the i-force
    // contributions are reduced across lanes into (fxtmp,fytmp,fztmp).  (A 2-way
    // ILP unroll here regressed the A100 ~30%: it is throughput/atomic-bound, not
    // latency-bound, so extra in-flight gathers only deepen the memory queues.)
    Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(team, jnum),
      [&] (const int jj, KK_ACC_FLOAT& fx, KK_ACC_FLOAT& fy, KK_ACC_FLOAT& fz) {
      // LayoutRight transpose makes these consecutive-jj reads coalesce across the
      // warp's vector lanes; LayoutLeft strides them by nmax (uncoalesced).
      const int jraw = (DUAL && inner) ? d_inbr(i,jj)
                     : (use_transpose ? d_neighbors_t(i,jj) : d_neighbors(i,jj));
      pair_contrib(i, xtmp, ytmp, ztmp, qtmp, itype, jraw, fx, fy, fz, store);
    }, fxtmp, fytmp, fztmp);

    // one lane per atom writes the reduced i-force to the global array.  On a full
    // list this atom is the sole writer of f(i) AND the framework does not zero
    // f(i) on this path (it relies on the pair kernel to own it, cf. the
    // "NEIGHFLAG == FULL && ZEROFLAG" store in pair_kokkos.h), so assign rather
    // than add -- pair runs first among force styles, so nothing is clobbered.
    // On a half list other atoms' j-scatter also hit f(i), so it must be atomic.
    Kokkos::single(Kokkos::PerThread(team), [&] () {
      if (FULL) {
        this->f(i,0) = fxtmp;
        this->f(i,1) = fytmp;
        this->f(i,2) = fztmp;
      } else {
        Kokkos::atomic_add(&this->f(i,0), fxtmp);
        Kokkos::atomic_add(&this->f(i,1), fytmp);
        Kokkos::atomic_add(&this->f(i,2), fztmp);
      }
    });
  }
};

/* ----------------------------------------------------------------------
   Dual-cutoff control, entirely on the device: three trivial kernels run
   before the force kernel each step.  The refresh decision is never read
   back to the host, so the scheme adds no per-step fence.

   The criterion is LAMMPS's own: refresh once any atom has moved more than
   half the inner skin, since two atoms can then have closed by the full
   skin.  The reduction covers locals AND ghosts -- ghost coordinates are
   forward-communicated every step, so a ghost's displacement is seen
   directly and the check stays exact on any number of ranks with no
   communication.
------------------------------------------------------------------------- */

template<class DeviceType>
struct LJCL2DualDisp {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef KK_FLOAT value_type;

  typename AT::t_kkfloat_1d_3_lr_randomread x;
  typename AT::t_kkfloat_1d_3_lr xhold;

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, KK_FLOAT& m) const {
    const KK_FLOAT dx = x(i,0) - xhold(i,0);
    const KK_FLOAT dy = x(i,1) - xhold(i,1);
    const KK_FLOAT dz = x(i,2) - xhold(i,2);
    const KK_FLOAT r2 = dx*dx + dy*dy + dz*dz;
    if (r2 > m) m = r2;
  }
};

template<class DeviceType>
struct LJCL2DualDecide {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  typename AT::t_int_1d d_ictl;
  Kokkos::View<KK_FLOAT, DeviceType> d_maxsq;
  KK_FLOAT thresh2;
  int every, force_refresh;

  KOKKOS_INLINE_FUNCTION
  void operator()(const int) const {
    const int since = d_ictl(1);
    const int over = (d_maxsq() > thresh2) ? 1 : 0;
    const int go = (force_refresh || (since + 1 >= every) || over) ? 1 : 0;
    d_ictl(0) = go;
    d_ictl(1) = go ? 0 : since + 1;
    if (go) {
      d_ictl(2) += 1;
      if (over && !force_refresh) d_ictl(3) += 1;
    }
  }
};

// Correctness check for the dual cutoff (LMP_LJCL2_DUAL_CHECK=N, debug only).
// Counts, at the CURRENT positions, the in-force-cutoff pairs of the master row
// and of the inner row of every atom.  The inner list is a subset of the master
// list, so the two counts are equal if and only if the inner list is missing no
// pair that the force kernel would have to compute.  Integer counts, so the
// test is exact and immune to the atomic-summation nondeterminism that makes
// thermo comparisons only ~6 digits meaningful.
template<class DeviceType>
struct LJCL2DualCheck {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  typename AT::t_kkfloat_1d_3_lr_randomread x;
  typename AT::t_int_1d_randomread type;
  typename AT::t_int_1d d_ilist, d_numneigh, d_innum, d_ictl;
  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_neighbors_2d_lr d_inbr;
  KK_FLOAT m_cutsq[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];

  KOKKOS_INLINE_FUNCTION
  int in_cutoff(const int i, const int j) const {
    const KK_FLOAT dx = x(i,0) - x(j,0);
    const KK_FLOAT dy = x(i,1) - x(j,1);
    const KK_FLOAT dz = x(i,2) - x(j,2);
    return (dx*dx + dy*dy + dz*dz < m_cutsq[type(i)][type(j)]) ? 1 : 0;
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const int ii) const {
    const int i = d_ilist(ii);
    int cm = 0, ci = 0;
    const int jnum = d_numneigh(i);
    for (int jj = 0; jj < jnum; jj++) cm += in_cutoff(i, d_neighbors(i,jj) & NEIGHMASK);
    const int inum_i = d_innum(i);
    for (int jj = 0; jj < inum_i; jj++) ci += in_cutoff(i, d_inbr(i,jj) & NEIGHMASK);
    Kokkos::atomic_add(&d_ictl(4), cm);
    Kokkos::atomic_add(&d_ictl(5), ci);
  }
};

template<class DeviceType>
struct LJCL2SumInt {   // diagnostics only (thermo steps)
  typedef DeviceType device_type;
  typename ArrayTypes<DeviceType>::t_int_1d v;
  KOKKOS_INLINE_FUNCTION void operator()(const int i, int& s) const { s += v(i); }
};

template<class DeviceType>
struct LJCL2DualReset {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  typename AT::t_kkfloat_1d_3_lr_randomread x;
  typename AT::t_kkfloat_1d_3_lr xhold;
  typename AT::t_int_1d d_innum;
  typename AT::t_int_1d d_ictl;
  int nlocal;

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i) const {
    if (!d_ictl(0)) return;
    xhold(i,0) = x(i,0);
    xhold(i,1) = x(i,1);
    xhold(i,2) = x(i,2);
    if (i < nlocal) d_innum(i) = 0;
  }
};

/* ----------------------------------------------------------------------
   Union-list force kernel: one team slot per GROUP of CI i atoms, the
   group's merged entry list split across the vector lanes.

   Per entry the lane loads j's index once, gathers x/q/type once, and runs
   the pairs of the (up to CI) member atoms whose mask bit is set, summing
   their reaction force in registers -- so ONE atomic triple lands on f(j)
   instead of one per pair.  Measured mean popcount is the divisor on the
   index load, the gather and the atomic; the distance tests and the
   in-cutoff LJ/Ewald work are unchanged (the mask reproduces the flat list
   exactly, hence no lane inflation of the kind that sank the cluster
   kernel, kokkos_neigh.md section 13).
------------------------------------------------------------------------- */

template<int N>
struct UAccum {
  KK_ACC_FLOAT v[N];
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION UAccum() { for (int k = 0; k < N; k++) v[k] = 0.0; }
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION void operator+=(const UAccum &rhs) {
    for (int k = 0; k < N; k++) v[k] += rhs.v[k];
  }
};

template<class DeviceType, int CI>
struct LJCL2UnionForce : public LJCL2Common<DeviceType> {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef Kokkos::TeamPolicy<DeviceType, Kokkos::IndexType<int>> policy_type;
  typedef typename policy_type::member_type member_type;

  typename AT::t_neighbors_2d_lr d_unbr;
  typename AT::t_int_1d d_unum;
  int ngroup, groups_per_team;

  KOKKOS_INLINE_FUNCTION
  void operator()(const member_type& team) const {
    const int g = team.league_rank()*groups_per_team + team.team_rank();
    if (g >= ngroup) return;
    const int p0 = g*CI;
    const int nmem = (this->inum - p0 < CI) ? (this->inum - p0) : CI;

    int ii[CI], ti[CI];
    KK_FLOAT xi[CI], yi[CI], zi[CI], qi[CI];
#pragma unroll
    for (int k = 0; k < CI; k++) {
      const int i = this->d_ilist((k < nmem) ? p0 + k : p0);
      ii[k] = i;
      xi[k] = this->x(i,0); yi[k] = this->x(i,1); zi[k] = this->x(i,2);
      qi[k] = this->q(i);   ti[k] = this->type(i);
    }

    UAccum<3*CI> fi;
    const int unum = d_unum(g);

    Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(team, unum),
      [&] (const int e, UAccum<3*CI>& acc) {
      // LayoutRight: consecutive e across the vector lanes are contiguous ints
      const int w = d_unbr(g,e);
      const int j = w & UJ_MASK;
      const int mask = (w >> UMASK_SHIFT) & 0xFF;
      const int sb = (w >> USB_SHIFT) & 3;
      const KK_FLOAT xj = this->x(j,0);
      const KK_FLOAT yj = this->x(j,1);
      const KK_FLOAT zj = this->x(j,2);
      const KK_FLOAT qj = this->q(j);
      const int jtype = this->type(j);
      KK_ACC_FLOAT fjx = 0.0, fjy = 0.0, fjz = 0.0;

#pragma unroll
      for (int k = 0; k < CI; k++) {
        if (!((mask >> k) & 1)) continue;
        const KK_FLOAT delx = xi[k] - xj;
        const KK_FLOAT dely = yi[k] - yj;
        const KK_FLOAT delz = zi[k] - zj;
        const KK_FLOAT rsq = delx*delx + dely*dely + delz*delz;
        if (rsq >= this->m_cutsq[ti[k]][jtype]) continue;
        const KK_FLOAT fpair = this->pair_fpair(rsq, ti[k], jtype, qi[k]*qj, sb);
        const KK_ACC_FLOAT fx = static_cast<KK_ACC_FLOAT>(delx*fpair);
        const KK_ACC_FLOAT fy = static_cast<KK_ACC_FLOAT>(dely*fpair);
        const KK_ACC_FLOAT fz = static_cast<KK_ACC_FLOAT>(delz*fpair);
        acc.v[3*k+0] += fx; acc.v[3*k+1] += fy; acc.v[3*k+2] += fz;
        fjx += fx; fjy += fy; fjz += fz;
      }

      if (this->newton_pair || j < this->nlocal) {
        if (fjx != 0.0 || fjy != 0.0 || fjz != 0.0) {
          Kokkos::atomic_add(&this->f(j,0), -fjx);
          Kokkos::atomic_add(&this->f(j,1), -fjy);
          Kokkos::atomic_add(&this->f(j,2), -fjz);
        }
      }
    }, fi);

    Kokkos::single(Kokkos::PerThread(team), [&] () {
      for (int k = 0; k < nmem; k++) {
        Kokkos::atomic_add(&this->f(ii[k],0), fi.v[3*k+0]);
        Kokkos::atomic_add(&this->f(ii[k],1), fi.v[3*k+1]);
        Kokkos::atomic_add(&this->f(ii[k],2), fi.v[3*k+2]);
      }
    });
  }
};

/* ----------------------------------------------------------------------
   Union-list build: merge the CI flat half-list rows of one group into one
   entry per distinct j.  One team per group; the team's threads stream the
   rows and insert into a shared-memory open-addressing table keyed on
   j & (nslot-1) -- the low bits, so runs of consecutive j (which is what a
   spatially sorted list produces) land in consecutive slots and the
   slot-ordered output keeps the flat list's gather locality.

   Correctness does not depend on the merge succeeding: an entry that misses
   the table (probe limit) or that carries special-bond bits is emitted as a
   single-member entry, so every flat-list pair appears in exactly one entry
   either way.  That is why the table needs no fatal saturation check.
------------------------------------------------------------------------- */

template<class DeviceType, int CI>
struct LJCL2UnionBuild {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef Kokkos::TeamPolicy<DeviceType, Kokkos::IndexType<int>> policy_type;
  typedef typename policy_type::member_type member_type;
  typedef Kokkos::View<int*, typename DeviceType::scratch_memory_space,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>> t_sh_int_1d;

  typename AT::t_int_1d d_ilist;
  typename AT::t_int_1d d_numneigh;
  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_neighbors_2d_lr d_neighbors_t;
  typename AT::t_neighbors_2d_lr d_unbr;
  typename AT::t_int_1d d_unum;
  typename AT::t_int_1d d_ustat;
  int inum, maxu, nslot;
  bool use_transpose;

  KOKKOS_INLINE_FUNCTION
  void operator()(const member_type& team) const {
    const int g = team.league_rank();
    const int p0 = g*CI;
    const int nmem = (inum - p0 < CI) ? (inum - p0) : CI;
    const int slot_mask = nslot - 1;

    t_sh_int_1d s_tab(team.team_scratch(0), nslot);
    t_sh_int_1d s_ctr(team.team_scratch(0), 1);

    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nslot),
                         [&] (const int s) { s_tab(s) = 0; });
    Kokkos::single(Kokkos::PerTeam(team), [&] () { s_ctr(0) = 0; });
    team.team_barrier();

    // ---- insert phase: all table operations are atomic, so the member rows
    // need no barrier between them
    for (int k = 0; k < nmem; k++) {
      const int i = d_ilist(p0+k);
      const int jnum = d_numneigh(i);
      const int kbit = 1 << (UMASK_SHIFT + k);
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, jnum), [&] (const int jj) {
        const int jraw = use_transpose ? d_neighbors_t(i,jj) : d_neighbors(i,jj);
        const int sb = (jraw >> SBBITS) & 3;
        const int j = jraw & NEIGHMASK;
        int placed = 0;
        if (!sb) {
          int slot = j & slot_mask;
          for (int p = 0; p < UPROBE; p++) {
            int w = s_tab(slot);
            if (w == 0) {
              const int old = Kokkos::atomic_compare_exchange(&s_tab(slot), 0, (j+1) | kbit);
              if (old == 0) { placed = 1; break; }   // claimed the slot, bit is set
              w = old;
            }
            if ((w & UJ_MASK) == j+1) {
              Kokkos::atomic_fetch_or(&s_tab(slot), kbit);
              placed = 1; break;
            }
            slot = (slot+1) & slot_mask;
          }
        }
        if (!placed) {
          // single-member entry; fill the group's capacity from the top so the
          // merged entries can still be written in slot order from index 0
          const int c = Kokkos::atomic_fetch_add(&s_ctr(0), 1);
          if (c < maxu) d_unbr(g, maxu-1-c) = j | kbit | (sb << USB_SHIFT);
        }
      });
    }
    team.team_barrier();

    const int nspec = s_ctr(0);

    // ---- count occupied slots (and the pair total, for the stats line)
    int nuniq = 0, npair = 0;
    Kokkos::parallel_reduce(Kokkos::TeamThreadRange(team, nslot),
      [&] (const int s, int& cnt, int& pc) {
      const int w = s_tab(s);
      if (w) { cnt++; pc += upopcount((w & UMASK_ALL) >> UMASK_SHIFT); }
    }, nuniq, npair);

    const int total = nuniq + nspec;
    // the specials still sit at the top; they must not be overwritten before
    // they are moved down, hence capacity for both copies
    const int need = nuniq + 2*nspec;

    if (need <= maxu) {
      // ---- write the merged entries in slot order (locality) ...
      Kokkos::parallel_scan(Kokkos::TeamThreadRange(team, nslot),
        [&] (const int s, int& upd, const bool final) {
        const int w = s_tab(s);
        const int occ = w ? 1 : 0;
        if (final && occ) d_unbr(g, upd) = ((w & UJ_MASK) - 1) | (w & UMASK_ALL);
        upd += occ;
      });
      team.team_barrier();
      // ... then move the specials in behind them
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nspec),
        [&] (const int c) { d_unbr(g, nuniq+c) = d_unbr(g, maxu-1-c); });
    }

    Kokkos::single(Kokkos::PerTeam(team), [&] () {
      d_unum(g) = (need <= maxu) ? total : 0;
      if (need > maxu) Kokkos::atomic_fetch_or(&d_ustat(0), 1);
      Kokkos::atomic_fetch_max(&d_ustat(1), need);
      Kokkos::atomic_add(&d_ustat(2), total);
      Kokkos::atomic_add(&d_ustat(3), npair + nspec);
      Kokkos::atomic_add(&d_ustat(4), nspec);
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
    // thermo steps walk the master list in the base kernel, so the inner list
    // simply goes unused; report its counters here, where a sync is free
    if (dual_stats && dual_on == 1 && d_innum.extent(0) > 0) {
      k_ictl.template modify<DeviceType>();
      k_ictl.template sync<LMPHostType>();
      LJCL2SumInt<DeviceType> fs;
      fs.v = d_innum;
      int total = 0;
      Kokkos::parallel_reduce("PairLJCutCoulLong2::dual_stats",
        Kokkos::RangePolicy<DeviceType>(0, this->atom->nlocal), fs, total);
      if (this->comm->me == 0)
        printf("DUAL_STATS: step " BIGINT_FORMAT " refreshes %d displacement-triggered %d "
               "inner entries/atom %.1f\n", this->update->ntimestep,
               k_ictl.view_host()(2), k_ictl.view_host()(3),
               (double) total / this->atom->nlocal);
    }
    PairLJCutCoulLongKokkos<DeviceType>::compute(eflag_in, vflag_in);
    return;
  }

  this->eflag = eflag_in;
  this->vflag = vflag_in;
  this->ev_init(eflag_in, vflag_in, 0);

  this->atomKK->sync(this->execution_space, X_MASK | F_MASK | TYPE_MASK | Q_MASK);
  this->atomKK->modified(this->execution_space, F_MASK);

  auto* k_list = static_cast<NeighListKokkos<DeviceType>*>(this->list);

  // ---- experimental i-batched union list (half list only).  LMP_LJCL2_UNION
  // gives CI, the number of consecutive i atoms that share one entry per
  // distinct j; 0 (default) keeps the per-atom flat walk.
  if (union_ci < 0) {
    union_ci = 0;
    const char *e = std::getenv("LMP_LJCL2_UNION");
    if (e) {
      const int v = std::atoi(e);
      if (v == 2 || v == 4 || v == 8) union_ci = v;
      else if (v != 0 && this->comm->me == 0)
        this->error->warning(FLERR, "LMP_LJCL2_UNION must be 2, 4 or 8 -- ignored");
    }
    union_stats = (std::getenv("LMP_UNION_STATS") != nullptr);
    // the entry format carries j in 22 bits
    if (union_ci && (this->atom->nlocal + this->atom->nghost) >= UJ_MASK) {
      if (this->comm->me == 0)
        this->error->warning(FLERR, "too many atoms for the union list format -- disabled");
      union_ci = 0;
    }
  }

  if (use_half && union_ci) {
    if (this->neighbor->lastcall > union_built_step) {
      if (union_ci == 2) union_build<2>();
      else if (union_ci == 4) union_build<4>();
      else union_build<8>();
      union_built_step = this->neighbor->lastcall;
    }
    if (union_ci == 2) union_compute<2>();
    else if (union_ci == 4) union_compute<4>();
    else union_compute<8>();
    return;
  }

  // ---- experimental fused dual-cutoff inner list, half list only
  if (dual_on < 0) {
    dual_on = 0;
    if (std::getenv("LMP_LJCL2_DUAL")) dual_on = 1;
    if (const char *e = std::getenv("LMP_LJCL2_DUAL_EVERY")) dual_every = std::atoi(e);
    if (const char *e = std::getenv("LMP_LJCL2_DUAL_SKIN")) dual_skin = std::atof(e);
    dual_stats = (std::getenv("LMP_DUAL_STATS") != nullptr);
    if (const char *e = std::getenv("LMP_LJCL2_DUAL_CHECK")) dual_check = std::atoi(e);
    if (dual_every < 1) dual_every = 1;
    if (dual_on && (dual_skin <= 0.0 || dual_skin >= this->neighbor->skin)) {
      if (this->comm->me == 0)
        this->error->warning(FLERR, "LMP_LJCL2_DUAL_SKIN must be > 0 and < the neighbor "
                                    "skin -- dual cutoff disabled");
      dual_on = 0;
    }
    if (dual_on && union_ci > 0) {
      if (this->comm->me == 0)
        this->error->warning(FLERR, "LMP_LJCL2_UNION and LMP_LJCL2_DUAL are exclusive "
                                    "-- dual cutoff disabled");
      dual_on = 0;
    }
  }

  bool force_refresh = false;
  const bool dual = use_half && (dual_on == 1);
  if (dual) {
    dual_setup(force_refresh);
    dual_decide(force_refresh);
  }

  // ---- experimental packed j-side gather
  if (pack_on < 0) {
    pack_on = 0;
    if (const char *e = std::getenv("LMP_LJCL2_PACK")) pack_on = std::atoi(e) ? 1 : 0;
  }
  if (pack_on) pack_refresh();

  // ---- kernel variant.  1-3 are diagnostics whose forces are WRONG by
  // construction; they are only ever selected explicitly, and the style warns
  // once so no benchmark log can be mistaken for a physical run.
  if (variant < 0) {
    variant = LJCL2_VAR_NORMAL;
    if (const char *e = std::getenv("LMP_LJCL2_VARIANT")) variant = std::atoi(e);
    if (variant < 0 || variant > LJCL2_VAR_FASTMATH) {
      if (this->comm->me == 0)
        this->error->warning(FLERR, "LMP_LJCL2_VARIANT out of range -- ignored");
      variant = LJCL2_VAR_NORMAL;
    }
    if (variant >= LJCL2_VAR_NOATOMIC && variant <= LJCL2_VAR_WALKONLY &&
        this->comm->me == 0)
      this->error->warning(FLERR, "lj/cut/coul/long2/kk is running a DIAGNOSTIC "
                                  "kernel variant: the forces it produces are "
                                  "incomplete and the trajectory is not physical");
  }

  if (dual) {
    if (pack_on) flat_compute<false,true,true>();
    else         flat_compute<false,true,false>();
  } else if (use_full) {
    if (pack_on) flat_compute<true,false,true>();
    else         flat_compute<true,false,false>();
  } else if (pack_on) {
    switch (variant) {
      case LJCL2_VAR_NOATOMIC: flat_compute<false,false,true,LJCL2_VAR_NOATOMIC>(); break;
      case LJCL2_VAR_NOCOUL:   flat_compute<false,false,true,LJCL2_VAR_NOCOUL>();   break;
      case LJCL2_VAR_WALKONLY: flat_compute<false,false,true,LJCL2_VAR_WALKONLY>(); break;
      case LJCL2_VAR_FASTMATH: flat_compute<false,false,true,LJCL2_VAR_FASTMATH>(); break;
      default:                 flat_compute<false,false,true>();                    break;
    }
  } else {
    switch (variant) {
      case LJCL2_VAR_NOATOMIC: flat_compute<false,false,false,LJCL2_VAR_NOATOMIC>(); break;
      case LJCL2_VAR_NOCOUL:   flat_compute<false,false,false,LJCL2_VAR_NOCOUL>();   break;
      case LJCL2_VAR_WALKONLY: flat_compute<false,false,false,LJCL2_VAR_WALKONLY>(); break;
      case LJCL2_VAR_FASTMATH: flat_compute<false,false,false,LJCL2_VAR_FASTMATH>(); break;
      default:                 flat_compute<false,false,false>();                    break;
    }
  }

  if (dual && dual_check && (this->update->ntimestep % dual_check) == 0) dual_verify();
}

/* ----------------------------------------------------------------------
   flat (per-atom) force kernel launch.  FULL/DUAL are compile time so that
   the plain half-list instantiation is identical to the pre-dual-cutoff code.
------------------------------------------------------------------------- */

template<class DeviceType>
template<bool FULL, bool DUAL, bool PACK, int VAR>
void PairLJCutCoulLong2Kokkos<DeviceType>::flat_compute()
{
  auto* k_list = static_cast<NeighListKokkos<DeviceType>*>(this->list);

  LJCL2Force<DeviceType,FULL,DUAL,PACK,VAR> ff;
  fill_common(ff);
  ff.d_numneigh = k_list->d_numneigh;
  ff.d_neighbors= k_list->d_neighbors;
  ff.d_neighbors_t = k_list->d_neighbors_transpose;  // empty view unless neigh/transpose on
  ff.use_transpose = this->lmp->kokkos->neigh_transpose;
  if (PACK) {
    ff.d_xq = d_xq;
    ff.d_type8 = d_type8;
  }
  if (DUAL) {
    ff.d_inbr  = d_inbr;
    ff.d_innum = d_innum;
    ff.d_ictl  = k_ictl.template view<DeviceType>();
    ff.inner_cutsq = static_cast<KK_FLOAT>(inner_cutsq);
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
  // full-list path has ~2x the neighbors per atom, so its optimum may differ; the
  // dual-cutoff inner list has ~0.4x, likewise.
  auto *kk = this->lmp->kokkos;
  const int vector_length = launch_vector_length();
  const int team_size     = kk->pair_team_size_set   ? kk->pair_team_size   : 64;
  const int atoms_per_team = (team_size >= vector_length) ? team_size / vector_length : 1;
  ff.atoms_per_team = atoms_per_team;
  const int nteams = (k_list->inum + atoms_per_team - 1) / atoms_per_team;

  using policy_t = Kokkos::TeamPolicy<DeviceType, Kokkos::IndexType<int>>;
  policy_t policy(nteams, atoms_per_team, vector_length);
  Kokkos::parallel_for("PairLJCutCoulLong2::force", policy, ff);
}

/* ----------------------------------------------------------------------
   vector lanes per atom, from the standard KOKKOS package options.  The
   dual-cutoff inner list is laid out per lane, so its allocation and the
   force kernel must agree on this value.
------------------------------------------------------------------------- */

template<class DeviceType>
int PairLJCutCoulLong2Kokkos<DeviceType>::launch_vector_length() const
{
  auto *kk = this->lmp->kokkos;
  return kk->threads_per_atom_set ? kk->threads_per_atom : 8;
}

/* ----------------------------------------------------------------------
   refresh the packed j-side gather arrays.  Runs once per step, before the
   force kernel, over locals AND ghosts (a neighbor j may be either).  Reads
   are sequential and coalesced, the write is 16 B/atom: ~16 MB of streaming
   traffic per step for 500k atoms, i.e. ~15 us against a 3 ms force kernel.

   Types are repacked along with the coordinates rather than only at a
   reneighbor: it costs one byte store per atom and removes any dependence on
   when the type array was last modified.
------------------------------------------------------------------------- */

template<class DeviceType>
void PairLJCutCoulLong2Kokkos<DeviceType>::pack_refresh()
{
  const int nall = this->atom->nlocal + this->atom->nghost;

  if ((int) d_xq.extent(0) < nall) {
    const int cap = nall + nall/10 + 8;
    d_xq = Kokkos::View<LJCL2XQ*, DeviceType>();      // free before reallocating
    d_type8 = Kokkos::View<unsigned char*, DeviceType>();
    d_xq = Kokkos::View<LJCL2XQ*, DeviceType>(
      Kokkos::view_alloc("ljcl2:xq", Kokkos::WithoutInitializing), cap);
    d_type8 = Kokkos::View<unsigned char*, DeviceType>(
      Kokkos::view_alloc("ljcl2:type8", Kokkos::WithoutInitializing), cap);
  }

  LJCL2Pack<DeviceType> fp;
  fp.x = this->atomKK->k_x.template view<DeviceType>();
  fp.q = this->atomKK->k_q.template view<DeviceType>();
  fp.type = this->atomKK->k_type.template view<DeviceType>();
  fp.d_xq = d_xq;
  fp.d_type8 = d_type8;
  Kokkos::parallel_for("PairLJCutCoulLong2::pack",
                       Kokkos::RangePolicy<DeviceType>(0, nall), fp);
}

/* ----------------------------------------------------------------------
   dual-cutoff allocation.  The inner row gets the SAME capacity as the
   master row: the inner list is a subset of the master list, so the append
   in the force kernel can never overflow and needs no bounds check, no
   overflow flag and no host readback.
------------------------------------------------------------------------- */

template<class DeviceType>
void PairLJCutCoulLong2Kokkos<DeviceType>::dual_setup(bool &force_refresh)
{
  auto* k_list = static_cast<NeighListKokkos<DeviceType>*>(this->list);
  const int nall = this->atom->nlocal + this->atom->nghost;
  const int nrow = static_cast<int>(k_list->d_numneigh.extent(0));
  const int maxn = k_list->maxneighs;

  if ((int) k_ictl.view_device().extent(0) != 6) {
    k_ictl = DAT::tdual_int_1d("ljcl2:ictl", 6);
    d_maxsq = Kokkos::View<KK_FLOAT, DeviceType>("ljcl2:maxsq");
    force_refresh = true;
  }
  if ((int) d_inbr.extent(0) != nrow || (int) d_inbr.extent(1) != maxn) {
    d_inbr = typename AT::t_neighbors_2d_lr();   // free before reallocating
    d_inbr = typename AT::t_neighbors_2d_lr(
      Kokkos::view_alloc("ljcl2:inbr", Kokkos::WithoutInitializing), nrow, maxn);
    force_refresh = true;
  }
  if ((int) d_innum.extent(0) != nrow) {
    d_innum = typename AT::t_int_1d("ljcl2:innum", nrow);
    force_refresh = true;
  }
  if ((int) d_xhold.extent(0) < nall) {
    d_xhold = typename AT::t_kkfloat_1d_3_lr("ljcl2:xhold", nall + nall/10 + 8);
    force_refresh = true;
  }
  dual_nall = nall;

  // inner cutoff = the largest force cutoff over all type pairs, plus the inner
  // skin.  Taking the max (rather than a per-type inner cutoff) keeps the inner
  // list a superset of what any type pair needs.
  double cutmax = 0.0;
  const int nt = this->atom->ntypes;
  for (int i = 1; i <= nt; i++)
    for (int j = 1; j <= nt; j++)
      if (this->m_cutsq[i][j] > cutmax) cutmax = this->m_cutsq[i][j];
  const double ic = sqrt(cutmax) + dual_skin;
  inner_cutsq = ic*ic;

  // a master rebuild renumbers atoms and rebuilds the ghost shell, so the
  // displacement baseline and the inner list are both stale
  if (this->neighbor->lastcall > dual_built_step) {
    force_refresh = true;
    dual_built_step = this->neighbor->lastcall;
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairLJCutCoulLong2Kokkos<DeviceType>::dual_decide(bool force_refresh)
{
  const KK_FLOAT half_skin = static_cast<KK_FLOAT>(0.5*dual_skin);

  LJCL2DualDisp<DeviceType> fd;
  fd.x = this->atomKK->k_x.template view<DeviceType>();
  fd.xhold = d_xhold;
  Kokkos::parallel_reduce("PairLJCutCoulLong2::dual_disp",
    Kokkos::RangePolicy<DeviceType>(0, force_refresh ? 0 : dual_nall), fd,
    Kokkos::Max<KK_FLOAT, typename DeviceType::memory_space>(d_maxsq));

  LJCL2DualDecide<DeviceType> fc;
  fc.d_ictl = k_ictl.template view<DeviceType>();
  fc.d_maxsq = d_maxsq;
  fc.thresh2 = half_skin*half_skin;
  fc.every = dual_every;
  fc.force_refresh = force_refresh ? 1 : 0;
  Kokkos::parallel_for("PairLJCutCoulLong2::dual_decide",
                       Kokkos::RangePolicy<DeviceType>(0, 1), fc);

  LJCL2DualReset<DeviceType> fr;
  fr.x = this->atomKK->k_x.template view<DeviceType>();
  fr.xhold = d_xhold;
  fr.d_innum = d_innum;
  fr.d_ictl = k_ictl.template view<DeviceType>();
  fr.nlocal = this->atom->nlocal;
  Kokkos::parallel_for("PairLJCutCoulLong2::dual_reset",
                       Kokkos::RangePolicy<DeviceType>(0, dual_nall), fr);
}

/* ----------------------------------------------------------------------
   dual-cutoff verification (debug).  Fences, so it is for validation runs
   only, never for benchmarks.
------------------------------------------------------------------------- */

template<class DeviceType>
void PairLJCutCoulLong2Kokkos<DeviceType>::dual_verify()
{
  auto* k_list = static_cast<NeighListKokkos<DeviceType>*>(this->list);

  k_ictl.template modify<DeviceType>();
  k_ictl.template sync<LMPHostType>();
  k_ictl.view_host()(4) = 0;
  k_ictl.view_host()(5) = 0;
  k_ictl.template modify<LMPHostType>();
  k_ictl.template sync<DeviceType>();

  LJCL2DualCheck<DeviceType> fv;
  fv.x = this->atomKK->k_x.template view<DeviceType>();
  fv.type = this->atomKK->k_type.template view<DeviceType>();
  fv.d_ilist = k_list->d_ilist;
  fv.d_numneigh = k_list->d_numneigh;
  fv.d_neighbors = k_list->d_neighbors;
  fv.d_innum = d_innum;
  fv.d_inbr = d_inbr;
  fv.d_ictl = k_ictl.template view<DeviceType>();
  const int nt = this->atom->ntypes;
  for (int i = 1; i <= nt; i++)
    for (int j = 1; j <= nt; j++) fv.m_cutsq[i][j] = this->m_cutsq[i][j];

  Kokkos::parallel_for("PairLJCutCoulLong2::dual_verify",
                       Kokkos::RangePolicy<DeviceType>(0, k_list->inum), fv);

  k_ictl.template modify<DeviceType>();
  k_ictl.template sync<LMPHostType>();
  const int master = k_ictl.view_host()(4);
  const int inner  = k_ictl.view_host()(5);
  if (this->comm->me == 0)
    printf("DUAL_CHECK: step " BIGINT_FORMAT " refresh %d in-cutoff pairs master %d "
           "inner %d -> %s\n", this->update->ntimestep, k_ictl.view_host()(0),
           master, inner, (master == inner) ? "OK" : "*** MISSED PAIRS ***");
  if (master != inner)
    this->error->warning(FLERR, "dual-cutoff inner list is missing in-cutoff pairs");
}

/* ----------------------------------------------------------------------
   per-atom views, constants and stack parameters shared by both kernels
------------------------------------------------------------------------- */

template<class DeviceType>
template<class FunctorT>
void PairLJCutCoulLong2Kokkos<DeviceType>::fill_common(FunctorT &ff)
{
  auto* k_list = static_cast<NeighListKokkos<DeviceType>*>(this->list);

  ff.x    = this->atomKK->k_x.template view<DeviceType>();
  ff.f    = this->atomKK->k_f.template view<DeviceType>();
  ff.q    = this->atomKK->k_q.template view<DeviceType>();
  ff.type = this->atomKK->k_type.template view<DeviceType>();
  ff.d_ilist = k_list->d_ilist;
  ff.inum    = k_list->inum;

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
}

/* ----------------------------------------------------------------------
   merge the flat half list into the i-batched union list.  Runs once per
   reneighbor.  The per-group capacity is self-tuning: the kernel always
   reports the capacity it needed (whether or not it fit), so one retry
   converges, and no capacity guess is baked in.
------------------------------------------------------------------------- */

template<class DeviceType>
template<int CI>
void PairLJCutCoulLong2Kokkos<DeviceType>::union_build()
{
  auto* k_list = static_cast<NeighListKokkos<DeviceType>*>(this->list);
  const int inum = k_list->inum;
  const int ngroup = (inum + CI - 1) / CI;
  const int nslot = (CI > 4) ? 8192 : 4096;   // power of 2, >~2x the entry count
  union_ngroup = ngroup;
  if (union_maxu < 64) union_maxu = 64;

  if ((int) k_ustat.view_device().extent(0) != 5)
    k_ustat = DAT::tdual_int_1d("ljcl2:ustat", 5);

  using policy_t = Kokkos::TeamPolicy<DeviceType, Kokkos::IndexType<int>>;
  typedef Kokkos::View<int*, typename DeviceType::scratch_memory_space,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>> t_sh_int_1d;
  const int build_team = 128;

  int attempt = 0;
  for (; attempt < 5; attempt++) {
    if ((int) d_unbr.extent(0) != ngroup || (int) d_unbr.extent(1) != union_maxu) {
      d_unbr = typename AT::t_neighbors_2d_lr();   // free before reallocating
      d_unbr = typename AT::t_neighbors_2d_lr(
        Kokkos::view_alloc("ljcl2:unbr", Kokkos::WithoutInitializing), ngroup, union_maxu);
    }
    if ((int) d_unum.extent(0) != ngroup)
      d_unum = typename AT::t_int_1d("ljcl2:unum", ngroup);

    for (int k = 0; k < 5; k++) k_ustat.view_host()(k) = 0;
    k_ustat.template modify<LMPHostType>();
    k_ustat.template sync<DeviceType>();

    LJCL2UnionBuild<DeviceType,CI> ub;
    ub.d_ilist    = k_list->d_ilist;
    ub.d_numneigh = k_list->d_numneigh;
    ub.d_neighbors = k_list->d_neighbors;
    ub.d_neighbors_t = k_list->d_neighbors_transpose;
    ub.use_transpose = this->lmp->kokkos->neigh_transpose;
    ub.d_unbr = d_unbr;
    ub.d_unum = d_unum;
    ub.d_ustat = k_ustat.template view<DeviceType>();
    ub.inum  = inum;
    ub.maxu  = union_maxu;
    ub.nslot = nslot;

    policy_t policy(ngroup, build_team, 1);
    policy.set_scratch_size(0, Kokkos::PerTeam(t_sh_int_1d::shmem_size(nslot) +
                                              t_sh_int_1d::shmem_size(1)));
    Kokkos::parallel_for("PairLJCutCoulLong2::union_build", policy, ub);

    k_ustat.template modify<DeviceType>();
    k_ustat.template sync<LMPHostType>();

    if (!k_ustat.view_host()(0)) break;
    union_maxu = static_cast<int>(1.15 * k_ustat.view_host()(1)) + 16;
  }
  if (attempt == 5)
    this->error->one(FLERR, "union neighbor list capacity did not converge");

  if (union_stats && this->comm->me == 0) {
    const double entries = k_ustat.view_host()(2);
    const double pairs   = k_ustat.view_host()(3);
    printf("UNION_STATS: CI %d groups %d entries %.0f pairs %.0f "
           "entries/atom %.1f share %.2f specials %.0f max_need %d maxu %d retries %d\n",
           CI, ngroup, entries, pairs, entries/inum, pairs/entries,
           (double) k_ustat.view_host()(4), k_ustat.view_host()(1), union_maxu, attempt);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int CI>
void PairLJCutCoulLong2Kokkos<DeviceType>::union_compute()
{
  LJCL2UnionForce<DeviceType,CI> ff;
  fill_common(ff);
  ff.d_unbr  = d_unbr;
  ff.d_unum  = d_unum;
  ff.ngroup  = union_ngroup;

  // launch shape from the same package options as the flat path, but the unit
  // of work is a GROUP of CI atoms: pair/team/size / threads/per/atom groups
  // share a team, each group's entry list split over threads/per/atom lanes.
  auto *kk = this->lmp->kokkos;
  const int vector_length = kk->threads_per_atom_set ? kk->threads_per_atom : 8;
  const int team_size     = kk->pair_team_size_set   ? kk->pair_team_size   : 64;
  const int groups_per_team = (team_size >= vector_length) ? team_size / vector_length : 1;
  ff.groups_per_team = groups_per_team;
  const int nteams = (union_ngroup + groups_per_team - 1) / groups_per_team;

  using policy_t = Kokkos::TeamPolicy<DeviceType, Kokkos::IndexType<int>>;
  policy_t policy(nteams, groups_per_team, vector_length);
  Kokkos::parallel_for("PairLJCutCoulLong2::force_union", policy, ff);
}

namespace LAMMPS_NS {
template class PairLJCutCoulLong2Kokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairLJCutCoulLong2Kokkos<LMPHostType>;
#endif
}
