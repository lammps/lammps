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
#include "math_const.h"
#include "neigh_list_kokkos.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "update.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <type_traits>
#include <vector>

using namespace LAMMPS_NS;
using namespace EwaldConst;
using MathConst::MY_PI;

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
  LJCL2_VAR_FASTMATH = 4,  // normal, with __expf / __fdividef for the Ewald block
  LJCL2_VAR_RATIONAL = 5   // normal, with the rational-polynomial Ewald force (below)
};

// ---------------------------------------------------------------------------
// Rational-polynomial Ewald real-space force (LMP_LJCL2_VARIANT=5, and the
// compact kernel's RAT option).  Borrowed from GROMACS (pmecorrF in the nbnxm
// kernels): with z2 = (g_ewald r)^2, the erfc-screened Coulomb force is
//
//   F/r = qqrd2e qi qj [ 1/r^3 + g_ewald^3 P(z2) ],
//   P(z2) = -( erf(z)/z^3 - 2 exp(-z2)/(sqrt(pi) z^2) )  ~  N(z2)/D(z2)
//
// a 6th/4th-order rational fit with no exp, no erfc and one reciprocal.  The
// special-bond correction of the analytical path, -(1-factor_coul) qq/r^3,
// folds into the first term as factor_coul/r^3, so the branch disappears too.
//
// Accuracy, evaluated in fp32 against the exact erfc force over r in [1,12]
// at this input's g_ewald (0.1678): max |dF| per unit qq 2.5e-7, against
// 5.4e-7 for the Abramowitz-Stegun path in fp32.  The relative error near the
// cutoff is larger (1.4e-5 vs 3.6e-6) because 1/r^3 and the correction
// nearly cancel there, but the force itself is ~3e-4 of its r=1 value.  The
// fit is a single-precision one, so a double build keeps the A&S path.
// ---------------------------------------------------------------------------
KOKKOS_INLINE_FUNCTION float ljcl2_pmecorrF(const float z2) {
  const float z4 = z2*z2;
  float d0 = 0.0011193462567257629232f*z4 + 0.11583842382862377919f;
  const float d1 = 0.014866955030185295499f*z4 + 0.50736591960530292870f;
  d0 = d0*z4 + 1.0f;
  d0 = d1*z2 + d0;
  float n0 = -1.7357322914161492954e-8f*z4 - 0.000053401640219807709149f;
  float n1 = 1.4703624142580877519e-6f*z4 + 0.0010054721316683106153f;
  n0 = n0*z4 - 0.019278317264888380590f;
  n1 = n1*z4 + 0.069670166153766424023f;
  n0 = n0*z4 - 0.75225204789749321333f;
  n0 = n1*z2 + n0;
#if defined(__CUDA_ARCH__)
  return __fdividef(n0, d0);   // d0 >= 1 and modest (6.1 at z = 2): 2 ulp, no slow path
#else
  return n0 / d0;
#endif
}

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
  template<bool FAST = false, bool NOCOUL = false, bool RAT = false>
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT pair_fpair(const KK_FLOAT rsq, const int itype, const int jtype,
                      const KK_FLOAT qiqj, const int sb) const {
    return pair_fpair_p<FAST,NOCOUL,RAT>(rsq, qiqj, sb,
             m_cut_ljsq[itype][jtype], m_cut_coulsq[itype][jtype],
             m_params[itype][jtype].lj1, m_params[itype][jtype].lj2);
  }

  // same, with the per-type-pair parameters supplied by the caller (the compact
  // kernel reads them from a shared-memory table instead of constant memory)
  template<bool FAST = false, bool NOCOUL = false, bool RAT = false>
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT pair_fpair_p(const KK_FLOAT rsq, const KK_FLOAT qiqj, const int sb,
                        const KK_FLOAT cut_ljsq, const KK_FLOAT cut_coulsq,
                        const KK_FLOAT lj1, const KK_FLOAT lj2) const {
    // one reciprocal-sqrt supplies rinv, r2inv, and r (was 2 divides + 1 sqrt)
    const KK_FLOAT rinv = Kokkos::rsqrt(rsq);
    const KK_FLOAT r2inv = rinv*rinv;
    KK_FLOAT fpair = static_cast<KK_FLOAT>(0.0);

    if (rsq < cut_ljsq) {
      const KK_FLOAT r6inv = r2inv*r2inv*r2inv;
      fpair += special_lj[sb] * r6inv*(lj1*r6inv - lj2) * r2inv;
    }

    if (RAT && !NOCOUL && (sizeof(KK_FLOAT) == sizeof(float)) && rsq < cut_coulsq) {
      const KK_FLOAT g2 = g_ewald*g_ewald;
      const KK_FLOAT pc = static_cast<KK_FLOAT>(ljcl2_pmecorrF(static_cast<float>(g2*rsq)));
      fpair += qqrd2e * qiqj * (special_coul[sb]*rinv*r2inv + g2*g_ewald*pc);
    } else if (!NOCOUL && rsq < cut_coulsq) {
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

template<class DeviceType, bool FULL, bool DUAL, bool PACK, int VAR = LJCL2_VAR_NORMAL,
         bool BALLOT = false>
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
  // The append is a contended per-atom atomic, and it is expensive: a refresh
  // step measures 5.072 ms against 3.000 for a plain flat step, and the master
  // walk explains only the 3.000 (kokkos_neigh.md 10.7).  With BALLOT it is
  // replaced by a warp-aggregated append -- see the header.  An earlier per-lane
  // SEGMENTED variant (lane L writes keepers at L + vl*t) was built and measured
  // 10% SLOWER: splitting the loop into a refresh branch and a tight branch
  // inlines this body twice (REG 43 -> 55) and gives the tight loop a per-lane
  // trip count that nvcc can no longer unroll, which costs more than the atomics
  // saved.  BALLOT avoids that trap: one body, one branch, and the inner row
  // stays compact.
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

    if (DUAL && store) {
      const bool keep = (rsq < inner_cutsq);
      if (BALLOT) {
#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ >= 700)
        // Warp-aggregated append.  __match_any_sync isolates the lanes working
        // on THIS atom, so nothing is assumed about how Kokkos maps vector lanes
        // onto the warp, and __ballot_sync gives that group's keeper mask.  One
        // lane reserves the group's whole block with a single atomic and
        // broadcasts the base, so the atomics per row drop from one per keeper
        // (~360 at inner skin 2.5) to one per iteration that keeps anything
        // (~95), and the lanes of a group stop serializing on one address.
        //
        // Correctness must NOT depend on which lanes __activemask() reports --
        // it carries no convergence guarantee, and this is why the base comes
        // from an atomic instead of from a count carried in registers.  Were the
        // group ever split across two masks, a register-carried count would be
        // short in both halves and the halves would overwrite each other,
        // silently dropping pairs from the inner list.  With a per-subset
        // reservation any partition of the group stays correct; only the order
        // within the row changes, and the row is a set (the special-bond bits
        // ride inside jraw), so order is immaterial.  A truly atomic-free form
        // would need every lane of the group to reach the vote on every
        // iteration, i.e. padding the vector range to a multiple of the vector
        // length and clamping the index -- worth trying only if the measurement
        // says these remaining atomics still cost something.
        const unsigned active = __activemask();
        const unsigned gmask  = __match_any_sync(active, i);
        const unsigned pred   = __ballot_sync(active, keep) & gmask;
        if (pred) {   // uniform across gmask: every lane of it computed `pred`
          unsigned laneid;
          asm volatile("mov.u32 %0, %%laneid;" : "=r"(laneid));
          const int leader = __ffs(pred) - 1;
          int base = 0;
          if ((int) laneid == leader)
            base = Kokkos::atomic_fetch_add(&d_innum(i), __popc(pred));
          base = __shfl_sync(gmask, base, leader);
          if (keep) d_inbr(i, base + __popc(pred & ((1u << laneid) - 1u))) = jraw;
        }
#else
        // Host, pre-Volta and non-CUDA backends have no warp vote, so fall back
        // to the plain per-keeper atomic.  A per-lane running count would be
        // WRONG here whenever vector_length > 1: every lane of the group would
        // write at its own offset and keepers would overwrite each other.
        if (keep) d_inbr(i, Kokkos::atomic_fetch_add(&d_innum(i), 1)) = jraw;
#endif
      } else {
        if (keep) d_inbr(i, Kokkos::atomic_fetch_add(&d_innum(i), 1)) = jraw;
      }
    }

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
        this->template pair_fpair<VAR == LJCL2_VAR_FASTMATH, VAR == LJCL2_VAR_NOCOUL,
                                  VAR == LJCL2_VAR_RATIONAL>(
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
   Warp-compacted force kernel (LMP_LJCL2_COMPACT=1), CUDA only.

   The idea is the "wavefront"/stream-compaction trick of GPU ray tracers.
   In the flat kernel ~42% of the entries of this input pass the cutoff test,
   which means that essentially EVERY warp iteration contains at least one
   passing lane, so the expensive body (rsqrt, LJ, Ewald, three atomics) is
   issued on every iteration with ~42% of its lanes doing useful work.

   Here one warp owns one atom i at a time and splits it in two phases:
     1. test phase: 32 lanes take 32 consecutive entries, gather x(j),
        type(j) and test the cutoff.  Passing pairs are appended to a
        per-warp shared-memory queue (ballot + popc prefix, no atomics).
     2. physics phase: whenever the queue holds >= 32 pairs, all 32 lanes
        pop one each and run the body fully converged.  The tail (< 32
        pairs) is drained once per atom.
   The queue holds (delx, dely, delz, jraw); type(j) and q(j) are re-read in
   the physics phase (L1 hits: the line was touched in phase 1).

   The per-type-pair parameters can come from a shared-memory table (SMEM),
   loaded once per team from a small device View.  In the flat kernel they
   come from the functor in constant memory, where 32 lanes indexing with 32
   different (itype,jtype) serialize into up to 32 constant-cache replays per
   load; the GPU package keeps them in shared memory for the same reason.
   Each warp walks atoms_per_warp atoms so that the table fill is amortized.

   One atom per warp is also what makes the compaction cheap: every queued
   pair belongs to the same i, so the i-force needs no per-owner bookkeeping.

   Measured on melt2 (A100, single, neigh/transpose on): force kernel 4.56 ->
   3.97 ms with constant-memory tables, 4.13 ms with SMEM (the table costs ~12
   registers and the constant-cache replays turned out to be cheap), 3.85 ms
   with RAT as well.  Without neigh/transpose the 32 lanes read 32 different
   rows of the LayoutLeft list and the kernel is 3x slower.

   PACK takes x, q and type from the packed (x,y,z,q) records and carries q(j)
   and type(j) in the queue, so the physics phase issues no gathers at all.
   DUAL combines the kernel with the dual-cutoff inner list: since one warp
   owns a whole row, the refresh-step append is a ballot prefix plus a
   register counter -- none of the contended atomics of the flat-kernel
   append (see the BALLOT notes in the header).
------------------------------------------------------------------------- */

template<class DeviceType, bool FULL, bool RAT, bool SMEM, bool PACK = false, bool DUAL = false>
struct LJCL2Compact : public LJCL2Common<DeviceType> {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef Kokkos::TeamPolicy<DeviceType, Kokkos::IndexType<int>> policy_type;
  typedef typename policy_type::member_type member_type;

  static constexpr int QCAP = 64;   // max queue fill is 31 + 32 = 63
  enum { T_CUTSQ = 0, T_CUTLJSQ = 1, T_CUTCOULSQ = 2, T_LJ1 = 3, T_LJ2 = 4, T_N = 5 };
  // compile-time stride between the table sections, so the section offsets
  // are LDS immediates rather than registers
  static constexpr int TSTRIDE = (MAX_TYPES_STACKPARAMS+1)*(MAX_TYPES_STACKPARAMS+1);
  static constexpr int TABLEN = SMEM ? T_N*TSTRIDE : 0;
  // queue record per slot: del x/y/z (+ q(j) with PACK), jraw (+ type(j) with PACK)
  static constexpr int QF = PACK ? 4 : 3;
  static constexpr int QI = PACK ? 2 : 1;

  typename AT::t_int_1d d_numneigh;
  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_neighbors_2d_lr d_neighbors_t;
  Kokkos::View<const KK_FLOAT*, DeviceType> d_tab;   // [T_N][TSTRIDE], index i*ntp1+j
  Kokkos::View<const LJCL2XQ*, DeviceType> d_xq;     // PACK only
  Kokkos::View<const unsigned char*, DeviceType> d_type8;
  typename AT::t_neighbors_2d_lr d_inbr;             // DUAL only
  typename AT::t_int_1d d_innum;
  typename AT::t_int_1d d_ictl;
  KK_FLOAT inner_cutsq;
  int ntp1, warps_per_team, atoms_per_warp;
  bool use_transpose;

  KOKKOS_INLINE_FUNCTION
  static int shmem_bytes(int warps) {
    return TABLEN*sizeof(KK_FLOAT) + warps*QCAP*(QF*sizeof(KK_FLOAT) + QI*sizeof(int)) + 64;
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const member_type& team) const {
#if defined(__CUDA_ARCH__)
    const int nwarp = warps_per_team;
    char *sh = static_cast<char *>(team.team_shmem().get_shmem(shmem_bytes(nwarp)));
    KK_FLOAT *tab = reinterpret_cast<KK_FLOAT *>(sh);
    const int w = team.team_rank();
    const int lane = threadIdx.x;
    KK_FLOAT *qx = tab + TABLEN + w*QF*QCAP;
    KK_FLOAT *qy = qx + QCAP;
    KK_FLOAT *qz = qy + QCAP;
    KK_FLOAT *qq = qz + QCAP;                            // PACK only
    int *qj = reinterpret_cast<int *>(tab + TABLEN + nwarp*QF*QCAP) + w*QI*QCAP;
    int *qt = qj + QCAP;                                 // PACK only

    if (SMEM) {
      for (int k = w*32 + lane; k < T_N*TSTRIDE; k += nwarp*32) tab[k] = d_tab(k);
      team.team_barrier();
    }

    const unsigned FULLMASK = 0xffffffffu;
    const unsigned lanelt = (1u << lane) - 1u;
    const int a0 = (team.league_rank()*nwarp + w)*atoms_per_warp;

    // dual cutoff: a refresh step walks the master list and rebuilds the inner
    // row; the other steps walk the inner row.  d_ictl(0) is launch-uniform.
    const int refresh = DUAL ? d_ictl(0) : 1;
    const bool walk_inner = DUAL && !refresh;
    const bool store = DUAL && refresh;

    for (int a = 0; a < atoms_per_warp; a++) {
      const int ii = a0 + a;
      if (ii >= this->inum) break;                       // warp-uniform
      const int i = this->d_ilist(ii);
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
      const int ibase = itype*ntp1;
      const int jnum = walk_inner ? d_innum(i) : d_numneigh(i);

      KK_ACC_FLOAT fx = 0.0, fy = 0.0, fz = 0.0;
      int qn = 0;                                        // warp-uniform
      int icnt = 0;                                      // inner-row length (refresh)

      // physics for queue slot p (all lanes of the warp, or the tail)
      auto physics = [&] (const int p) {
        const KK_FLOAT delx = qx[p], dely = qy[p], delz = qz[p];
        const int jraw = qj[p];
        const int j = jraw & NEIGHMASK;
        const int sb = (jraw >> SBBITS) & 3;
        const int jtype = PACK ? qt[p] : this->type(j);
        const KK_FLOAT qj_ = PACK ? qq[p] : this->q(j);
        const KK_FLOAT rsq = delx*delx + dely*dely + delz*delz;
        KK_FLOAT fpair;
        if (SMEM) {
          const int t = ibase + jtype;
          fpair = this->template pair_fpair_p<false,false,RAT>(rsq, qtmp*qj_, sb,
                    tab[T_CUTLJSQ*TSTRIDE + t], tab[T_CUTCOULSQ*TSTRIDE + t],
                    tab[T_LJ1*TSTRIDE + t], tab[T_LJ2*TSTRIDE + t]);
        } else {
          fpair = this->template pair_fpair<false,false,RAT>(rsq, itype, jtype, qtmp*qj_, sb);
        }
        fx += static_cast<KK_ACC_FLOAT>(delx*fpair);
        fy += static_cast<KK_ACC_FLOAT>(dely*fpair);
        fz += static_cast<KK_ACC_FLOAT>(delz*fpair);
        if (!FULL && (this->newton_pair || j < this->nlocal)) {
          Kokkos::atomic_add(&this->f(j,0), -static_cast<KK_ACC_FLOAT>(delx*fpair));
          Kokkos::atomic_add(&this->f(j,1), -static_cast<KK_ACC_FLOAT>(dely*fpair));
          Kokkos::atomic_add(&this->f(j,2), -static_cast<KK_ACC_FLOAT>(delz*fpair));
        }
      };

      // One physics call site serves both the full pops and the tail drain:
      // inlining the body twice costs ~15 registers.  After any pop qn <= 31,
      // so on the final (base >= jnum) pass the queue holds only the tail.
      for (int base = 0; ; base += 32) {
        const bool more = (base < jnum);                 // warp-uniform
        if (more) {
          const int jj = base + lane;
          bool pass = false, keep = false;
          KK_FLOAT delx = 0, dely = 0, delz = 0, qj_ = 0;
          int jraw = 0, jtype = 0;
          if (jj < jnum) {
            jraw = walk_inner ? d_inbr(i,jj)
                 : (use_transpose ? d_neighbors_t(i,jj) : d_neighbors(i,jj));
            const int j = jraw & NEIGHMASK;
            if (PACK) {
              const LJCL2XQ j4 = d_xq(j);
              delx = xtmp - j4.x; dely = ytmp - j4.y; delz = ztmp - j4.z; qj_ = j4.q;
              jtype = static_cast<int>(d_type8(j));
            } else {
              delx = xtmp - this->x(j,0);
              dely = ytmp - this->x(j,1);
              delz = ztmp - this->x(j,2);
              jtype = this->type(j);
            }
            const KK_FLOAT rsq = delx*delx + dely*dely + delz*delz;
            const KK_FLOAT cutsq = SMEM ? tab[T_CUTSQ*TSTRIDE + ibase + jtype]
                                        : this->m_cutsq[itype][jtype];
            pass = (rsq < cutsq);
            if (store) keep = (rsq < inner_cutsq);
          }
          // Inner-row rebuild.  The warp owns row i outright, so the append is
          // a ballot prefix and a register counter: no atomics at all, and the
          // row is written in order, 32 contiguous ints at a time.
          if (store) {
            const unsigned kb = __ballot_sync(FULLMASK, keep);
            if (keep) d_inbr(i, icnt + __popc(kb & lanelt)) = jraw;
            icnt += __popc(kb);
          }
          const unsigned b = __ballot_sync(FULLMASK, pass);
          if (pass) {
            const int p = qn + __popc(b & lanelt);
            qx[p] = delx; qy[p] = dely; qz[p] = delz; qj[p] = jraw;
            if (PACK) { qq[p] = qj_; qt[p] = jtype; }
          }
          qn += __popc(b);
        }
        if (qn >= 32 || (!more && qn > 0)) {
          __syncwarp();
          const int take = (qn >= 32) ? 32 : qn;
          if (lane < take) physics(qn - take + lane);    // pop the top `take` (LIFO)
          qn -= take;
          __syncwarp();                                  // pops done before next pushes
        }
        if (!more) break;
      }

      for (int o = 16; o > 0; o >>= 1) {
        fx += __shfl_xor_sync(FULLMASK, fx, o);
        fy += __shfl_xor_sync(FULLMASK, fy, o);
        fz += __shfl_xor_sync(FULLMASK, fz, o);
      }
      if (lane == 0) {
        if (store) d_innum(i) = icnt;
        if (FULL) {
          this->f(i,0) = fx; this->f(i,1) = fy; this->f(i,2) = fz;
        } else {
          Kokkos::atomic_add(&this->f(i,0), fx);
          Kokkos::atomic_add(&this->f(i,1), fy);
          Kokkos::atomic_add(&this->f(i,2), fz);
        }
      }
    }
#else
    (void) team;
    Kokkos::abort("lj/cut/coul/long2/kk: the compact kernel requires the CUDA backend");
#endif
  }
};

/* ----------------------------------------------------------------------
   Shared-memory tiled kernel (LMP_LJCL2_LOCAL=1|2), CUDA + single only.

   Borrowed from histogram privatization and tiled SpMV/SpMM rather than
   from MD codes.  A group of NI consecutive (spatially sorted) i atoms has
   a compact j-SET: the union of their neighbor rows.  At reneighbor time a
   build kernel merges the rows in a shared-memory hash table, emits the
   j-set, and rewrites every row entry as a 16-bit LOCAL index into it
   (14 bits index + 2 special-bond bits).  The force kernel then runs one
   group per block:
     1. gather the j-set (x,y,z,q,type) into shared memory, once;
     2. warp-per-atom compaction loop (as in LJCL2Compact), but every j
        coordinate read is an LDS and the index stream is 2 bytes/entry;
     3. j-forces accumulate in SHARED memory, then one flush per distinct j
        replaces one global atomic triple per pair.
   With ~10 pairs per distinct j, the global gathers and atomics drop ~10x.

   A100 has no native shared-memory float atomic add (it compiles to a CAS
   loop, ATOMS.CAST.SPIN), and none for 64-bit integers either, but int32 is
   native.  LOCAL=1 uses the float CAS; LOCAL=2 uses exact fixed point,
   v = llrint(f * 2^24) split into hi = v >> 16 and lo = v & 0xffff, summed
   with two native int32 atomics each (lo cannot overflow for < 32768
   contributions, and there are at most NI per j).

   The kernels are raw CUDA because Kokkos caps team scratch at the 48 KB
   default; a 64-atom tile needs 60-110 KB (A100 opts in up to 163 KB).
------------------------------------------------------------------------- */

#if defined(KOKKOS_ENABLE_CUDA) && defined(__CUDACC__)
#define LJCL2_HAVE_LOCAL 1

enum { LJCL2_LMAX = 16383, LJCL2_LSB_SHIFT = 14, LJCL2_NTAB = 13*13 };

struct LJCL2LocalArgs {
  const float4 *xq;
  const unsigned char *type8;
  const int *ilist, *numneigh;
  const unsigned short *lnbr;
  long long lnbr_ld;
  const int *ljset, *lm;
  unsigned short *linr;        // DUAL: inner local rows (same shape as lnbr)
  int *linnum;
  const int *ictl;             // DUAL: [0] = refresh flag (device-side)
  float inner_cutsq;
  float *f;
  long long fs0, fs1;
  int jcap, inum, ni, ms, ntp1;
  float cutsq_max, g_ewald, qqrd2e;
  float special_lj[4], special_coul[4];
  float lj1[LJCL2_NTAB], lj2[LJCL2_NTAB], cutlj[LJCL2_NTAB], cutcoul[LJCL2_NTAB];
};

struct LJCL2LocalBuildArgs {
  const int *ilist, *numneigh, *nbr;
  long long ns0, ns1;          // nbr(i,jj) = nbr[i*ns0 + jj*ns1]
  unsigned short *lnbr;
  long long lnbr_ld;
  int *ljset, *lm, *stat;      // stat: [0] max M, [1] overflow, [2] sum M, [3] sum entries
  int jcap, inum, ni, nslot;
};

// per-warp pair queue: one test pass pushes up to 32*VEC entries on top of
// at most 31 left over, so 160 covers VEC <= 4
enum { LJCL2_LQCAP = 160 };

// bytes of dynamic shared memory of the force kernel
inline int ljcl2_local_shmem(int ms, int nwarp, int fix) {
  return ms*16 + ms*(fix ? 24 : 12) + ((ms + 15) & ~15) + nwarp*LJCL2_LQCAP*4 + nwarp*64*4;
}

__global__ void __launch_bounds__(256) ljcl2_local_build(const LJCL2LocalBuildArgs a)
{
  extern __shared__ int shb[];
  const int nslot = a.nslot, mask = nslot - 1;
  int *keys = shb;              // j+1, 0 = empty
  int *rank = shb + nslot;
  __shared__ int s_wsum[32];
  __shared__ int s_ovf;
  const int tid = threadIdx.x, nthr = blockDim.x;
  const int g = blockIdx.x, p0 = g*a.ni;
  const int nrow = (a.inum - p0 < a.ni) ? (a.inum - p0) : a.ni;

  for (int s = tid; s < nslot; s += nthr) keys[s] = 0;
  if (tid == 0) s_ovf = 0;
  __syncthreads();

  int nent = 0;
  for (int r = 0; r < nrow; r++) {
    const int i = a.ilist[p0 + r];
    const int jnum = a.numneigh[i];
    nent += (tid == 0) ? jnum : 0;
    for (int jj = tid; jj < jnum; jj += nthr) {
      const int j = a.nbr[i*a.ns0 + jj*a.ns1] & NEIGHMASK;
      int s = j & mask;
      for (int p = 0; p < nslot; p++) {
        const int old = atomicCAS(&keys[s], 0, j+1);
        if (old == 0 || old == j+1) break;
        s = (s + 1) & mask;
        if (p == nslot-1) s_ovf = 1;
      }
    }
  }
  __syncthreads();

  // block-wide exclusive scan over occupied slots, contiguous chunk per thread
  // so that the local index order follows the slot order (= runs of j)
  const int chunk = nslot / nthr;
  const int s0 = tid*chunk;
  int cnt = 0;
  for (int s = s0; s < s0 + chunk; s++) cnt += (keys[s] != 0);
  const int lane = tid & 31, w = tid >> 5;
  int incl = cnt;
  for (int o = 1; o < 32; o <<= 1) {
    const int v = __shfl_up_sync(0xffffffffu, incl, o);
    if (lane >= o) incl += v;
  }
  if (lane == 31) s_wsum[w] = incl;
  __syncthreads();
  if (w == 0) {
    int v = (lane < (nthr >> 5)) ? s_wsum[lane] : 0;
    for (int o = 1; o < 32; o <<= 1) {
      const int u = __shfl_up_sync(0xffffffffu, v, o);
      if (lane >= o) v += u;
    }
    s_wsum[lane] = v;   // inclusive warp totals
  }
  __syncthreads();
  const int total = s_wsum[(nthr >> 5) - 1];
  int base = incl - cnt + (w ? s_wsum[w-1] : 0);
  const bool fits = (total <= a.jcap) && (total <= LJCL2_LMAX) && !s_ovf;
  int *js = a.ljset + (long long) g*a.jcap;
  for (int s = s0; s < s0 + chunk; s++) {
    const int k = keys[s];
    if (k) { rank[s] = base; if (fits) js[base] = k - 1; base++; }
  }
  if (tid == 0) {
    a.lm[g] = fits ? total : 0;
    atomicMax(&a.stat[0], total);
    if (!fits) atomicOr(&a.stat[1], 1);
    atomicAdd(&a.stat[2], total);
    atomicAdd(&a.stat[3], nent);
  }
  if (!fits) return;   // uniform: total and s_ovf are block-wide
  __syncthreads();

  for (int r = 0; r < nrow; r++) {
    const int i = a.ilist[p0 + r];
    const int jnum = a.numneigh[i];
    unsigned short *row = a.lnbr + (long long) i*a.lnbr_ld;
    for (int jj = tid; jj < jnum; jj += nthr) {
      const int jraw = a.nbr[i*a.ns0 + jj*a.ns1];
      const int j = jraw & NEIGHMASK;
      const int sb = (jraw >> SBBITS) & 3;
      int s = j & mask;
      while (keys[s] != j+1) s = (s + 1) & mask;
      row[jj] = static_cast<unsigned short>(rank[s] | (sb << LJCL2_LSB_SHIFT));
    }
  }
}

// VEC = row entries per lane per test pass, read with one 2/4/8-byte load.
// DUAL: a refresh step walks the master local row and rebuilds the inner row
// (pairs within cut + inner skin) with a ballot prefix -- the warp owns the
// row, so no atomics; the other steps walk the inner row only.
template<bool RAT, int FIX, int VEC, bool DUAL>
__global__ void __launch_bounds__(256, 2) ljcl2_local_force(const LJCL2LocalArgs a)
{
  extern __shared__ float4 sh4[];
  const int MS = a.ms;
  const int nwarp = blockDim.x >> 5;
  const int w = threadIdx.x >> 5, lane = threadIdx.x & 31;
  float4 *s_xq = sh4;
  float *s_f = reinterpret_cast<float *>(s_xq + MS);            // FIX==0: [3][MS] float
  int *s_fi = reinterpret_cast<int *>(s_xq + MS);               // FIX!=0: [6][MS] int
  unsigned char *s_t = reinterpret_cast<unsigned char *>(s_f + (FIX ? 6 : 3)*MS);
  int *s_q = reinterpret_cast<int *>(s_t + ((MS + 15) & ~15)) + w*LJCL2_LQCAP;
  float *s_row = reinterpret_cast<float *>(s_t + ((MS + 15) & ~15) + nwarp*LJCL2_LQCAP*4) + w*64;

  const int g = blockIdx.x;
  const int M = a.lm[g];
  const int *js = a.ljset + (long long) g*a.jcap;
  for (int k = threadIdx.x; k < M; k += blockDim.x) {
    const int j = js[k];
    s_xq[k] = a.xq[j];
    s_t[k] = a.type8[j];
    if (FIX) {
      for (int c = 0; c < 6; c++) s_fi[c*MS + k] = 0;
    } else {
      s_f[k] = 0.0f; s_f[MS + k] = 0.0f; s_f[2*MS + k] = 0.0f;
    }
  }
  __syncthreads();

  const unsigned FULLMASK = 0xffffffffu;
  const unsigned lanelt = (1u << lane) - 1u;
  const float g2 = a.g_ewald*a.g_ewald, g3 = g2*a.g_ewald;
  const int p0 = g*a.ni;
  const int refresh = DUAL ? a.ictl[0] : 1;                     // launch-uniform
  const bool walk_inner = DUAL && !refresh;
  const bool store = DUAL && refresh;

  for (int r = w; r < a.ni; r += nwarp) {
    const int ii = p0 + r;
    if (ii >= a.inum) break;                                   // warp-uniform
    const int i = a.ilist[ii];
    const float4 xi = a.xq[i];
    const int itype = a.type8[i];
    __syncwarp();
    if (lane < a.ntp1) {
      const int t = itype*a.ntp1 + lane;
      s_row[lane] = a.lj1[t]; s_row[16+lane] = a.lj2[t];
      s_row[32+lane] = a.cutlj[t]; s_row[48+lane] = a.cutcoul[t];
    }
    __syncwarp();
    const int jnum = walk_inner ? a.linnum[i] : a.numneigh[i];
    const unsigned short *row = (walk_inner ? a.linr : a.lnbr) + (long long) i*a.lnbr_ld;
    unsigned short *irow = a.linr + (long long) i*a.lnbr_ld;
    float fx = 0.0f, fy = 0.0f, fz = 0.0f;
    int qn = 0, icnt = 0;

    for (int base = 0; ; base += 32*VEC) {
      const bool more = (base < jnum);
      if (more) {
        // lane owns VEC consecutive entries: one vector load (the row stride
        // is padded to a multiple of 8 entries, so it never crosses a row)
        const int j0 = base + lane*VEC;
        int ev[VEC];
#pragma unroll
        for (int v = 0; v < VEC; v++) ev[v] = 0;
        if (j0 < jnum) {
          if (VEC == 1) {
            ev[0] = row[j0];
          } else if (VEC == 2) {
            const unsigned u = *reinterpret_cast<const unsigned *>(row + j0);
            ev[0] = u & 0xffff; ev[VEC-1] = u >> 16;
          } else {
            const uint2 u = *reinterpret_cast<const uint2 *>(row + j0);
            ev[0] = u.x & 0xffff; ev[1 % VEC] = u.x >> 16;
            ev[2 % VEC] = u.y & 0xffff; ev[3 % VEC] = u.y >> 16;
          }
        }
#pragma unroll
        for (int v = 0; v < VEC; v++) {
          bool pass = false, keep = false;
          if (j0 + v < jnum) {
            const float4 p = s_xq[ev[v] & LJCL2_LMAX];
            const float dx = xi.x - p.x, dy = xi.y - p.y, dz = xi.z - p.z;
            const float rsq = dx*dx + dy*dy + dz*dz;
            pass = (rsq < a.cutsq_max);
            if (store) keep = (rsq < a.inner_cutsq);
          }
          if (store) {
            const unsigned kb = __ballot_sync(FULLMASK, keep);
            if (keep) irow[icnt + __popc(kb & lanelt)] = static_cast<unsigned short>(ev[v]);
            icnt += __popc(kb);
          }
          const unsigned b = __ballot_sync(FULLMASK, pass);
          if (pass) s_q[qn + __popc(b & lanelt)] = ev[v];
          qn += __popc(b);
        }
      }
      while (qn >= 32 || (!more && qn > 0)) {
        __syncwarp();
        const int take = (qn >= 32) ? 32 : qn;
        if (lane < take) {
          const int e = s_q[qn - take + lane];
          const int loc = e & LJCL2_LMAX;
          const int sb = e >> LJCL2_LSB_SHIFT;
          const float4 p = s_xq[loc];
          const int jt = s_t[loc];
          const float dx = xi.x - p.x, dy = xi.y - p.y, dz = xi.z - p.z;
          const float rsq = dx*dx + dy*dy + dz*dz;
          const float rinv = rsqrtf(rsq);
          const float r2inv = rinv*rinv;
          float fpair = 0.0f;
          if (rsq < s_row[32+jt]) {
            const float r6inv = r2inv*r2inv*r2inv;
            fpair += a.special_lj[sb] * r6inv*(s_row[jt]*r6inv - s_row[16+jt]) * r2inv;
          }
          if (rsq < s_row[48+jt]) {
            const float qiqj = xi.w*p.w;
            if (RAT) {
              fpair += a.qqrd2e*qiqj*(a.special_coul[sb]*rinv*r2inv + g3*ljcl2_pmecorrF(g2*rsq));
            } else {
              const float r = rsq*rinv;
              const float grij = a.g_ewald*r;
              const float expm2 = expf(-grij*grij);
              const float t = 1.0f/(1.0f + static_cast<float>(EWALD_P)*grij);
              const float erfc = t*(static_cast<float>(A1)+t*(static_cast<float>(A2)+
                                 t*(static_cast<float>(A3)+t*(static_cast<float>(A4)+
                                 t*static_cast<float>(A5)))))*expm2;
              const float prefactor = a.qqrd2e*qiqj*rinv;
              float forcecoul = prefactor*(erfc + static_cast<float>(EWALD_F)*grij*expm2);
              const float factor_coul = a.special_coul[sb];
              if (factor_coul < 1.0f) forcecoul -= (1.0f - factor_coul)*prefactor;
              fpair += forcecoul*r2inv;
            }
          }
          const float px = dx*fpair, py = dy*fpair, pz = dz*fpair;
          fx += px; fy += py; fz += pz;
          if (FIX) {
            const long long vx = __float2ll_rn(-px*16777216.0f);
            const long long vy = __float2ll_rn(-py*16777216.0f);
            const long long vz = __float2ll_rn(-pz*16777216.0f);
            atomicAdd(&s_fi[loc],        static_cast<int>(vx & 0xffff));
            atomicAdd(&s_fi[MS + loc],   static_cast<int>(vx >> 16));
            atomicAdd(&s_fi[2*MS + loc], static_cast<int>(vy & 0xffff));
            atomicAdd(&s_fi[3*MS + loc], static_cast<int>(vy >> 16));
            atomicAdd(&s_fi[4*MS + loc], static_cast<int>(vz & 0xffff));
            atomicAdd(&s_fi[5*MS + loc], static_cast<int>(vz >> 16));
          } else {
            atomicAdd(&s_f[loc], -px);
            atomicAdd(&s_f[MS + loc], -py);
            atomicAdd(&s_f[2*MS + loc], -pz);
          }
        }
        qn -= take;
        __syncwarp();
      }
      if (!more) break;
    }

    for (int o = 16; o > 0; o >>= 1) {
      fx += __shfl_xor_sync(FULLMASK, fx, o);
      fy += __shfl_xor_sync(FULLMASK, fy, o);
      fz += __shfl_xor_sync(FULLMASK, fz, o);
    }
    if (lane == 0) {
      if (store) a.linnum[i] = icnt;
      float *fi = a.f + i*a.fs0;
      atomicAdd(fi, fx); atomicAdd(fi + a.fs1, fy); atomicAdd(fi + 2*a.fs1, fz);
    }
  }
  __syncthreads();

  // flush: one global atomic triple per distinct j that received a force
  for (int k = threadIdx.x; k < M; k += blockDim.x) {
    float vx, vy, vz;
    if (FIX) {
      const double sc = 1.0/16777216.0;
      vx = static_cast<float>(sc*static_cast<double>((static_cast<long long>(s_fi[MS+k]) << 16) + s_fi[k]));
      vy = static_cast<float>(sc*static_cast<double>((static_cast<long long>(s_fi[3*MS+k]) << 16) + s_fi[2*MS+k]));
      vz = static_cast<float>(sc*static_cast<double>((static_cast<long long>(s_fi[5*MS+k]) << 16) + s_fi[4*MS+k]));
    } else {
      vx = s_f[k]; vy = s_f[MS+k]; vz = s_f[2*MS+k];
    }
    if (vx != 0.0f || vy != 0.0f || vz != 0.0f) {
      float *fj = a.f + js[k]*a.fs0;
      atomicAdd(fj, vx); atomicAdd(fj + a.fs1, vy); atomicAdd(fj + 2*a.fs1, vz);
    }
  }
}
#endif

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

/* ======================================================================
   Cluster-pair list path (LMP_LJCL2_CLUSTER=1)

   GROMACS-style (nbnxm) layout on one rank, without ghost atoms:
   - list build (each LAMMPS reneighbor): local atoms are binned into x-y
     columns of width ~(8/density)^(1/3), sorted by z within a column, and
     each column is padded to a multiple of 8 slots.  i-cluster ci = slots
     [8ci, 8ci+8), j-cluster cj = slots [4cj, 4cj+4).  One warp per i-cluster
     searches the neighboring columns with bounding boxes at rlist = rc + skin
     and stores, per (ci, cj, periodic shift), the exact 32-bit pair mask.
     Each unordered pair is kept once: shift > 13 (the zero shift), or zero
     shift and a higher owning i-cluster, or the same i-cluster and a higher
     slot.  Tiles containing special pairs carry 2-bit codes in a side array.
   - pruning: an inner list keeps the slot pairs within rc + delta.  It is
     redone on the device whenever an atom moved more than delta/2 since the
     last prune (no host synchronization; the prune kernel checks a flag).
   - force: one warp per i-cluster, lane = kj*8 + ki.  The i-force stays in
     registers, the j-force is reduced across the 8 lanes sharing kj with 3
     shuffle steps, and one atomic per j atom goes to a cluster-ordered force
     array that is scattered to f once per step.
   ====================================================================== */

static constexpr int CL_CI = 8, CL_CJ = 4;
static constexpr int CL_ZBITS = 20;

template<class DeviceType>
struct ClKey {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typename AT::t_kkfloat_1d_3_lr_randomread x;
  Kokkos::View<unsigned int*, DeviceType> key;
  Kokkos::View<int*, DeviceType> perm, colcnt;
  double blo[3], icw[2], izq;
  int ncx, ncy;
  KOKKOS_INLINE_FUNCTION void operator()(const int i) const {
    int cx = static_cast<int>((x(i,0) - blo[0]) * icw[0]);
    int cy = static_cast<int>((x(i,1) - blo[1]) * icw[1]);
    cx = cx < 0 ? 0 : (cx >= ncx ? ncx-1 : cx);
    cy = cy < 0 ? 0 : (cy >= ncy ? ncy-1 : cy);
    double zf = (x(i,2) - blo[2]) * izq;
    const double zmax = static_cast<double>((1u << CL_ZBITS) - 1u);
    zf = zf < 0.0 ? 0.0 : (zf > zmax ? zmax : zf);
    const unsigned int col = static_cast<unsigned int>(cx + ncx*cy);
    key(i) = (col << CL_ZBITS) | static_cast<unsigned int>(zf);
    perm(i) = i;
    Kokkos::atomic_add(&colcnt(col), 1);
  }
};

template<class DeviceType>
struct ClSlot {
  typedef DeviceType device_type;
  Kokkos::View<unsigned int*, DeviceType> key;
  Kokkos::View<int*, DeviceType> perm, colfirst, colslot, atom;
  KOKKOS_INLINE_FUNCTION void operator()(const int p) const {
    const int col = static_cast<int>(key(p) >> CL_ZBITS);
    atom(colslot(col) + p - colfirst(col)) = perm(p);
  }
};

template<class DeviceType>
struct ClIcol {
  typedef DeviceType device_type;
  Kokkos::View<int*, DeviceType> colslot, icol;
  KOKKOS_INLINE_FUNCTION void operator()(const int c) const {
    for (int ic = colslot(c)/CL_CI; ic < colslot(c+1)/CL_CI; ic++) icol(ic) = c;
  }
};

template<class DeviceType>
struct ClGather {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typename AT::t_kkfloat_1d_3_lr_randomread x;
  typename AT::t_kkfloat_1d_randomread q;
  typename AT::t_int_1d_randomread type;
  Kokkos::View<int*, DeviceType> atom;
  Kokkos::View<LJCL2XQ*, DeviceType> xq;
  Kokkos::View<unsigned char*, DeviceType> t8;
  KOKKOS_INLINE_FUNCTION void operator()(const int s) const {
    const int a = atom(s);
    LJCL2XQ r;
    if (a >= 0) {
      r.x = x(a,0); r.y = x(a,1); r.z = x(a,2); r.q = q(a);
      t8(s) = static_cast<unsigned char>(type(a));
    } else {
      r.x = r.y = r.z = r.q = 0;
      t8(s) = 1;
    }
    xq(s) = r;
  }
};

template<class DeviceType>
struct ClBBox {
  typedef DeviceType device_type;
  Kokkos::View<int*, DeviceType> atom;
  Kokkos::View<LJCL2XQ*, DeviceType> xq;
  Kokkos::View<LJCL2BB*, DeviceType> jbb;
  KOKKOS_INLINE_FUNCTION void operator()(const int jc) const {
    LJCL2BB b;
    for (int d = 0; d < 3; d++) { b.lo[d] = 1.0e30f; b.hi[d] = -1.0e30f; }
    for (int k = 0; k < CL_CJ; k++) {
      const int s = jc*CL_CJ + k;
      if (atom(s) < 0) continue;
      const LJCL2XQ r = xq(s);
      const float c[3] = {static_cast<float>(r.x), static_cast<float>(r.y),
                          static_cast<float>(r.z)};
      for (int d = 0; d < 3; d++) {
        b.lo[d] = c[d] < b.lo[d] ? c[d] : b.lo[d];
        b.hi[d] = c[d] > b.hi[d] ? c[d] : b.hi[d];
      }
    }
    jbb(jc) = b;
  }
};

// bounding-box distance^2 between i box a and j box b shifted by (sx,sy,sz)
KOKKOS_INLINE_FUNCTION float cl_bbdist2(const LJCL2BB &a, const LJCL2BB &b,
                                        const float sx, const float sy, const float sz) {
  const float sh[3] = {sx, sy, sz};
  float d2 = 0.0f;
  for (int d = 0; d < 3; d++) {
    const float lo = b.lo[d] + sh[d], hi = b.hi[d] + sh[d];
    float g = lo - a.hi[d];
    const float g2 = a.lo[d] - hi;
    g = g > g2 ? g : g2;
    if (g > 0.0f) d2 += g*g;
  }
  return d2;
}

template<class DeviceType>
struct ClSearch {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef Kokkos::TeamPolicy<DeviceType, Kokkos::IndexType<int>> policy_type;
  typedef typename policy_type::member_type member_type;

  Kokkos::View<int*, DeviceType> atom, colslot, icol, nout, ctl;
  Kokkos::View<LJCL2XQ*, DeviceType> xq;
  Kokkos::View<LJCL2BB*, DeviceType> jbb;
  Kokkos::View<LJCL2CTile**, Kokkos::LayoutRight, DeviceType> out;
  Kokkos::View<unsigned int*[2], Kokkos::LayoutRight, DeviceType> sc;
  typename AT::t_tagint_1d_randomread tag;
  typename AT::t_tagint_2d_randomread special;
  typename AT::t_int_2d_randomread nspecial;
  double blo[3], prd[3], cw[2];
  float rlist2;
  int ncx, ncy, nx, ny, maxout, maxsc;

  KOKKOS_INLINE_FUNCTION
  static int shmem_bytes() { return 8*(3*sizeof(float) + 4*sizeof(int) + 2*sizeof(tagint)) + 64; }

  KOKKOS_INLINE_FUNCTION
  void operator()(const member_type& team) const {
#if defined(__CUDA_ARCH__)
    char *sh = static_cast<char *>(team.team_shmem().get_shmem(shmem_bytes()));
    float *sx_i = reinterpret_cast<float *>(sh);          // [3][8]
    int *s_ai = reinterpret_cast<int *>(sx_i + 24);       // atom index, -1 = pad
    tagint *s_lo = reinterpret_cast<tagint *>(sh + ((24*sizeof(float) + 8*sizeof(int) + 15) & ~15));
    tagint *s_hi = s_lo + 8;
    const int lane = threadIdx.x;
    const unsigned FULLMASK = 0xffffffffu;
    const unsigned lanelt = (1u << lane) - 1u;
    const int ci = team.league_rank();

    if (lane < CL_CI) {
      const int s = ci*CL_CI + lane;
      const int a = atom(s);
      const LJCL2XQ r = xq(s);
      sx_i[lane] = r.x; sx_i[8+lane] = r.y; sx_i[16+lane] = r.z;
      s_ai[lane] = a;
      tagint lo = MAXTAGINT, hi = 0;
      if (a >= 0) {
        const int n3 = nspecial(a,2);
        for (int k = 0; k < n3; k++) {
          const tagint t = special(a,k);
          lo = t < lo ? t : lo; hi = t > hi ? t : hi;
        }
      }
      s_lo[lane] = lo; s_hi[lane] = hi;
    }
    __syncwarp();

    LJCL2BB ib;
    {
      const LJCL2BB b0 = jbb(2*ci), b1 = jbb(2*ci+1);
      for (int d = 0; d < 3; d++) {
        ib.lo[d] = b0.lo[d] < b1.lo[d] ? b0.lo[d] : b1.lo[d];
        ib.hi[d] = b0.hi[d] > b1.hi[d] ? b0.hi[d] : b1.hi[d];
      }
    }
    const int c0 = icol(ci);
    const int cx0 = c0 % ncx, cy0 = c0 / ncx;
    int cnt = 0;                                          // warp-uniform

    for (int dy = -ny; dy <= ny; dy++) {
      int cy = cy0 + dy, sy = 0;
      if (cy < 0) { cy += ncy; sy = -1; } else if (cy >= ncy) { cy -= ncy; sy = 1; }
      const float shy = static_cast<float>(sy*prd[1]);
      const float ylo = static_cast<float>(blo[1] + cy*cw[1]) + shy, yhi = ylo + static_cast<float>(cw[1]);
      float gy = ylo - ib.hi[1]; { const float g2 = ib.lo[1] - yhi; gy = gy > g2 ? gy : g2; }
      gy = gy > 0.0f ? gy : 0.0f;
      if (gy*gy >= rlist2) continue;
      for (int dx = -nx; dx <= nx; dx++) {
        int cx = cx0 + dx, sx = 0;
        if (cx < 0) { cx += ncx; sx = -1; } else if (cx >= ncx) { cx -= ncx; sx = 1; }
        const float shx = static_cast<float>(sx*prd[0]);
        const float xlo = static_cast<float>(blo[0] + cx*cw[0]) + shx, xhi = xlo + static_cast<float>(cw[0]);
        float gx = xlo - ib.hi[0]; { const float g2 = ib.lo[0] - xhi; gx = gx > g2 ? gx : g2; }
        gx = gx > 0.0f ? gx : 0.0f;
        if (gx*gx + gy*gy >= rlist2) continue;
        const int col = cx + ncx*cy;
        const int jc0 = colslot(col)/CL_CJ, jc1 = colslot(col+1)/CL_CJ;
        for (int jb = jc0; jb < jc1; jb += 32) {
          const int jc = jb + lane;
          for (int sz = -1; sz <= 1; sz++) {
            const float shz = static_cast<float>(sz*prd[2]);
            const int shift = (sx+1) + 3*(sy+1) + 9*(sz+1);
            unsigned int mask = 0;
            if (jc < jc1) {
              const int own = jc >> 1;
              const bool keep = (shift > 13) || ((shift == 13) && (own >= ci));
              if (keep && (cl_bbdist2(ib, jbb(jc), shx, shy, shz) < rlist2)) {
                for (int kj = 0; kj < CL_CJ; kj++) {
                  const int js = jc*CL_CJ + kj;
                  if (atom(js) < 0) continue;
                  const LJCL2XQ r = xq(js);
                  const float xj = static_cast<float>(r.x) + shx;
                  const float yj = static_cast<float>(r.y) + shy;
                  const float zj = static_cast<float>(r.z) + shz;
                  for (int ki = 0; ki < CL_CI; ki++) {
                    if (s_ai[ki] < 0) continue;
                    if ((shift == 13) && (own == ci) && (ci*CL_CI + ki >= js)) continue;
                    const float ddx = sx_i[ki] - xj, ddy = sx_i[8+ki] - yj, ddz = sx_i[16+ki] - zj;
                    if (ddx*ddx + ddy*ddy + ddz*ddz < rlist2) mask |= 1u << (kj*8 + ki);
                  }
                }
              }
            }
            // special codes for the slot pairs of this tile (window prefilter)
            int sidx = -1;
            if (mask) {
              unsigned int w0 = 0, w1 = 0;
              for (int kj = 0; kj < CL_CJ; kj++) {
                if (!((mask >> (kj*8)) & 0xffu)) continue;
                const int ja = atom(jc*CL_CJ + kj);
                const tagint tj = tag(ja);
                for (int ki = 0; ki < CL_CI; ki++) {
                  const int b = kj*8 + ki;
                  if (!((mask >> b) & 1u)) continue;
                  if ((tj < s_lo[ki]) || (tj > s_hi[ki])) continue;
                  const int ia = s_ai[ki];
                  const int n1 = nspecial(ia,0), n2 = nspecial(ia,1), n3 = nspecial(ia,2);
                  int code = 0;
                  for (int k = 0; k < n3; k++)
                    if (special(ia,k) == tj) { code = (k < n1) ? 1 : ((k < n2) ? 2 : 3); break; }
                  if (code) {
                    if (b < 16) w0 |= static_cast<unsigned int>(code) << (2*b);
                    else w1 |= static_cast<unsigned int>(code) << (2*(b-16));
                  }
                }
              }
              if (w0 | w1) {
                sidx = Kokkos::atomic_fetch_add(&ctl(2), 1);
                if (sidx < maxsc) { sc(sidx,0) = w0; sc(sidx,1) = w1; }
                else { ctl(3) = 1; sidx = -1; }
              }
            }
            const unsigned bal = __ballot_sync(FULLMASK, mask != 0u);
            if (mask) {
              const int idx = cnt + __popc(bal & lanelt);
              if (idx < maxout) {
                LJCL2CTile t; t.cj = jc; t.mask = mask; t.sidx = sidx; t.shift = shift;
                out(ci, idx) = t;
              }
            }
            cnt += __popc(bal);
          }
        }
      }
    }
    if (lane == 0) {
      nout(ci) = cnt < maxout ? cnt : maxout;
      if (cnt > maxout) Kokkos::atomic_max(&ctl(0), cnt);
    }
#else
    (void) team;
    Kokkos::abort("cluster search requires CUDA");
#endif
  }
};

// per step: has any atom moved more than delta/2 since the last prune?
template<class DeviceType>
struct ClDisp {
  typedef DeviceType device_type;
  Kokkos::View<int*, DeviceType> atom, ctl;
  Kokkos::View<LJCL2XQ*, DeviceType> xq;
  Kokkos::View<float*[4], Kokkos::LayoutRight, DeviceType> xp;
  float half2;
  KOKKOS_INLINE_FUNCTION void operator()(const int s) const {
    if (atom(s) < 0) return;
    const LJCL2XQ r = xq(s);
    const float dx = r.x - xp(s,0), dy = r.y - xp(s,1), dz = r.z - xp(s,2);
    if (dx*dx + dy*dy + dz*dz > half2) ctl(4) = 1;      // benign race: all write 1
  }
};

template<class DeviceType>
struct ClFlagReset {
  typedef DeviceType device_type;
  Kokkos::View<int*, DeviceType> ctl;
  KOKKOS_INLINE_FUNCTION void operator()(const int) const { ctl(4) = 0; ctl(6) = 0; }
};

KOKKOS_INLINE_FUNCTION void cl_shift(const int shift, const float prd0, const float prd1,
                                     const float prd2, float &sx, float &sy, float &sz) {
  sx = static_cast<float>(shift % 3 - 1) * prd0;
  sy = static_cast<float>((shift / 3) % 3 - 1) * prd1;
  sz = static_cast<float>(shift / 9 - 1) * prd2;
}

// staged tile header: each lane of a warp decodes one tile of a 32-tile batch
// (shift vector, special-code words) so the per-tile loop only issues two
// broadcast shared-memory loads
struct alignas(16) ClHdr {
  int cj; unsigned int mask; int sidx; int shift;
  float sx, sy, sz; unsigned int sc;        // sc: packed for the force kernel, see below
};

template<class DeviceType>
struct ClPrune {
  typedef DeviceType device_type;
  typedef Kokkos::TeamPolicy<DeviceType, Kokkos::IndexType<int>> policy_type;
  typedef typename policy_type::member_type member_type;
  Kokkos::View<int*, DeviceType> nout, nin, ctl;
  Kokkos::View<LJCL2XQ*, DeviceType> xq;
  Kokkos::View<float*[4], Kokkos::LayoutRight, DeviceType> xp;
  Kokkos::View<LJCL2CTile**, Kokkos::LayoutRight, DeviceType> out, in;
  float rin2, prd0, prd1, prd2;
  int nicl, warps_per_team;

  KOKKOS_INLINE_FUNCTION
  static int shmem_bytes(int warps) { return warps*32*(int)sizeof(ClHdr) + 64; }

  KOKKOS_INLINE_FUNCTION
  void operator()(const member_type& team) const {
#if defined(__CUDA_ARCH__)
    if (!ctl(4)) return;                                 // launch-uniform
    const int w = team.team_rank();
    const int ci = team.league_rank()*warps_per_team + w;
    if (ci >= nicl) return;                              // warp-uniform
    ClHdr *hdr = static_cast<ClHdr *>(
        team.team_shmem().get_shmem(shmem_bytes(warps_per_team))) + w*32;
    const int lane = threadIdx.x;
    const unsigned FULLMASK = 0xffffffffu;
    const unsigned lanelt = (1u << lane) - 1u;
    const int ki = lane & 7, kj = lane >> 3;
    const LJCL2XQ i4 = xq(ci*CL_CI + ki);
    const int n = nout(ci);
    int cnt = 0;                                         // warp-uniform
    for (int t0 = 0; t0 < n; t0 += 32) {
      const int nt = (n - t0 < 32) ? (n - t0) : 32;
      LJCL2CTile mine = {0, 0u, -1, 13};
      __syncwarp();
      if (lane < nt) {
        mine = out(ci, t0 + lane);
        ClHdr hh;
        hh.cj = mine.cj; hh.mask = mine.mask; hh.sidx = mine.sidx; hh.shift = mine.shift;
        cl_shift(mine.shift, prd0, prd1, prd2, hh.sx, hh.sy, hh.sz);
        hh.sc = 0u;
        hdr[lane] = hh;
      }
      __syncwarp();
      unsigned int newmask = 0u;                         // lane s keeps tile s's result
      for (int s = 0; s < nt; s++) {
        const ClHdr h = hdr[s];
        const LJCL2XQ j4 = xq(h.cj*CL_CJ + kj);
        const float dx = i4.x - (j4.x + h.sx), dy = i4.y - (j4.y + h.sy), dz = i4.z - (j4.z + h.sz);
        const bool keep = ((h.mask >> lane) & 1u) && (dx*dx + dy*dy + dz*dz < rin2);
        const unsigned b = __ballot_sync(FULLMASK, keep);
        if (lane == s) newmask = b;
      }
      // one coalesced write of the kept tiles of this batch
      const unsigned kept = __ballot_sync(FULLMASK, newmask != 0u);
      if (newmask) {
        LJCL2CTile t = mine; t.mask = newmask;
        in(ci, cnt + __popc(kept & lanelt)) = t;
      }
      cnt += __popc(kept);
    }
    if (lane == 0) nin(ci) = cnt;
    if (lane < CL_CI) {
      const LJCL2XQ r = xq(ci*CL_CI + lane);
      xp(ci*CL_CI + lane, 0) = r.x; xp(ci*CL_CI + lane, 1) = r.y; xp(ci*CL_CI + lane, 2) = r.z;
    }
#else
    (void) team;
    Kokkos::abort("cluster prune requires CUDA");
#endif
  }
};

template<class DeviceType, bool RAT, bool EV, bool UNICUT>
struct ClForce : public LJCL2Common<DeviceType> {
  typedef DeviceType device_type;
  typedef Kokkos::TeamPolicy<DeviceType, Kokkos::IndexType<int>> policy_type;
  typedef typename policy_type::member_type member_type;
  Kokkos::View<int*, DeviceType> atom, nin;
  Kokkos::View<LJCL2XQ*, DeviceType> xq;
  Kokkos::View<unsigned char*, DeviceType> t8;
  Kokkos::View<LJCL2CTile**, Kokkos::LayoutRight, DeviceType> in;
  Kokkos::View<unsigned int*[2], Kokkos::LayoutRight, DeviceType> sc;
  Kokkos::View<KK_ACC_FLOAT*[3], Kokkos::LayoutRight, DeviceType> fcl;
  Kokkos::View<double*, DeviceType> ev;
  Kokkos::View<int*, DeviceType> ctl;                   // ctl(6): persistent work counter
  float prd0, prd1, prd2;
  KK_FLOAT cutsq_uni;                                   // UNICUT: the one cutoff^2
  int nicl, warps_per_team, persist = 0;

  // per warp: 32 staged headers + the two special-code words of each tile
  KOKKOS_INLINE_FUNCTION
  static int shmem_bytes(int warps) { return warps*32*((int)sizeof(ClHdr) + 8) + 64; }

  KOKKOS_INLINE_FUNCTION
  void operator()(const member_type& team) const {
#if defined(__CUDA_ARCH__)
    const int w = team.team_rank();
    char *sh = static_cast<char *>(team.team_shmem().get_shmem(shmem_bytes(warps_per_team)));
    ClHdr *hdr = reinterpret_cast<ClHdr *>(sh) + w*32;
    unsigned int *scw = reinterpret_cast<unsigned int *>(
        reinterpret_cast<ClHdr *>(sh) + warps_per_team*32) + w*64;
    const int lane = threadIdx.x;
    const unsigned FULLMASK = 0xffffffffu;
    const int ki = lane & 7, kj = lane >> 3;

    // persistent mode: a grid that fits on the GPU at once, warps fetch
    // i-clusters from a counter, so no blocks queue ahead of kernels on other
    // streams; otherwise one i-cluster per warp
    int ci = team.league_rank()*warps_per_team + w;
    if (persist) {
      int c = 0;
      if (lane == 0) c = Kokkos::atomic_fetch_add(&ctl(6), 1);
      ci = __shfl_sync(FULLMASK, c, 0);
    }
    while (ci < nicl) {                                  // warp-uniform
      const int islot = ci*CL_CI + ki;
      const LJCL2XQ i4 = xq(islot);
      const int itype = static_cast<int>(t8(islot));
      const int half = lane >> 4, sbshift = 2*(lane & 15);

      KK_ACC_FLOAT fix = 0, fiy = 0, fiz = 0;
      double e_vdwl = 0.0, e_coul = 0.0, v[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

      const int n = nin(ci);
      for (int t0 = 0; t0 < n; t0 += 32) {
        const int nt = (n - t0 < 32) ? (n - t0) : 32;
        __syncwarp();
        if (lane < nt) {
          const LJCL2CTile t = in(ci, t0 + lane);
          ClHdr hh;
          hh.cj = t.cj; hh.mask = t.mask; hh.sidx = t.sidx; hh.shift = t.shift;
          cl_shift(t.shift, prd0, prd1, prd2, hh.sx, hh.sy, hh.sz);
          hh.sc = 0u;
          hdr[lane] = hh;
          unsigned int w0 = 0u, w1 = 0u;
          if (t.sidx >= 0) { w0 = sc(t.sidx,0); w1 = sc(t.sidx,1); }
          scw[2*lane] = w0; scw[2*lane+1] = w1;
        }
        __syncwarp();
        for (int s = 0; s < nt; s++) {
          const ClHdr h = hdr[s];
          const int jslot = h.cj*CL_CJ + kj;
          const LJCL2XQ j4 = xq(jslot);
          const KK_FLOAT delx = i4.x - (j4.x + h.sx);
          const KK_FLOAT dely = i4.y - (j4.y + h.sy);
          const KK_FLOAT delz = i4.z - (j4.z + h.sz);
          const KK_FLOAT rsq = delx*delx + dely*dely + delz*delz;
          KK_ACC_FLOAT fx = 0, fy = 0, fz = 0;
          bool m = ((h.mask >> lane) & 1u);
          if (UNICUT) m = m && (rsq < cutsq_uni);
          if (m) {
            const int jtype = static_cast<int>(t8(jslot));
            if (UNICUT || (rsq < this->m_cutsq[itype][jtype])) {
              const int sb = (scw[2*s + half] >> sbshift) & 3;
              const KK_FLOAT qiqj = i4.q*j4.q;
              const KK_FLOAT fpair = this->template pair_fpair<false,false,RAT>(rsq, itype, jtype, qiqj, sb);
              fx = static_cast<KK_ACC_FLOAT>(delx*fpair);
              fy = static_cast<KK_ACC_FLOAT>(dely*fpair);
              fz = static_cast<KK_ACC_FLOAT>(delz*fpair);
              if (EV) {
                const KK_FLOAT r2inv = static_cast<KK_FLOAT>(1.0)/rsq;
                if (rsq < this->m_cut_ljsq[itype][jtype]) {
                  const KK_FLOAT r6inv = r2inv*r2inv*r2inv;
                  e_vdwl += this->special_lj[sb] *
                    (r6inv*(this->m_params[itype][jtype].lj3*r6inv - this->m_params[itype][jtype].lj4)
                     - this->m_params[itype][jtype].offset);
                }
                if (rsq < this->m_cut_coulsq[itype][jtype]) {
                  const double r = sqrt(static_cast<double>(rsq));
                  const double grij = this->g_ewald * r;
                  const double expm2 = exp(-grij*grij);
                  const double t = 1.0 / (1.0 + EWALD_P*grij);
                  const double erfc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2;
                  const double prefactor = this->qqrd2e * qiqj / r;
                  double ec = prefactor*erfc;
                  if (sb) ec -= (1.0 - this->special_coul[sb])*prefactor;
                  e_coul += ec;
                }
                v[0] += delx*delx*fpair; v[1] += dely*dely*fpair; v[2] += delz*delz*fpair;
                v[3] += delx*dely*fpair; v[4] += delx*delz*fpair; v[5] += dely*delz*fpair;
              }
            }
          }
          fix += fx; fiy += fy; fiz += fz;
          // j-force: reduce over the 8 lanes that share kj
          for (int o = 1; o < 8; o <<= 1) {
            fx += __shfl_xor_sync(FULLMASK, fx, o);
            fy += __shfl_xor_sync(FULLMASK, fy, o);
            fz += __shfl_xor_sync(FULLMASK, fz, o);
          }
          if ((ki == 0) && ((h.mask >> (kj*8)) & 0xffu)) {
            Kokkos::atomic_add(&fcl(jslot,0), -fx);
            Kokkos::atomic_add(&fcl(jslot,1), -fy);
            Kokkos::atomic_add(&fcl(jslot,2), -fz);
          }
        }
      }
      // i-force: reduce over the 4 lanes that share ki
      for (int o = 8; o < 32; o <<= 1) {
        fix += __shfl_xor_sync(FULLMASK, fix, o);
        fiy += __shfl_xor_sync(FULLMASK, fiy, o);
        fiz += __shfl_xor_sync(FULLMASK, fiz, o);
      }
      if ((lane < CL_CI) && (atom(islot) >= 0)) {
        Kokkos::atomic_add(&fcl(islot,0), fix);
        Kokkos::atomic_add(&fcl(islot,1), fiy);
        Kokkos::atomic_add(&fcl(islot,2), fiz);
      }
      if (EV) {
        double acc[8] = {e_vdwl, e_coul, v[0], v[1], v[2], v[3], v[4], v[5]};
        for (int k = 0; k < 8; k++) {
          for (int o = 16; o > 0; o >>= 1) acc[k] += __shfl_xor_sync(FULLMASK, acc[k], o);
          if (lane == 0) Kokkos::atomic_add(&ev(k), acc[k]);
        }
      }
      if (!persist) break;
      int c = 0;
      if (lane == 0) c = Kokkos::atomic_fetch_add(&ctl(6), 1);
      ci = __shfl_sync(FULLMASK, c, 0);
    }
#else
    (void) team;
    Kokkos::abort("cluster force requires CUDA");
#endif
  }
};

template<class DeviceType>
struct ClScatter {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typename AT::t_kkacc_1d_3 f;
  Kokkos::View<int*, DeviceType> atom;
  Kokkos::View<KK_ACC_FLOAT*[3], Kokkos::LayoutRight, DeviceType> fcl;
  KOKKOS_INLINE_FUNCTION void operator()(const int s) const {
    const int a = atom(s);
    if (a >= 0) { f(a,0) += fcl(s,0); f(a,1) += fcl(s,1); f(a,2) += fcl(s,2); }
    fcl(s,0) = 0; fcl(s,1) = 0; fcl(s,2) = 0;
  }
};

}  // anonymous namespace

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairLJCutCoulLong2Kokkos<DeviceType>::PairLJCutCoulLong2Kokkos(LAMMPS *lmp) :
  PairLJCutCoulLongKokkos<DeviceType>(lmp) {}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairLJCutCoulLong2Kokkos<DeviceType>::~PairLJCutCoulLong2Kokkos()
{
  if (this->copymode) return;
#if defined(KOKKOS_ENABLE_CUDA)
  if (cl_stream) {
    cl_space = typename DeviceType::execution_space();   // drop the instance first
    cudaStreamSynchronize(cl_stream);
    cudaEventDestroy(cl_ev_in);
    cudaEventDestroy(cl_ev_out);
    cudaStreamDestroy(cl_stream);
  }
#endif
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairLJCutCoulLong2Kokkos<DeviceType>::init_style()
{
  PairLJCutCoulLongKokkos<DeviceType>::init_style();

  if (cl_on < 0) {
    const char *e = std::getenv("LMP_LJCL2_CLUSTER");
    cl_on = (e && std::atoi(e)) ? 1 : 0;
    if (const char *d = std::getenv("LMP_CL_DELTA")) cl_delta = std::atof(d);
    if (const char *d = std::getenv("LMP_CL_COLSCALE")) cl_colscale = std::atof(d);
    cl_stats = (std::getenv("LMP_CL_STATS") != nullptr);
    if (const char *a = std::getenv("LMP_CL_ASYNC")) cl_async = std::atoi(a);
    if (const char *p = std::getenv("LMP_CL_PERSIST")) cl_persist = std::atoi(p);
    if (const char *k = std::getenv("LMP_CL_KFIRST")) cl_kfirst = std::atoi(k);
  }
  if (!cl_on) return;

  bool cuda = false;
#if defined(KOKKOS_ENABLE_CUDA)
  cuda = std::is_same<typename DeviceType::execution_space, Kokkos::Cuda>::value;
#endif
  if (!cuda)
    this->error->all(FLERR, "LMP_LJCL2_CLUSTER requires lj/cut/coul/long2/kk on a CUDA GPU");
  if (this->comm->nprocs != 1)
    this->error->all(FLERR, "LMP_LJCL2_CLUSTER currently supports a single MPI rank only");
  if (this->domain->triclinic ||
      !(this->domain->xperiodic && this->domain->yperiodic && this->domain->zperiodic))
    this->error->all(FLERR, "LMP_LJCL2_CLUSTER requires an orthogonal, fully periodic box");
  if (!this->force->newton_pair || (this->neighflag != HALFTHREAD))
    this->error->all(FLERR, "LMP_LJCL2_CLUSTER requires newton on and 'neigh half'");
  if (this->lmp->kokkos->neigh_cluster)
    this->error->all(FLERR, "LMP_LJCL2_CLUSTER and 'neigh/cluster yes' are exclusive");
  if (cl_delta <= 0.0)
    this->error->all(FLERR, "LMP_CL_DELTA must be > 0");

  // the cluster list replaces the flat pair list: build that only on demand
  // (never, in practice), and compute the virial explicitly since pair forces
  // no longer go to ghost atoms
  this->neighbor->find_request(this)->set_occasional(1);
  this->no_virial_fdotr_compute = 1;

  // a second, non-blocking stream for the force-only steps.  Non-blocking, so
  // that synchronous copies on the legacy default stream elsewhere do not
  // serialize against it.
#if defined(KOKKOS_ENABLE_CUDA)
  if (cl_persist > 0 && cl_nsm == 0) {
    int dev = 0;
    cudaGetDevice(&dev);
    cudaDeviceGetAttribute(&cl_nsm, cudaDevAttrMultiProcessorCount, dev);
  }
  if (cl_async && !cl_stream) {
    if constexpr (std::is_same<typename DeviceType::execution_space, Kokkos::Cuda>::value) {
      cudaStreamCreateWithFlags(&cl_stream, cudaStreamNonBlocking);
      cudaEventCreateWithFlags(&cl_ev_in, cudaEventDisableTiming);
      cudaEventCreateWithFlags(&cl_ev_out, cudaEventDisableTiming);
      cl_space = Kokkos::Cuda(cl_stream);
    }
  }
#endif
}

/* ----------------------------------------------------------------------
   join the cluster force kernels left running on cl_stream, then add the
   cluster-ordered forces into f on the default stream
------------------------------------------------------------------------- */

template<class DeviceType>
void PairLJCutCoulLong2Kokkos<DeviceType>::finish_async_compute()
{
  if (!cl_pending) return;
  cl_pending = false;
#if defined(KOKKOS_ENABLE_CUDA)
  if constexpr (std::is_same<typename DeviceType::execution_space, Kokkos::Cuda>::value)
    cudaStreamWaitEvent(Kokkos::Cuda().cuda_stream(), cl_ev_out, 0);
#endif
  cl_scatter();
}

template<class DeviceType>
void PairLJCutCoulLong2Kokkos<DeviceType>::cl_scatter()
{
  ClScatter<DeviceType> fs;
  fs.f = this->atomKK->k_f.template view<DeviceType>();
  fs.atom = d_cl_atom; fs.fcl = d_cl_f;
  Kokkos::parallel_for("PairLJCutCoulLong2::cl_scatter", Kokkos::RangePolicy<DeviceType>(0, cl_nslot), fs);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairLJCutCoulLong2Kokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
  if (cl_on == 1) {
    if constexpr (std::is_same<typename DeviceType::execution_space,
                               LMPDeviceType::execution_space>::value) {
      this->eflag = eflag_in;
      this->vflag = vflag_in;
      this->ev_init(eflag_in, vflag_in, 0);
      if (this->eflag_atom || this->vflag_atom)
        this->error->all(FLERR, "LMP_LJCL2_CLUSTER does not support per-atom energy/virial");
      this->atomKK->sync(this->execution_space,
                         X_MASK | F_MASK | TYPE_MASK | Q_MASK | TAG_MASK | SPECIAL_MASK);
      cl_step(eflag_in, vflag_in);
      this->atomKK->modified(this->execution_space, F_MASK);
      return;
    }
  }

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

  // 'package kokkos neigh/cluster yes' morphs the request to a FULL-style list
  // (cluster-level Newton) while neighflag stays HALFTHREAD, so none of the
  // half-list kernels below may walk it: the base-class cluster path does.
  if (this->lmp->kokkos->neigh_cluster) {
    PairLJCutCoulLongKokkos<DeviceType>::compute(eflag_in, vflag_in);
    return;
  }

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

  // ---- experimental shared-memory tiled kernel (needs the packed arrays)
  if (local_on < 0) {
    local_on = 0;
    if (const char *e = std::getenv("LMP_LJCL2_LOCAL")) local_on = std::atoi(e);
    if (const char *e = std::getenv("LMP_LJCL2_LNI")) local_ni = std::atoi(e);
    if (const char *e = std::getenv("LMP_LJCL2_LVEC")) local_vec = std::atoi(e);
    if (local_vec != 2 && local_vec != 4) local_vec = 1;
    local_stats = (std::getenv("LMP_LOCAL_STATS") != nullptr);
    if (local_on < 0 || local_on > 2) local_on = 0;
    if (local_ni < 1) local_ni = 1;
    bool ok = false;
#if defined(LJCL2_HAVE_LOCAL)
    ok = std::is_same<typename DeviceType::execution_space, Kokkos::Cuda>::value &&
         (sizeof(KK_FLOAT) == sizeof(float)) && (sizeof(KK_ACC_FLOAT) == sizeof(float));
#endif
    if (local_on && (!ok || !use_half || union_ci)) {
      if (this->comm->me == 0)
        this->error->warning(FLERR, "LMP_LJCL2_LOCAL needs CUDA, single precision and a "
                                    "half list without UNION -- ignored");
      local_on = 0;
    }
  }

  // ---- experimental packed j-side gather
  if (pack_on < 0) {
    pack_on = 0;
    if (const char *e = std::getenv("LMP_LJCL2_PACK")) pack_on = std::atoi(e) ? 1 : 0;
  }
  if (local_on) pack_on = 1;
  if (pack_on) pack_refresh();

  if (local_on && use_half) {
    if (this->neighbor->lastcall > local_built_step) {
      local_build();
      local_built_step = this->neighbor->lastcall;
    }
    if (local_on) {     // the build may have disabled it
      local_compute();
      return;
    }
  }

  // ---- atomic-free warp-aggregated inner-list append (dual cutoff only)
  if (ballot_on < 0) {
    ballot_on = 0;
    if (const char *e = std::getenv("LMP_LJCL2_BALLOT")) ballot_on = std::atoi(e) ? 1 : 0;
  }

  // ---- kernel variant.  1-3 are diagnostics whose forces are WRONG by
  // construction; they are only ever selected explicitly, and the style warns
  // once so no benchmark log can be mistaken for a physical run.
  if (variant < 0) {
    variant = LJCL2_VAR_NORMAL;
    if (const char *e = std::getenv("LMP_LJCL2_VARIANT")) variant = std::atoi(e);
    if (variant < 0 || variant > LJCL2_VAR_RATIONAL) {
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

  // ---- warp-compacted kernel (CUDA only; exclusive with dual and pack)
  if (compact_on < 0) {
    compact_on = 0;
    compact_smem = 0;
    if (const char *e = std::getenv("LMP_LJCL2_COMPACT")) compact_on = std::atoi(e) ? 1 : 0;
    if (const char *e = std::getenv("LMP_LJCL2_SMEM")) compact_smem = std::atoi(e) ? 1 : 0;
    if (const char *e = std::getenv("LMP_LJCL2_APW")) compact_apw = std::atoi(e);
    if (compact_apw < 1) compact_apw = 1;
    bool cuda = false;
#if defined(KOKKOS_ENABLE_CUDA)
    cuda = std::is_same<typename DeviceType::execution_space, Kokkos::Cuda>::value;
#endif
    if (compact_on && !cuda) {
      if (this->comm->me == 0)
        this->error->warning(FLERR, "LMP_LJCL2_COMPACT needs the CUDA backend -- ignored");
      compact_on = 0;
    }
    if (compact_on && compact_smem && (dual || pack_on)) {
      if (this->comm->me == 0)
        this->error->warning(FLERR, "LMP_LJCL2_SMEM is not combined with LMP_LJCL2_DUAL/PACK "
                                    "-- using constant-memory tables");
      compact_smem = 0;
    }
  }

  if (compact_on) {
    const bool rat = (variant == LJCL2_VAR_RATIONAL);
    // instantiated combinations only: full list and SMEM without pack/dual
    if (use_full) {
      if (rat) compact_compute<true,true,false>();
      else     compact_compute<true,false,false>();
    } else if (compact_smem) {
      if (rat) compact_compute<false,true,true>();
      else     compact_compute<false,false,true>();
    } else if (rat) {
      if (dual) { if (pack_on) compact_compute<false,true,false,true,true>();
                  else         compact_compute<false,true,false,false,true>(); }
      else      { if (pack_on) compact_compute<false,true,false,true,false>();
                  else         compact_compute<false,true,false,false,false>(); }
    } else {
      if (dual) { if (pack_on) compact_compute<false,false,false,true,true>();
                  else         compact_compute<false,false,false,false,true>(); }
      else      { if (pack_on) compact_compute<false,false,false,true,false>();
                  else         compact_compute<false,false,false,false,false>(); }
    }
    if (dual && dual_check && (this->update->ntimestep % dual_check) == 0) dual_verify();
    return;
  }

  if (dual) {
    if (ballot_on) {
      if (pack_on) flat_compute<false,true,true,LJCL2_VAR_NORMAL,true>();
      else         flat_compute<false,true,false,LJCL2_VAR_NORMAL,true>();
    } else {
      if (pack_on) flat_compute<false,true,true>();
      else         flat_compute<false,true,false>();
    }
  } else if (use_full) {
    if (variant == LJCL2_VAR_RATIONAL) flat_compute<true,false,false,LJCL2_VAR_RATIONAL>();
    else if (pack_on) flat_compute<true,false,true>();
    else              flat_compute<true,false,false>();
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
      case LJCL2_VAR_RATIONAL: flat_compute<false,false,false,LJCL2_VAR_RATIONAL>(); break;
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
template<bool FULL, bool DUAL, bool PACK, int VAR, bool BALLOT>
void PairLJCutCoulLong2Kokkos<DeviceType>::flat_compute()
{
  auto* k_list = static_cast<NeighListKokkos<DeviceType>*>(this->list);

  LJCL2Force<DeviceType,FULL,DUAL,PACK,VAR,BALLOT> ff;
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

/* ----------------------------------------------------------------------
   packed per-type-pair table for the compact kernel's shared-memory copy:
   [cutsq | cut_ljsq | cut_coulsq | lj1 | lj2], each section (MAX+1)^2 long
   and indexed i*(ntypes+1)+j.  Rebuilt on
   the host every call (a few hundred values) but uploaded only when it has
   changed, so a run without fix adapt pays one upload in total.
------------------------------------------------------------------------- */

template<class DeviceType>
void PairLJCutCoulLong2Kokkos<DeviceType>::tab_refresh()
{
  // stride must match LJCL2Compact::TSTRIDE; ntypes <= MAX_TYPES_STACKPARAMS
  // is guaranteed by the caller (larger systems take the base-class path)
  const int ntp1 = this->atom->ntypes + 1;
  const int ntab = (MAX_TYPES_STACKPARAMS+1)*(MAX_TYPES_STACKPARAMS+1);
  std::vector<KK_FLOAT> t(5*ntab, static_cast<KK_FLOAT>(0.0));
  for (int i = 1; i < ntp1; i++)
    for (int j = 1; j < ntp1; j++) {
      const int k = i*ntp1 + j;
      t[0*ntab + k] = this->m_cutsq[i][j];
      t[1*ntab + k] = this->m_cut_ljsq[i][j];
      t[2*ntab + k] = this->m_cut_coulsq[i][j];
      t[3*ntab + k] = this->m_params[i][j].lj1;
      t[4*ntab + k] = this->m_params[i][j].lj2;
    }
  if (t == h_tab) return;
  h_tab = t;
  d_tab = Kokkos::View<KK_FLOAT*, DeviceType>("ljcl2:tab", t.size());
  auto h = Kokkos::create_mirror_view(d_tab);
  for (size_t k = 0; k < t.size(); k++) h(k) = t[k];
  Kokkos::deep_copy(d_tab, h);
}

/* ----------------------------------------------------------------------
   compact kernel launch: one warp per atom (vector_length 32), warps_per_team
   warps per team from pair/team/size (default 128 threads = 4 warps), each
   warp walking compact_apw consecutive atoms of d_ilist.
------------------------------------------------------------------------- */

template<class DeviceType>
template<bool FULL, bool RAT, bool SMEM, bool PACK, bool DUAL>
void PairLJCutCoulLong2Kokkos<DeviceType>::compact_compute()
{
  auto* k_list = static_cast<NeighListKokkos<DeviceType>*>(this->list);
  if (SMEM) tab_refresh();

  using FunctorT = LJCL2Compact<DeviceType,FULL,RAT,SMEM,PACK,DUAL>;
  FunctorT ff;
  fill_common(ff);
  ff.d_numneigh = k_list->d_numneigh;
  ff.d_neighbors = k_list->d_neighbors;
  ff.d_neighbors_t = k_list->d_neighbors_transpose;
  ff.use_transpose = this->lmp->kokkos->neigh_transpose;
  ff.d_tab = d_tab;
  ff.ntp1 = this->atom->ntypes + 1;
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

  auto *kk = this->lmp->kokkos;
  const int team_threads = kk->pair_team_size_set ? kk->pair_team_size : 128;
  const int nwarp = (team_threads >= 32) ? team_threads/32 : 1;
  ff.warps_per_team = nwarp;
  ff.atoms_per_warp = compact_apw;
  const int atoms_per_team = nwarp*compact_apw;
  const int nteams = (k_list->inum + atoms_per_team - 1) / atoms_per_team;

  using policy_t = Kokkos::TeamPolicy<DeviceType, Kokkos::IndexType<int>>;
  policy_t policy(nteams > 0 ? nteams : 1, nwarp, 32);
  policy.set_scratch_size(0, Kokkos::PerTeam(FunctorT::shmem_bytes(nwarp)));
  Kokkos::parallel_for("PairLJCutCoulLong2::force_compact", policy, ff);
}

/* ----------------------------------------------------------------------
   tiled kernel: build the j-sets and the 16-bit local rows (once per
   reneighbor).  The j-set capacity self-tunes: the kernel always reports the
   largest j-set it met, so one retry converges.  One host readback per
   reneighbor, none per step.
------------------------------------------------------------------------- */

template<class DeviceType>
void PairLJCutCoulLong2Kokkos<DeviceType>::local_build()
{
#if defined(LJCL2_HAVE_LOCAL)
  if constexpr (std::is_same<typename DeviceType::execution_space, Kokkos::Cuda>::value) {
    auto* k_list = static_cast<NeighListKokkos<DeviceType>*>(this->list);
    const int inum = k_list->inum;
    const int ni = local_ni;
    const int ngroup = (inum + ni - 1) / ni;
    const int nrow = static_cast<int>(k_list->d_numneigh.extent(0));
    const int maxn = (k_list->maxneighs + 7) & ~7;   // 16-byte aligned rows
    const int nslot = (ni <= 32) ? 8192 : 16384;
    local_ngroup = ngroup;
    if (local_jcap < 256) local_jcap = 2048;

    if ((int) d_lnbr.extent(0) != nrow || (int) d_lnbr.extent(1) != maxn) {
      d_lnbr = decltype(d_lnbr)();
      d_lnbr = decltype(d_lnbr)(Kokkos::view_alloc("ljcl2:lnbr", Kokkos::WithoutInitializing),
                                nrow, maxn);
    }
    if ((int) d_lm.extent(0) != ngroup) d_lm = Kokkos::View<int*, DeviceType>("ljcl2:lm", ngroup);
    if (dual_on == 1) {
      if ((int) d_linr.extent(0) != nrow || (int) d_linr.extent(1) != maxn) {
        d_linr = decltype(d_linr)();
        d_linr = decltype(d_linr)(Kokkos::view_alloc("ljcl2:linr", Kokkos::WithoutInitializing),
                                  nrow, maxn);
      }
      if ((int) d_linnum.extent(0) != nrow)
        d_linnum = Kokkos::View<int*, DeviceType>("ljcl2:linnum", nrow);
    }
    if (d_lstat.extent(0) != 4) d_lstat = Kokkos::View<int*, DeviceType>("ljcl2:lstat", 4);
    auto h_stat = Kokkos::create_mirror_view(d_lstat);

    const bool tr = this->lmp->kokkos->neigh_transpose;
    LJCL2LocalBuildArgs b;
    b.ilist = k_list->d_ilist.data();
    b.numneigh = k_list->d_numneigh.data();
    if (tr) {
      b.nbr = k_list->d_neighbors_transpose.data();
      b.ns0 = k_list->d_neighbors_transpose.stride(0);
      b.ns1 = k_list->d_neighbors_transpose.stride(1);
    } else {
      b.nbr = k_list->d_neighbors.data();
      b.ns0 = k_list->d_neighbors.stride(0);
      b.ns1 = k_list->d_neighbors.stride(1);
    }
    b.lnbr = d_lnbr.data();
    b.lnbr_ld = d_lnbr.stride(0);
    b.inum = inum;
    b.ni = ni;
    b.nslot = nslot;

    const int bthr = 256;
    const int bshm = 2*nslot*static_cast<int>(sizeof(int));
    static bool battr = false;
    if (!battr) {
      cudaFuncSetAttribute(ljcl2_local_build, cudaFuncAttributeMaxDynamicSharedMemorySize, 2*16384*4);
      battr = true;
    }
    cudaStream_t st = Kokkos::Cuda().cuda_stream();

    int attempt = 0;
    for (; attempt < 4; attempt++) {
      if ((long long) d_ljset.extent(0) != (long long) ngroup*local_jcap) {
        d_ljset = Kokkos::View<int*, DeviceType>();
        d_ljset = Kokkos::View<int*, DeviceType>(
          Kokkos::view_alloc("ljcl2:ljset", Kokkos::WithoutInitializing), (size_t) ngroup*local_jcap);
      }
      Kokkos::deep_copy(d_lstat, 0);
      Kokkos::fence();
      b.ljset = d_ljset.data();
      b.lm = d_lm.data();
      b.stat = d_lstat.data();
      b.jcap = local_jcap;
      if (ngroup > 0) ljcl2_local_build<<<ngroup, bthr, bshm, st>>>(b);
      Kokkos::deep_copy(h_stat, d_lstat);   // fences
      if (!h_stat(1)) break;
      if (h_stat(0) > LJCL2_LMAX || h_stat(0) > nslot/2) {
        if (this->comm->me == 0)
          this->error->warning(FLERR, fmt::format("LMP_LJCL2_LOCAL: a group of {} atoms has {} "
                               "distinct neighbors (max {}) -- tiled kernel disabled; try a "
                               "smaller LMP_LJCL2_LNI", ni, h_stat(0), static_cast<int>(LJCL2_LMAX)));
        local_on = 0;
        return;
      }
      local_jcap = ((static_cast<int>(1.1*h_stat(0)) + 64) + 31) & ~31;
    }
    if (attempt == 4) this->error->one(FLERR, "LMP_LJCL2_LOCAL: j-set capacity did not converge");
    local_maxm = h_stat(0);

    // the force kernel's shared memory must fit the device's opt-in limit
    auto *kk = this->lmp->kokkos;
    int nthr = kk->pair_team_size_set ? kk->pair_team_size : 256;
    nthr = (nthr < 32) ? 32 : ((nthr > 256) ? 256 : (nthr & ~31));
    const int ms = (local_maxm + 31) & ~31;
    int optin = 0, dev = 0;
    cudaGetDevice(&dev);
    cudaDeviceGetAttribute(&optin, cudaDevAttrMaxSharedMemoryPerBlockOptin, dev);
    const int need = ljcl2_local_shmem(ms, nthr/32, local_on == 2);
    if (need > optin) {
      if (this->comm->me == 0)
        this->error->warning(FLERR, fmt::format("LMP_LJCL2_LOCAL: tile needs {} bytes of shared "
                             "memory (max {}) -- tiled kernel disabled; try a smaller "
                             "LMP_LJCL2_LNI", need, optin));
      local_on = 0;
      return;
    }
    if (local_stats && this->comm->me == 0)
      printf("LOCAL_STATS: step " BIGINT_FORMAT " ni %d groups %d mean M %.0f max M %d "
             "entries/M %.1f shmem %d B jcap %d retries %d\n", this->update->ntimestep, ni,
             ngroup, (double) h_stat(2)/ngroup, h_stat(0), (double) h_stat(3)/h_stat(2),
             need, local_jcap, attempt);
  }
#endif
}

/* ----------------------------------------------------------------------
   tiled kernel launch: one block per group, pair/team/size threads
   (default 256 = 8 warps), on the Kokkos stream so ordering is preserved.
------------------------------------------------------------------------- */

template<class DeviceType>
void PairLJCutCoulLong2Kokkos<DeviceType>::local_compute()
{
#if defined(LJCL2_HAVE_LOCAL)
  if constexpr (std::is_same<typename DeviceType::execution_space, Kokkos::Cuda>::value) {
    auto* k_list = static_cast<NeighListKokkos<DeviceType>*>(this->list);
    auto f = this->atomKK->k_f.template view<DeviceType>();

    LJCL2LocalArgs a;
    a.xq = reinterpret_cast<const float4 *>(d_xq.data());
    a.type8 = d_type8.data();
    a.ilist = k_list->d_ilist.data();
    a.numneigh = k_list->d_numneigh.data();
    a.lnbr = d_lnbr.data();
    a.lnbr_ld = d_lnbr.stride(0);
    a.ljset = d_ljset.data();
    a.lm = d_lm.data();
    const bool ldual = (dual_on == 1);
    a.linr = ldual ? d_linr.data() : nullptr;
    a.linnum = ldual ? d_linnum.data() : nullptr;
    a.ictl = ldual ? k_ictl.template view<DeviceType>().data() : nullptr;
    a.inner_cutsq = static_cast<float>(inner_cutsq);
    a.f = reinterpret_cast<float *>(f.data());
    a.fs0 = f.stride(0);
    a.fs1 = f.stride(1);
    a.jcap = local_jcap;
    a.inum = k_list->inum;
    a.ni = local_ni;
    a.ms = (local_maxm + 31) & ~31;
    const int nt = this->atom->ntypes;
    a.ntp1 = nt + 1;
    a.g_ewald = static_cast<float>(this->g_ewald);
    a.qqrd2e = static_cast<float>(this->force->qqrd2e);
    for (int k = 0; k < 4; k++) {
      a.special_lj[k] = static_cast<float>(this->force->special_lj[k]);
      a.special_coul[k] = static_cast<float>(this->force->special_coul[k]);
    }
    float cmax = 0.0f;
    for (int k = 0; k < LJCL2_NTAB; k++) a.lj1[k] = a.lj2[k] = a.cutlj[k] = a.cutcoul[k] = 0.0f;
    for (int i = 1; i <= nt; i++)
      for (int j = 1; j <= nt; j++) {
        const int t = i*a.ntp1 + j;
        a.lj1[t] = static_cast<float>(this->m_params[i][j].lj1);
        a.lj2[t] = static_cast<float>(this->m_params[i][j].lj2);
        a.cutlj[t] = static_cast<float>(this->m_cut_ljsq[i][j]);
        a.cutcoul[t] = static_cast<float>(this->m_cut_coulsq[i][j]);
        if (this->m_cutsq[i][j] > cmax) cmax = static_cast<float>(this->m_cutsq[i][j]);
      }
    a.cutsq_max = cmax;

    auto *kk = this->lmp->kokkos;
    int nthr = kk->pair_team_size_set ? kk->pair_team_size : 256;
    nthr = (nthr < 32) ? 32 : ((nthr > 256) ? 256 : (nthr & ~31));
    const int shm = ljcl2_local_shmem(a.ms, nthr/32, local_on == 2);
    cudaStream_t st = Kokkos::Cuda().cuda_stream();
    const bool rat = (variant == LJCL2_VAR_RATIONAL);
    // one attribute high-water mark per instantiation
    static int attr[24] = {0};
    auto launch = [&](void (*kern)(const LJCL2LocalArgs), int slot) {
      if (shm > attr[slot]) {
        cudaFuncSetAttribute(kern, cudaFuncAttributeMaxDynamicSharedMemorySize, shm);
        attr[slot] = shm;
      }
      if (local_ngroup > 0) kern<<<local_ngroup, nthr, shm, st>>>(a);
    };
    const int vs = (local_vec == 4) ? 2 : ((local_vec == 2) ? 1 : 0);
    const int slot = 12*(ldual ? 1 : 0) + 4*vs + ((local_on == 2) ? 2 : 0) + (rat ? 1 : 0);
    void (*tabk[24])(const LJCL2LocalArgs) = {
      ljcl2_local_force<false,0,1,false>, ljcl2_local_force<true,0,1,false>,
      ljcl2_local_force<false,1,1,false>, ljcl2_local_force<true,1,1,false>,
      ljcl2_local_force<false,0,2,false>, ljcl2_local_force<true,0,2,false>,
      ljcl2_local_force<false,1,2,false>, ljcl2_local_force<true,1,2,false>,
      ljcl2_local_force<false,0,4,false>, ljcl2_local_force<true,0,4,false>,
      ljcl2_local_force<false,1,4,false>, ljcl2_local_force<true,1,4,false>,
      ljcl2_local_force<false,0,1,true>, ljcl2_local_force<true,0,1,true>,
      ljcl2_local_force<false,1,1,true>, ljcl2_local_force<true,1,1,true>,
      ljcl2_local_force<false,0,2,true>, ljcl2_local_force<true,0,2,true>,
      ljcl2_local_force<false,1,2,true>, ljcl2_local_force<true,1,2,true>,
      ljcl2_local_force<false,0,4,true>, ljcl2_local_force<true,0,4,true>,
      ljcl2_local_force<false,1,4,true>, ljcl2_local_force<true,1,4,true>};
    launch(tabk[slot], slot);
    const cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess)
      this->error->one(FLERR, fmt::format("LMP_LJCL2_LOCAL: kernel launch failed: {}",
                                          cudaGetErrorString(err)));
  }
#endif
}

/* ----------------------------------------------------------------------
   cluster-pair list path: build (once per LAMMPS reneighbor)
------------------------------------------------------------------------- */

template<class DeviceType>
void PairLJCutCoulLong2Kokkos<DeviceType>::cl_build()
{
  using exec = typename DeviceType::execution_space;
  const int nlocal = this->atom->nlocal;
  const double *prd = this->domain->prd;
  const double *blo = this->domain->boxlo;
  const double rlist = this->cutforce + this->neighbor->skin;

  // ---- column grid: ~cubic 8-atom i-clusters
  const double density = nlocal / (prd[0]*prd[1]*prd[2]);
  const double w = cl_colscale * cbrt(8.0/density);
  cl_ncx = MAX(1, static_cast<int>(prd[0]/w));
  cl_ncy = MAX(1, static_cast<int>(prd[1]/w));
  cl_ncol = cl_ncx*cl_ncy;
  if (cl_ncol >= (1 << (32 - CL_ZBITS)))
    this->error->one(FLERR, "LMP_LJCL2_CLUSTER: too many grid columns");
  cl_colw[0] = prd[0]/cl_ncx;
  cl_colw[1] = prd[1]/cl_ncy;
  const int nx = static_cast<int>(ceil(rlist/cl_colw[0])) + 1;
  const int ny = static_cast<int>(ceil(rlist/cl_colw[1])) + 1;
  if ((2*nx+1 > cl_ncx) || (2*ny+1 > cl_ncy) || (prd[2] < 2.0*rlist + 2.0*w))
    this->error->one(FLERR, "LMP_LJCL2_CLUSTER: box too small for the cluster search");

  auto d_x    = this->atomKK->k_x.template view<DeviceType>();
  auto d_q    = this->atomKK->k_q.template view<DeviceType>();
  auto d_type = this->atomKK->k_type.template view<DeviceType>();

  // ---- sort local atoms by (column, z)
  if ((int) d_cl_key.extent(0) != nlocal) {
    d_cl_key = Kokkos::View<unsigned int*, DeviceType>("ljcl2:cl_key", nlocal);
    d_cl_perm = Kokkos::View<int*, DeviceType>("ljcl2:cl_perm", nlocal);
  }
  d_cl_colcnt = Kokkos::View<int*, DeviceType>("ljcl2:cl_colcnt", cl_ncol + 1);
  d_cl_colfirst = Kokkos::View<int*, DeviceType>("ljcl2:cl_colfirst", cl_ncol + 1);
  d_cl_colslot = Kokkos::View<int*, DeviceType>("ljcl2:cl_colslot", cl_ncol + 1);
  {
    ClKey<DeviceType> fk;
    fk.x = d_x; fk.key = d_cl_key; fk.perm = d_cl_perm; fk.colcnt = d_cl_colcnt;
    for (int d = 0; d < 3; d++) fk.blo[d] = blo[d];
    fk.icw[0] = 1.0/cl_colw[0]; fk.icw[1] = 1.0/cl_colw[1];
    fk.izq = static_cast<double>(1u << CL_ZBITS) / prd[2];
    fk.ncx = cl_ncx; fk.ncy = cl_ncy;
    Kokkos::parallel_for("PairLJCutCoulLong2::cl_key", Kokkos::RangePolicy<DeviceType>(0, nlocal), fk);
  }
  Kokkos::Experimental::sort_by_key(exec{}, d_cl_key, d_cl_perm);
  {
    auto cnt = d_cl_colcnt; auto first = d_cl_colfirst; auto slot = d_cl_colslot;
    const int ncol = cl_ncol;
    Kokkos::parallel_scan("PairLJCutCoulLong2::cl_scan1", Kokkos::RangePolicy<DeviceType>(0, ncol + 1),
      KOKKOS_LAMBDA(const int c, int &upd, const bool final) {
        if (final) first(c) = upd;
        if (c < ncol) upd += cnt(c);
      });
    Kokkos::parallel_scan("PairLJCutCoulLong2::cl_scan2", Kokkos::RangePolicy<DeviceType>(0, ncol + 1),
      KOKKOS_LAMBDA(const int c, int &upd, const bool final) {
        if (final) slot(c) = upd;
        if (c < ncol) upd += ((cnt(c) + CL_CI - 1)/CL_CI)*CL_CI;
      });
  }
  {
    auto h = Kokkos::create_mirror_view(Kokkos::subview(d_cl_colslot, cl_ncol));
    Kokkos::deep_copy(h, Kokkos::subview(d_cl_colslot, cl_ncol));
    cl_nslot = h();
  }
  cl_nicl = cl_nslot / CL_CI;
  const int njcl = cl_nslot / CL_CJ;

  if ((int) d_cl_atom.extent(0) < cl_nslot) {
    const int cap = cl_nslot + cl_nslot/10 + 64;
    d_cl_atom = Kokkos::View<int*, DeviceType>("ljcl2:cl_atom", cap);
    d_cl_xq = Kokkos::View<LJCL2XQ*, DeviceType>("ljcl2:cl_xq", cap);
    d_cl_xp = Kokkos::View<float*[4], Kokkos::LayoutRight, DeviceType>("ljcl2:cl_xp", cap);
    d_cl_t8 = Kokkos::View<unsigned char*, DeviceType>("ljcl2:cl_t8", cap);
    d_cl_f = Kokkos::View<KK_ACC_FLOAT*[3], Kokkos::LayoutRight, DeviceType>("ljcl2:cl_f", cap);
    d_cl_jbb = Kokkos::View<LJCL2BB*, DeviceType>("ljcl2:cl_jbb", cap/CL_CJ + 1);
    d_cl_icol = Kokkos::View<int*, DeviceType>("ljcl2:cl_icol", cap/CL_CI + 1);
    d_cl_nout = Kokkos::View<int*, DeviceType>("ljcl2:cl_nout", cap/CL_CI + 1);
    d_cl_nin = Kokkos::View<int*, DeviceType>("ljcl2:cl_nin", cap/CL_CI + 1);
    d_cl_out = decltype(d_cl_out)();
    d_cl_in = decltype(d_cl_in)();
  }
  if (d_cl_ctl.extent(0) == 0) {
    d_cl_ctl = Kokkos::View<int*, DeviceType>("ljcl2:cl_ctl", 8);
    d_cl_ev = Kokkos::View<double*, DeviceType>("ljcl2:cl_ev", 8);
  }
  Kokkos::deep_copy(Kokkos::subview(d_cl_atom, std::make_pair(0, cl_nslot)), -1);
  {
    ClSlot<DeviceType> fs;
    fs.key = d_cl_key; fs.perm = d_cl_perm; fs.colfirst = d_cl_colfirst;
    fs.colslot = d_cl_colslot; fs.atom = d_cl_atom;
    Kokkos::parallel_for("PairLJCutCoulLong2::cl_slot", Kokkos::RangePolicy<DeviceType>(0, nlocal), fs);
    ClIcol<DeviceType> fi;
    fi.colslot = d_cl_colslot; fi.icol = d_cl_icol;
    Kokkos::parallel_for("PairLJCutCoulLong2::cl_icol", Kokkos::RangePolicy<DeviceType>(0, cl_ncol), fi);
  }
  {
    ClGather<DeviceType> fg;
    fg.x = d_x; fg.q = d_q; fg.type = d_type; fg.atom = d_cl_atom; fg.xq = d_cl_xq; fg.t8 = d_cl_t8;
    Kokkos::parallel_for("PairLJCutCoulLong2::cl_gather", Kokkos::RangePolicy<DeviceType>(0, cl_nslot), fg);
    ClBBox<DeviceType> fb;
    fb.atom = d_cl_atom; fb.xq = d_cl_xq; fb.jbb = d_cl_jbb;
    Kokkos::parallel_for("PairLJCutCoulLong2::cl_bbox", Kokkos::RangePolicy<DeviceType>(0, njcl), fb);
  }

  // ---- tile search, with self-tuning capacities (one host readback per try)
  if (cl_maxout == 0) {
    const double pairs = CL_CI * 0.5 * density * 4.0/3.0 * MY_PI * rlist*rlist*rlist;
    cl_maxout = static_cast<int>(pairs / (32.0*0.25)) + 64;
  }
  if (cl_maxsc == 0) cl_maxsc = 16*cl_nicl + 1024;
  if ((int) d_cl_sc.extent(0) < cl_maxsc)
    d_cl_sc = Kokkos::View<unsigned int*[2], Kokkos::LayoutRight, DeviceType>("ljcl2:cl_sc", cl_maxsc);

  auto tag = this->atomKK->k_tag.template view<DeviceType>();
  auto special = this->atomKK->k_special.template view<DeviceType>();
  auto nspecial = this->atomKK->k_nspecial.template view<DeviceType>();
  auto h_ctl = Kokkos::create_mirror_view(d_cl_ctl);
  int attempt = 0;
  for (;; attempt++) {
    if ((int) d_cl_out.extent(0) < cl_nicl || (int) d_cl_out.extent(1) < cl_maxout) {
      d_cl_out = decltype(d_cl_out)();
      d_cl_in = decltype(d_cl_in)();
      d_cl_out = Kokkos::View<LJCL2CTile**, Kokkos::LayoutRight, DeviceType>(
        Kokkos::view_alloc("ljcl2:cl_out", Kokkos::WithoutInitializing), d_cl_nout.extent(0), cl_maxout);
      d_cl_in = Kokkos::View<LJCL2CTile**, Kokkos::LayoutRight, DeviceType>(
        Kokkos::view_alloc("ljcl2:cl_in", Kokkos::WithoutInitializing), d_cl_nout.extent(0), cl_maxout);
    }
    Kokkos::deep_copy(d_cl_ctl, 0);
    ClSearch<DeviceType> fs;
    fs.atom = d_cl_atom; fs.colslot = d_cl_colslot; fs.icol = d_cl_icol; fs.nout = d_cl_nout;
    fs.ctl = d_cl_ctl; fs.xq = d_cl_xq; fs.jbb = d_cl_jbb; fs.out = d_cl_out; fs.sc = d_cl_sc;
    fs.tag = tag; fs.special = special; fs.nspecial = nspecial;
    for (int d = 0; d < 3; d++) { fs.blo[d] = blo[d]; fs.prd[d] = prd[d]; }
    fs.cw[0] = cl_colw[0]; fs.cw[1] = cl_colw[1];
    fs.rlist2 = static_cast<float>(rlist*rlist);
    fs.ncx = cl_ncx; fs.ncy = cl_ncy; fs.nx = nx; fs.ny = ny;
    fs.maxout = cl_maxout; fs.maxsc = cl_maxsc;
    using policy_t = Kokkos::TeamPolicy<DeviceType, Kokkos::IndexType<int>>;
    policy_t policy(cl_nicl, 1, 32);
    policy.set_scratch_size(0, Kokkos::PerTeam(ClSearch<DeviceType>::shmem_bytes()));
    Kokkos::parallel_for("PairLJCutCoulLong2::cl_search", policy, fs);
    Kokkos::deep_copy(h_ctl, d_cl_ctl);
    bool retry = false;
    if (h_ctl(0) > 0) { cl_maxout = static_cast<int>(1.15*h_ctl(0)) + 32; retry = true; }
    if (h_ctl(3)) {
      cl_maxsc = static_cast<int>(1.5*h_ctl(2)) + 1024;
      d_cl_sc = Kokkos::View<unsigned int*[2], Kokkos::LayoutRight, DeviceType>("ljcl2:cl_sc", cl_maxsc);
      retry = true;
    }
    if (!retry) break;
    if (attempt == 4) this->error->one(FLERR, "LMP_LJCL2_CLUSTER: list capacity did not converge");
  }

  // force a prune on the first step after the build
  Kokkos::deep_copy(Kokkos::subview(d_cl_ctl, 4), 1);
  cl_built_step = this->neighbor->lastcall;

  if (cl_stats && this->comm->me == 0) {
    auto nout = d_cl_nout; auto out = d_cl_out;
    long long ntile = 0, npair = 0;
    Kokkos::parallel_reduce("PairLJCutCoulLong2::cl_stats", Kokkos::RangePolicy<DeviceType>(0, cl_nicl),
      KOKKOS_LAMBDA(const int ci, long long &t, long long &p) {
        t += nout(ci);
        for (int k = 0; k < nout(ci); k++) {
          unsigned int m = out(ci,k).mask;
          int c = 0; while (m) { m &= m - 1u; c++; }
          p += c;
        }
      }, ntile, npair);
    printf("CL_STATS: build step " BIGINT_FORMAT " cols %dx%d (%.2f x %.2f A) slots %d (pad %.1f%%) "
           "iclusters %d outer tiles/ic %.1f fill %.3f maxout %d specials %d tries %d\n",
           this->update->ntimestep, cl_ncx, cl_ncy, cl_colw[0], cl_colw[1], cl_nslot,
           100.0*(cl_nslot - nlocal)/cl_nslot, cl_nicl, (double) ntile/cl_nicl,
           (double) npair/(32.0*ntile), cl_maxout, h_ctl(2), attempt + 1);
  }
}

/* ----------------------------------------------------------------------
   cluster-pair list path: one force evaluation
------------------------------------------------------------------------- */

template<class DeviceType>
void PairLJCutCoulLong2Kokkos<DeviceType>::cl_step(int eflag_in, int vflag_in)
{
  // a call that was never joined (no VerletKokkos hook): complete it now
  if (cl_pending) finish_async_compute();

  const bool rebuild = (cl_built_step < 0) || (this->neighbor->lastcall > cl_built_step);
  if (rebuild) cl_build();

  const bool ev = (eflag_in || vflag_in);

  // overlap only force-only steps of a verlet/kk run, which joins the stream
  // through finish_async_compute() after the kspace compute
  bool async = false;
#if defined(KOKKOS_ENABLE_CUDA)
  async = cl_stream && !ev && !cl_stats && !this->update->setupflag &&
          (this->update->whichflag == 1) &&
          (strncmp(this->update->integrate_style, "verlet/kk", 9) == 0);
#endif
  const typename DeviceType::execution_space space =
    async ? cl_space : typename DeviceType::execution_space();
#if defined(KOKKOS_ENABLE_CUDA)
  if (async) {
    if constexpr (std::is_same<typename DeviceType::execution_space, Kokkos::Cuda>::value) {
      // everything queued so far on the default stream (integrate, comm, the
      // list build) must finish before the cluster kernels read x
      cudaEventRecord(cl_ev_in, Kokkos::Cuda().cuda_stream());
      cudaStreamWaitEvent(cl_stream, cl_ev_in, 0);
    }
  }
#endif

  auto d_x    = this->atomKK->k_x.template view<DeviceType>();
  auto d_q    = this->atomKK->k_q.template view<DeviceType>();
  auto d_type = this->atomKK->k_type.template view<DeviceType>();
  const double *prd = this->domain->prd;

  if (!rebuild) {
    ClGather<DeviceType> fg;
    fg.x = d_x; fg.q = d_q; fg.type = d_type; fg.atom = d_cl_atom; fg.xq = d_cl_xq; fg.t8 = d_cl_t8;
    Kokkos::parallel_for("PairLJCutCoulLong2::cl_gather", Kokkos::RangePolicy<DeviceType>(space, 0, cl_nslot), fg);
    ClDisp<DeviceType> fd;
    fd.atom = d_cl_atom; fd.ctl = d_cl_ctl; fd.xq = d_cl_xq; fd.xp = d_cl_xp;
    fd.half2 = static_cast<float>(0.25*cl_delta*cl_delta);
    Kokkos::parallel_for("PairLJCutCoulLong2::cl_disp", Kokkos::RangePolicy<DeviceType>(space, 0, cl_nslot), fd);
  }

  const int nwarp = 4;
  using policy_t = Kokkos::TeamPolicy<DeviceType, Kokkos::IndexType<int>>;
  const int nteams = (cl_nicl + nwarp - 1)/nwarp;
  {
    ClPrune<DeviceType> fp;
    fp.nout = d_cl_nout; fp.nin = d_cl_nin; fp.ctl = d_cl_ctl; fp.xq = d_cl_xq; fp.xp = d_cl_xp;
    fp.out = d_cl_out; fp.in = d_cl_in;
    const double rin = this->cutforce + cl_delta;
    fp.rin2 = static_cast<float>(rin*rin);
    fp.prd0 = prd[0]; fp.prd1 = prd[1]; fp.prd2 = prd[2];
    fp.nicl = cl_nicl; fp.warps_per_team = nwarp;
    policy_t policy(space, nteams, nwarp, 32);
    policy.set_scratch_size(0, Kokkos::PerTeam(ClPrune<DeviceType>::shmem_bytes(nwarp)));
    Kokkos::parallel_for("PairLJCutCoulLong2::cl_prune", policy, fp);
  }
  if (cl_stats) {
    auto h = Kokkos::create_mirror_view(Kokkos::subview(d_cl_ctl, 4));
    Kokkos::deep_copy(h, Kokkos::subview(d_cl_ctl, 4));
    cl_nprune += h();
    cl_nsteps++;
  }
  {
    ClFlagReset<DeviceType> fr; fr.ctl = d_cl_ctl;
    Kokkos::parallel_for("PairLJCutCoulLong2::cl_flag", Kokkos::RangePolicy<DeviceType>(space, 0, 1), fr);
  }

  if (variant < 0) {
    variant = LJCL2_VAR_NORMAL;
    if (const char *e = std::getenv("LMP_LJCL2_VARIANT")) variant = std::atoi(e);
  }
  const bool rat = (variant == LJCL2_VAR_RATIONAL);
  if (ev) Kokkos::deep_copy(d_cl_ev, 0.0);
  cl_launch_space = space;
  if (rat) { if (ev) cl_force<true,true>(); else cl_force<true,false>(); }
  else     { if (ev) cl_force<false,true>(); else cl_force<false,false>(); }
  cl_launch_space = typename DeviceType::execution_space();

#if defined(KOKKOS_ENABLE_CUDA)
  if (async) {
    cudaEventRecord(cl_ev_out, cl_stream);
    cl_pending = true;
    return;
  }
#endif
  cl_scatter();

  if (ev) {
    auto h = Kokkos::create_mirror_view(d_cl_ev);
    Kokkos::deep_copy(h, d_cl_ev);
    if (this->eflag_global) { this->eng_vdwl += h(0); this->eng_coul += h(1); }
    if (this->vflag_global) for (int k = 0; k < 6; k++) this->virial[k] += h(2+k);
  }

  if (cl_stats && this->comm->me == 0 && (this->update->ntimestep % 200 == 0) && cl_nsteps) {
    long long tin = 0;
    auto nin = d_cl_nin;
    Kokkos::parallel_reduce("PairLJCutCoulLong2::cl_stats_in", Kokkos::RangePolicy<DeviceType>(0, cl_nicl),
      KOKKOS_LAMBDA(const int ci, long long &t) { t += nin(ci); }, tin);
    printf("CL_STATS: step " BIGINT_FORMAT " inner tiles/ic %.1f prunes %d / %d steps\n",
           this->update->ntimestep, (double) tin/cl_nicl, cl_nprune, cl_nsteps);
  }
}

template<class DeviceType>
template<bool RAT, bool EV>
void PairLJCutCoulLong2Kokkos<DeviceType>::cl_force()
{
  // all type pairs share one cutoff (the common case): test against a scalar
  // before the j type is even loaded
  const int ntypes = this->atom->ntypes;
  bool unicut = true;
  for (int i = 1; i <= ntypes; i++)
    for (int j = 1; j <= ntypes; j++)
      if (this->m_cutsq[i][j] != this->m_cutsq[1][1]) unicut = false;
  if (unicut) cl_force_launch<RAT,EV,true>();
  else        cl_force_launch<RAT,EV,false>();
}

template<class DeviceType>
template<bool RAT, bool EV, bool UNICUT>
void PairLJCutCoulLong2Kokkos<DeviceType>::cl_force_launch()
{
  ClForce<DeviceType,RAT,EV,UNICUT> ff;
  fill_common(ff);
  ff.cutsq_uni = static_cast<KK_FLOAT>(this->m_cutsq[1][1]);
  ff.atom = d_cl_atom; ff.nin = d_cl_nin; ff.xq = d_cl_xq; ff.t8 = d_cl_t8; ff.in = d_cl_in;
  ff.sc = d_cl_sc; ff.fcl = d_cl_f; ff.ev = d_cl_ev; ff.ctl = d_cl_ctl;
  const double *prd = this->domain->prd;
  ff.prd0 = prd[0]; ff.prd1 = prd[1]; ff.prd2 = prd[2];
  const int nwarp = 4;
  ff.nicl = cl_nicl; ff.warps_per_team = nwarp;
  int nteams = (cl_nicl + nwarp - 1)/nwarp;
  if (cl_persist > 0 && cl_nsm > 0 && cl_persist*cl_nsm < nteams) {
    ff.persist = 1;
    nteams = cl_persist*cl_nsm;
  }
  using policy_t = Kokkos::TeamPolicy<DeviceType, Kokkos::IndexType<int>>;
  policy_t policy(cl_launch_space, nteams, nwarp, 32);
  policy.set_scratch_size(0, Kokkos::PerTeam(ClForce<DeviceType,RAT,EV,UNICUT>::shmem_bytes(nwarp)));
  Kokkos::parallel_for(EV ? "PairLJCutCoulLong2::cl_force_ev" : "PairLJCutCoulLong2::cl_force",
                       policy, ff);
}

namespace LAMMPS_NS {
template class PairLJCutCoulLong2Kokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairLJCutCoulLong2Kokkos<LMPHostType>;
#endif
}
