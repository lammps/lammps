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
#include "domain.h"
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
#include <cstring>
#include <type_traits>

using namespace LAMMPS_NS;
using namespace EwaldConst;
using MathConst::MY_PI;

/* ----------------------------------------------------------------------
   Force-only fast paths of lj/cut/coul/long2/kk.  Two independent paths:

   1. Compact kernel (CUDA, half list + newton on): one warp per atom i
      over the half neighbor list.  A test phase appends the pairs inside
      the cutoff to a per-warp shared-memory queue, and the force phase
      runs on fully converged warps, 32 queued pairs at a time.  j-forces
      go to f with atomics.

   2. Cluster-pair path (LMP_LJCL2_CLUSTER=1; CUDA, one MPI rank,
      orthogonal periodic box): a GROMACS-style list of 8x4 atom tiles that
      replaces the LAMMPS neighbor list, with an inner list pruned on the
      device.  Optionally runs on a second CUDA stream that overlaps with
      the bonded and kspace kernels (LMP_CL_ASYNC=1, LMP_CL_PERSIST=N).

   Both evaluate the same per-pair physics (LJCL2Common).  In single
   precision the real-space Ewald force uses a rational-polynomial fit, in
   double precision the analytical erfc approximation of the stock style.
   Energy/virial steps and every other configuration of the compact path
   fall back to PairLJCutCoulLongKokkos::compute().
------------------------------------------------------------------------- */

namespace {

// per-step packed copy of (x,y,z,q) as one 16-byte record plus a byte-wide
// type array, so each j gather is one 128-bit load and one byte load
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
// Rational-polynomial Ewald real-space force, used in single precision.
// Borrowed from GROMACS (pmecorrF in the nbnxm kernels): with z2 = (g_ewald r)^2, the erfc-screened Coulomb force is
//
//   F/r = qqrd2e qi qj [ 1/r^3 + g_ewald^3 P(z2) ],
//   P(z2) = -( erf(z)/z^3 - 2 exp(-z2)/(sqrt(pi) z^2) )  ~  N(z2)/D(z2)
//
// a 6th/4th-order rational fit with no exp, no erfc and one reciprocal.  The
// special-bond correction of the analytical path, -(1-factor_coul) qq/r^3,
// folds into the first term as factor_coul/r^3, so the branch disappears too.
//
// Accuracy, evaluated in fp32 against the exact erfc force over r in [1,12]
// at g_ewald = 0.1678: max |dF| per unit qq 2.5e-7, against
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

// Per-pair physics, shared by the compact and the cluster kernel
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

  // fpair for one pair already known to satisfy rsq < m_cutsq[itype][jtype]
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT pair_fpair(const KK_FLOAT rsq, const int itype, const int jtype,
                      const KK_FLOAT qiqj, const int sb) const {
    // one reciprocal square root supplies rinv, r2inv and r
    const KK_FLOAT rinv = Kokkos::rsqrt(rsq);
    const KK_FLOAT r2inv = rinv*rinv;
    KK_FLOAT fpair = static_cast<KK_FLOAT>(0.0);

    if (rsq < m_cut_ljsq[itype][jtype]) {
      const KK_FLOAT r6inv = r2inv*r2inv*r2inv;
      fpair += special_lj[sb] * r6inv *
        (m_params[itype][jtype].lj1*r6inv - m_params[itype][jtype].lj2) * r2inv;
    }

    if (rsq < m_cut_coulsq[itype][jtype]) {
      if (sizeof(KK_FLOAT) == sizeof(float)) {
        // rational fit; the special-bond correction folds into the 1/r^3 term
        const KK_FLOAT g2 = g_ewald*g_ewald;
        const KK_FLOAT pc = static_cast<KK_FLOAT>(ljcl2_pmecorrF(static_cast<float>(g2*rsq)));
        fpair += qqrd2e * qiqj * (special_coul[sb]*rinv*r2inv + g2*g_ewald*pc);
      } else {
        // analytical Ewald real-space correction (Abramowitz-Stegun 7.1.26)
        const KK_FLOAT r = rsq*rinv;
        const KK_FLOAT grij = g_ewald * r;
        const KK_FLOAT expm2 = Kokkos::exp(-grij*grij);
        const KK_FLOAT t = static_cast<KK_FLOAT>(1.0) /
          (static_cast<KK_FLOAT>(1.0) + static_cast<KK_FLOAT>(EWALD_P)*grij);
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
    }
    return fpair;
  }
};

/* ----------------------------------------------------------------------
   Warp-compacted force kernel, CUDA only, half list + newton on.

   One warp owns one atom i at a time and splits its neighbor row in two
   phases:
     1. test phase: 32 lanes take 32 consecutive entries, gather the packed
        (x,y,z,q) record and type of j and test the cutoff.  Passing pairs
        are appended to a per-warp shared-memory queue (ballot + popc
        prefix, no atomics), together with q(j) and type(j).
     2. force phase: whenever the queue holds >= 32 pairs, all 32 lanes pop
        one each and evaluate the pair fully converged.  The tail (< 32
        pairs) is drained once per atom.
   Only ~1/4 of the half-list entries are inside the cutoff with the usual
   skin, so without compaction the expensive pair body would run on every
   iteration with most lanes idle.  Since all queued pairs belong to the
   same i, the i-force stays in registers.  Each warp walks ATOMS_PER_WARP
   consecutive atoms of d_ilist.

   The neighbor list must be read through the LayoutRight transpose
   ('package kokkos neigh/transpose on'): with the LayoutLeft list the 32
   lanes read 32 different rows and the kernel is ~3x slower.
------------------------------------------------------------------------- */

template<class DeviceType>
struct LJCL2Compact : public LJCL2Common<DeviceType> {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef Kokkos::TeamPolicy<DeviceType, Kokkos::IndexType<int>> policy_type;
  typedef typename policy_type::member_type member_type;

  static constexpr int QCAP = 64;   // max queue fill is 31 + 32 = 63
  static constexpr int ATOMS_PER_WARP = 8;

  typename AT::t_int_1d d_numneigh;
  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_neighbors_2d_lr d_neighbors_t;
  Kokkos::View<const LJCL2XQ*, DeviceType> d_xq;
  Kokkos::View<const unsigned char*, DeviceType> d_type8;
  int warps_per_team;
  bool use_transpose;

  // per warp and queue slot: del x/y/z and q(j), jraw and type(j)
  KOKKOS_INLINE_FUNCTION
  static int shmem_bytes(int warps) {
    return warps*QCAP*(4*sizeof(KK_FLOAT) + 2*sizeof(int)) + 64;
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const member_type& team) const {
#if defined(__CUDA_ARCH__)
    const int nwarp = warps_per_team;
    char *sh = static_cast<char *>(team.team_shmem().get_shmem(shmem_bytes(nwarp)));
    const int w = team.team_rank();
    const int lane = threadIdx.x;
    KK_FLOAT *qx = reinterpret_cast<KK_FLOAT *>(sh) + w*4*QCAP;
    KK_FLOAT *qy = qx + QCAP;
    KK_FLOAT *qz = qy + QCAP;
    KK_FLOAT *qq = qz + QCAP;
    int *qj = reinterpret_cast<int *>(reinterpret_cast<KK_FLOAT *>(sh) + nwarp*4*QCAP) + w*2*QCAP;
    int *qt = qj + QCAP;

    const unsigned FULLMASK = 0xffffffffu;
    const unsigned lanelt = (1u << lane) - 1u;
    const int a0 = (team.league_rank()*nwarp + w)*ATOMS_PER_WARP;

    for (int a = 0; a < ATOMS_PER_WARP; a++) {
      const int ii = a0 + a;
      if (ii >= this->inum) break;                       // warp-uniform
      const int i = this->d_ilist(ii);
      const LJCL2XQ i4 = d_xq(i);
      const KK_FLOAT xtmp = i4.x, ytmp = i4.y, ztmp = i4.z, qtmp = i4.q;
      const int itype = static_cast<int>(d_type8(i));
      const int jnum = d_numneigh(i);

      KK_ACC_FLOAT fx = 0.0, fy = 0.0, fz = 0.0;
      int qn = 0;                                        // warp-uniform

      // force phase for queue slot p (all lanes of the warp, or the tail)
      auto physics = [&] (const int p) {
        const KK_FLOAT delx = qx[p], dely = qy[p], delz = qz[p];
        const int jraw = qj[p];
        const int j = jraw & NEIGHMASK;
        const int sb = (jraw >> SBBITS) & 3;
        const KK_FLOAT rsq = delx*delx + dely*dely + delz*delz;
        const KK_FLOAT fpair = this->pair_fpair(rsq, itype, qt[p], qtmp*qq[p], sb);
        fx += static_cast<KK_ACC_FLOAT>(delx*fpair);
        fy += static_cast<KK_ACC_FLOAT>(dely*fpair);
        fz += static_cast<KK_ACC_FLOAT>(delz*fpair);
        if (this->newton_pair || j < this->nlocal) {
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
          bool pass = false;
          KK_FLOAT delx = 0, dely = 0, delz = 0, qj_ = 0;
          int jraw = 0, jtype = 0;
          if (jj < jnum) {
            jraw = use_transpose ? d_neighbors_t(i,jj) : d_neighbors(i,jj);
            const LJCL2XQ j4 = d_xq(jraw & NEIGHMASK);
            delx = xtmp - j4.x; dely = ytmp - j4.y; delz = ztmp - j4.z; qj_ = j4.q;
            jtype = static_cast<int>(d_type8(jraw & NEIGHMASK));
            const KK_FLOAT rsq = delx*delx + dely*dely + delz*delz;
            pass = (rsq < this->m_cutsq[itype][jtype]);
          }
          const unsigned b = __ballot_sync(FULLMASK, pass);
          if (pass) {
            const int p = qn + __popc(b & lanelt);
            qx[p] = delx; qy[p] = dely; qz[p] = delz; qj[p] = jraw;
            qq[p] = qj_; qt[p] = jtype;
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
        Kokkos::atomic_add(&this->f(i,0), fx);
        Kokkos::atomic_add(&this->f(i,1), fy);
        Kokkos::atomic_add(&this->f(i,2), fz);
      }
    }
#else
    (void) team;
    Kokkos::abort("lj/cut/coul/long2/kk: the compact kernel requires the CUDA backend");
#endif
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

template<class DeviceType, bool EV, bool UNICUT>
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
              const KK_FLOAT fpair = this->pair_fpair(rsq, itype, jtype, qiqj, sb);
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

/* ----------------------------------------------------------------------
   the cluster-pair path is a prototype and selected by environment
   variables: LMP_LJCL2_CLUSTER=1 enables it, LMP_CL_DELTA sets the
   inner-list buffer (default 2.0 distance units), LMP_CL_ASYNC=1 runs the
   force-only steps on a second stream, LMP_CL_PERSIST=N launches N blocks
   per SM for the force kernel (needed for the overlap to take effect)
------------------------------------------------------------------------- */

template<class DeviceType>
void PairLJCutCoulLong2Kokkos<DeviceType>::init_style()
{
  PairLJCutCoulLongKokkos<DeviceType>::init_style();

  if (cl_on < 0) {
    const char *e = std::getenv("LMP_LJCL2_CLUSTER");
    cl_on = (e && utils::inumeric(FLERR, e, true, this->lmp)) ? 1 : 0;
    if ((e = std::getenv("LMP_CL_DELTA"))) cl_delta = utils::numeric(FLERR, e, true, this->lmp);
    if ((e = std::getenv("LMP_CL_ASYNC"))) cl_async = utils::inumeric(FLERR, e, true, this->lmp);
    if ((e = std::getenv("LMP_CL_PERSIST"))) cl_persist = utils::inumeric(FLERR, e, true, this->lmp);
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

/* ---------------------------------------------------------------------- */

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

  // The compact kernel handles force-only steps with a half list and newton on
  // on a CUDA GPU.  Everything else -- energy/virial steps for thermo, other
  // list settings, more than MAX_TYPES_STACKPARAMS atom types, other backends,
  // 'neigh/cluster yes' (which turns the request into a full-style list) --
  // uses the base-class kernel unchanged.
  bool cuda = false;
#if defined(KOKKOS_ENABLE_CUDA)
  cuda = std::is_same<typename DeviceType::execution_space, Kokkos::Cuda>::value;
#endif
  if (!cuda || eflag_in || vflag_in || (this->neighflag != HALFTHREAD) ||
      !this->force->newton_pair || (this->atom->ntypes > MAX_TYPES_STACKPARAMS) ||
      this->lmp->kokkos->neigh_cluster) {
    PairLJCutCoulLongKokkos<DeviceType>::compute(eflag_in, vflag_in);
    return;
  }

  this->eflag = eflag_in;
  this->vflag = vflag_in;
  this->ev_init(eflag_in, vflag_in, 0);

  this->atomKK->sync(this->execution_space, X_MASK | F_MASK | TYPE_MASK | Q_MASK);
  this->atomKK->modified(this->execution_space, F_MASK);

  pack_refresh();
  compact_compute();
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
   per-atom views, constants and stack parameters shared by both paths
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
   compact kernel launch: one warp per atom (vector length 32),
   warps_per_team warps per team from pair/team/size (default 128 threads
   = 4 warps), each warp walking ATOMS_PER_WARP consecutive atoms of d_ilist
------------------------------------------------------------------------- */

template<class DeviceType>
void PairLJCutCoulLong2Kokkos<DeviceType>::compact_compute()
{
  auto* k_list = static_cast<NeighListKokkos<DeviceType>*>(this->list);

  using FunctorT = LJCL2Compact<DeviceType>;
  FunctorT ff;
  fill_common(ff);
  ff.d_numneigh = k_list->d_numneigh;
  ff.d_neighbors = k_list->d_neighbors;
  ff.d_neighbors_t = k_list->d_neighbors_transpose;  // empty view unless neigh/transpose on
  ff.use_transpose = this->lmp->kokkos->neigh_transpose;
  ff.d_xq = d_xq;
  ff.d_type8 = d_type8;

  auto *kk = this->lmp->kokkos;
  const int team_threads = kk->pair_team_size_set ? kk->pair_team_size : 128;
  const int nwarp = (team_threads >= 32) ? team_threads/32 : 1;
  ff.warps_per_team = nwarp;
  const int atoms_per_team = nwarp*FunctorT::ATOMS_PER_WARP;
  const int nteams = (k_list->inum + atoms_per_team - 1) / atoms_per_team;

  using policy_t = Kokkos::TeamPolicy<DeviceType, Kokkos::IndexType<int>>;
  policy_t policy(nteams > 0 ? nteams : 1, nwarp, 32);
  policy.set_scratch_size(0, Kokkos::PerTeam(FunctorT::shmem_bytes(nwarp)));
  Kokkos::parallel_for("PairLJCutCoulLong2::force_compact", policy, ff);
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
  const double w = cbrt(8.0/density);
  cl_ncx = MAX(1, static_cast<int>(prd[0]/w));
  cl_ncy = MAX(1, static_cast<int>(prd[1]/w));
  cl_ncol = cl_ncx*cl_ncy;
  if (cl_ncol >= (1 << (32 - CL_ZBITS)))
    this->error->one(FLERR, "Pair style lj/cut/coul/long2/kk cluster list: too many grid columns");
  cl_colw[0] = prd[0]/cl_ncx;
  cl_colw[1] = prd[1]/cl_ncy;
  const int nx = static_cast<int>(ceil(rlist/cl_colw[0])) + 1;
  const int ny = static_cast<int>(ceil(rlist/cl_colw[1])) + 1;
  if ((2*nx+1 > cl_ncx) || (2*ny+1 > cl_ncy) || (prd[2] < 2.0*rlist + 2.0*w))
    this->error->one(FLERR, "Pair style lj/cut/coul/long2/kk cluster list: box too small for the cluster search");

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
    if (attempt == 4) this->error->one(FLERR, "Pair style lj/cut/coul/long2/kk cluster list: list capacity did not converge");
  }

  // force a prune on the first step after the build
  Kokkos::deep_copy(Kokkos::subview(d_cl_ctl, 4), 1);
  cl_built_step = this->neighbor->lastcall;

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
  async = cl_stream && !ev && !this->update->setupflag &&
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
  {
    ClFlagReset<DeviceType> fr; fr.ctl = d_cl_ctl;
    Kokkos::parallel_for("PairLJCutCoulLong2::cl_flag", Kokkos::RangePolicy<DeviceType>(space, 0, 1), fr);
  }

  if (ev) Kokkos::deep_copy(d_cl_ev, 0.0);
  cl_launch_space = space;
  if (ev) cl_force<true>();
  else cl_force<false>();
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
}

/* ----------------------------------------------------------------------
   cluster force kernel launch
------------------------------------------------------------------------- */

template<class DeviceType>
template<bool EV>
void PairLJCutCoulLong2Kokkos<DeviceType>::cl_force()
{
  // all type pairs share one cutoff (the common case): test against a scalar
  // before the j type is even loaded
  const int ntypes = this->atom->ntypes;
  bool unicut = true;
  for (int i = 1; i <= ntypes; i++)
    for (int j = 1; j <= ntypes; j++)
      if (this->m_cutsq[i][j] != this->m_cutsq[1][1]) unicut = false;
  if (unicut) cl_force_launch<EV,true>();
  else        cl_force_launch<EV,false>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<bool EV, bool UNICUT>
void PairLJCutCoulLong2Kokkos<DeviceType>::cl_force_launch()
{
  ClForce<DeviceType,EV,UNICUT> ff;
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
  policy.set_scratch_size(0, Kokkos::PerTeam(ClForce<DeviceType,EV,UNICUT>::shmem_bytes(nwarp)));
  Kokkos::parallel_for(EV ? "PairLJCutCoulLong2::cl_force_ev" : "PairLJCutCoulLong2::cl_force",
                       policy, ff);
}

namespace LAMMPS_NS {
template class PairLJCutCoulLong2Kokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairLJCutCoulLong2Kokkos<LMPHostType>;
#endif
}
