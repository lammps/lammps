// clang-format off
/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   aE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Stan Moore (SNL)
------------------------------------------------------------------------- */

#include "pair_pace_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "kokkos.h"
#include "math_const.h"
#include "memory_kokkos.h"
#include "neighbor_kokkos.h"
#include "neigh_request.h"

#include "ace-evaluator/ace_version.h"
#include "ace-evaluator/ace_radial.h"

#include "ace-evaluator/ace_c_basis.h"
#include "ace-evaluator/ace_evaluator.h"
#include "ace-evaluator/ace_recursive.h"

#include <cstring>

// prototype repeated from base class implementation
namespace LAMMPS_NS {
struct ACEImpl {
  ACECTildeBasisSet *basis_set;
  ACERecursiveEvaluator *ace;
};
}    // namespace LAMMPS_NS

using namespace LAMMPS_NS;
using namespace MathConst;

enum{FS,FS_SHIFTEDSCALED};

/* ----------------------------------------------------------------------
   Host-only batched helpers: 8 neighbors of one atom per call, one value
   per lane, everything laid out [quantity][lane] so the serial (l,m)
   recurrences run with independent lanes innermost and vectorize.  The
   pattern (lane-minor batches, tail lanes padded with valid data, plain
   C++ lane loops, vector clones on the definitions) follows the batched
   ChIMES evaluator and pair_snap_intel.
------------------------------------------------------------------------- */

namespace {

constexpr int PACE_VLEN = 8;

// stack-buffer capacities for the batched path; potentials that exceed them
// (far beyond any published ACE basis) fall back to the scalar per-neighbor
// path selected in init_style
constexpr int PACE_BATCH_NRL_MAX = 160;   // (lmax+1)*nradmax
constexpr int PACE_BATCH_SPH_MAX = 512;   // idx_sph_max*nradmax
constexpr int PACE_BATCH_NRB_MAX = 64;    // nradbase

// Compile the batched kernels once per instruction-set level and let the
// loader pick at run time, so the stock (generic -O2) build still runs the
// lane loops at AVX2+FMA width on machines that have it.  GCC wants the
// attribute on the definition only; the guard keeps device and non-x86
// compilers away from it.
#if defined(__x86_64__) && defined(__gnu_linux__) && !defined(__AVX2__) && \
    !defined(__CUDACC__) && !defined(__HIP_DEVICE_COMPILE__) && \
    defined(__GNUC__) && !defined(__clang__) && !defined(__INTEL_COMPILER)
#define PACE_VECTOR_CLONES __attribute__((target_clones("arch=x86-64-v3", "default")))
#else
#define PACE_VECTOR_CLONES
#endif

// Force contraction for one batch: the Ylm+dYlm recurrence of the scalar
// compute_derivative_one turned inside out.  The weights of the shared atom
// are broadcast scalars; every lane keeps its own force accumulator, so no
// cross-lane reduction is needed.  w is the atom's weights block viewed as
// re,im pairs; strides are in complex elements.
PACE_VECTOR_CLONES
void pace_batched_derivative(const int lmax, const int nradmax, const int nradbase,
                             const KK_FLOAT *alm, const KK_FLOAT *blm,
                             const KK_FLOAT *cl, const KK_FLOAT *dl,
                             const KK_FLOAT *rx, const KK_FLOAT *ry, const KK_FLOAT *rz,
                             const KK_FLOAT *rinv,
                             const KK_FLOAT *fr_b, const KK_FLOAT *dfr_b,
                             const KK_FLOAT *dgr_b,
                             const KK_FLOAT *w, const int w_ss, const int w_sn,
                             const KK_FLOAT *w_r1, const int wr1_sn,
                             KK_FLOAT *fx, KK_FLOAT *fy, KK_FLOAT *fz)
{
  constexpr int V = PACE_VLEN;

  for (int lane = 0; lane < V; lane++) fx[lane] = fy[lane] = fz[lane] = 0.0;

  for (int n = 0; n < nradbase; ++n) {
    const KK_FLOAT wn = w_r1[n * wr1_sn] * Y00;
    for (int lane = 0; lane < V; lane++) {
      const KK_FLOAT DGR = dgr_b[n * V + lane] * wn;
      fx[lane] += DGR * rx[lane];
      fy[lane] += DGR * ry[lane];
      fz[lane] += DGR * rz[lane];
    }
  }

  KK_FLOAT plm[V], plm1[V], plm2[V];
  KK_FLOAT dplm[V], dplm1[V], dplm2[V];

  int idx_sph = 0;

  // m = 0: ylm and dylm are real
  for (int l = 0; l <= lmax; l++) {
    if (l == 0) {
      for (int lane = 0; lane < V; lane++) { plm[lane] = Y00; dplm[lane] = 0.0; }
    } else if (l == 1) {
      for (int lane = 0; lane < V; lane++) {
        plm[lane] = Y00 * sq3 * rz[lane];
        dplm[lane] = Y00 * sq3;
      }
    } else {
      const KK_FLOAT al = alm[idx_sph], bl = blm[idx_sph];
      for (int lane = 0; lane < V; lane++) {
        plm[lane] = al * (rz[lane] * plm1[lane] + bl * plm2[lane]);
        dplm[lane] = al * (plm1[lane] + rz[lane] * dplm1[lane] + bl * dplm2[lane]);
      }
    }

    // ylm and dylm are real for m = 0 and depend only on l
    KK_FLOAT d0[V], d1[V], d2[V];
    for (int lane = 0; lane < V; lane++) {
      const KK_FLOAT rdy = dplm[lane] * rz[lane];
      d0[lane] = -rdy * rx[lane];
      d1[lane] = -rdy * ry[lane];
      d2[lane] = dplm[lane] - rdy * rz[lane];
    }

    const KK_FLOAT *fr_l = fr_b + l * nradmax * V;
    const KK_FLOAT *dfr_l = dfr_b + l * nradmax * V;
    const KK_FLOAT *w_l = w + 2 * (idx_sph * w_ss);
    for (int n = 0; n < nradmax; n++) {
      const KK_FLOAT wre = w_l[2 * (n * w_sn)];
      for (int lane = 0; lane < V; lane++) {
        const KK_FLOAT R_over_r = fr_l[n * V + lane] * rinv[lane];
        const KK_FLOAT DR = dfr_l[n * V + lane];
        const KK_FLOAT ydr = plm[lane] * DR;
        fx[lane] += wre * (ydr * rx[lane] + d0[lane] * R_over_r);
        fy[lane] += wre * (ydr * ry[lane] + d1[lane] * R_over_r);
        fz[lane] += wre * (ydr * rz[lane] + d2[lane] * R_over_r);
      }
    }

    for (int lane = 0; lane < V; lane++) {
      plm2[lane] = plm1[lane]; dplm2[lane] = dplm1[lane];
      plm1[lane] = plm[lane]; dplm1[lane] = dplm[lane];
    }
    idx_sph++;
  }

  // m = 1
  for (int lane = 0; lane < V; lane++)
    plm1[lane] = plm2[lane] = dplm1[lane] = dplm2[lane] = 0.0;
  for (int l = 1; l <= lmax; l++) {
    if (l == 1) {
      for (int lane = 0; lane < V; lane++) { plm[lane] = -sq3o2 * Y00; dplm[lane] = 0.0; }
    } else if (l == 2) {
      const KK_FLOAT d2 = dl[l];
      for (int lane = 0; lane < V; lane++) {
        const KK_FLOAT t = d2 * plm1[lane];
        plm[lane] = t * rz[lane];
        dplm[lane] = t;
      }
    } else {
      const KK_FLOAT al = alm[idx_sph], bl = blm[idx_sph];
      for (int lane = 0; lane < V; lane++) {
        plm[lane] = al * (rz[lane] * plm1[lane] + bl * plm2[lane]);
        dplm[lane] = al * (plm1[lane] + rz[lane] * dplm1[lane] + bl * dplm2[lane]);
      }
    }

    // ylm and dylm depend only on (l, m): computed once per l, reused for
    // every n (as in the scalar kernel)
    KK_FLOAT ylm_re[V], ylm_im[V];
    KK_FLOAT d0_re[V], d0_im[V], d1_re[V], d1_im[V], d2_re[V], d2_im[V];
    for (int lane = 0; lane < V; lane++) {
      ylm_re[lane] = rx[lane] * plm[lane];
      ylm_im[lane] = ry[lane] * plm[lane];
      // dyx = (plm, 0), dyy = (0, plm), dyz = phase * dplm
      const KK_FLOAT dyz_re = rx[lane] * dplm[lane];
      const KK_FLOAT dyz_im = ry[lane] * dplm[lane];
      const KK_FLOAT rdy_re = rx[lane] * plm[lane] + rz[lane] * dyz_re;
      const KK_FLOAT rdy_im = ry[lane] * plm[lane] + rz[lane] * dyz_im;
      d0_re[lane] = plm[lane] - rdy_re * rx[lane];
      d0_im[lane] = -rdy_im * rx[lane];
      d1_re[lane] = -rdy_re * ry[lane];
      d1_im[lane] = plm[lane] - rdy_im * ry[lane];
      d2_re[lane] = dyz_re - rdy_re * rz[lane];
      d2_im[lane] = dyz_im - rdy_im * rz[lane];
    }

    const KK_FLOAT *fr_l = fr_b + l * nradmax * V;
    const KK_FLOAT *dfr_l = dfr_b + l * nradmax * V;
    const KK_FLOAT *w_l = w + 2 * (idx_sph * w_ss);
    for (int n = 0; n < nradmax; n++) {
      const KK_FLOAT wre = 2.0 * w_l[2 * (n * w_sn)];
      const KK_FLOAT wim = 2.0 * w_l[2 * (n * w_sn) + 1];
      for (int lane = 0; lane < V; lane++) {
        const KK_FLOAT R_over_r = fr_l[n * V + lane] * rinv[lane];
        const KK_FLOAT DR = dfr_l[n * V + lane];
        const KK_FLOAT ydr_re = ylm_re[lane] * DR;
        const KK_FLOAT ydr_im = ylm_im[lane] * DR;
        const KK_FLOAT gx_re = ydr_re * rx[lane] + d0_re[lane] * R_over_r;
        const KK_FLOAT gx_im = ydr_im * rx[lane] + d0_im[lane] * R_over_r;
        const KK_FLOAT gy_re = ydr_re * ry[lane] + d1_re[lane] * R_over_r;
        const KK_FLOAT gy_im = ydr_im * ry[lane] + d1_im[lane] * R_over_r;
        const KK_FLOAT gz_re = ydr_re * rz[lane] + d2_re[lane] * R_over_r;
        const KK_FLOAT gz_im = ydr_im * rz[lane] + d2_im[lane] * R_over_r;
        fx[lane] += wre * gx_re - wim * gx_im;
        fy[lane] += wre * gy_re - wim * gy_im;
        fz[lane] += wre * gz_re - wim * gz_im;
      }
    }

    for (int lane = 0; lane < V; lane++) {
      plm2[lane] = plm1[lane]; dplm2[lane] = dplm1[lane];
      plm1[lane] = plm[lane]; dplm1[lane] = dplm[lane];
    }
    idx_sph++;
  }

  // m > 1
  for (int lane = 0; lane < V; lane++)
    plm1[lane] = plm2[lane] = dplm1[lane] = dplm2[lane] = 0.0;
  KK_FLOAT plm_mm1_mm1 = -sq3o2 * Y00;
  KK_FLOAT phasem_re[V], phasem_im[V];
  for (int lane = 0; lane < V; lane++) {
    phasem_re[lane] = rx[lane];
    phasem_im[lane] = ry[lane];
  }
  for (int m = 2; m <= lmax; m++) {
    KK_FLOAT mph_re[V], mph_im[V];
    for (int lane = 0; lane < V; lane++) {
      mph_re[lane] = phasem_re[lane] * (KK_FLOAT) m;
      mph_im[lane] = phasem_im[lane] * (KK_FLOAT) m;
      const KK_FLOAT pre = phasem_re[lane] * rx[lane] - phasem_im[lane] * ry[lane];
      const KK_FLOAT pim = phasem_re[lane] * ry[lane] + phasem_im[lane] * rx[lane];
      phasem_re[lane] = pre;
      phasem_im[lane] = pim;
    }

    for (int l = m; l <= lmax; l++) {
      if (l == m) {
        const KK_FLOAT seed = cl[l] * plm_mm1_mm1;
        plm_mm1_mm1 = seed;
        for (int lane = 0; lane < V; lane++) { plm[lane] = seed; dplm[lane] = 0.0; }
      } else if (l == (m + 1)) {
        const KK_FLOAT t = dl[l] * plm_mm1_mm1;
        for (int lane = 0; lane < V; lane++) {
          plm[lane] = t * rz[lane];
          dplm[lane] = t;
        }
      } else {
        const KK_FLOAT al = alm[idx_sph], bl = blm[idx_sph];
        for (int lane = 0; lane < V; lane++) {
          plm[lane] = al * (rz[lane] * plm1[lane] + bl * plm2[lane]);
          dplm[lane] = al * (plm1[lane] + rz[lane] * dplm1[lane] + bl * dplm2[lane]);
        }
      }

      KK_FLOAT ylm_re[V], ylm_im[V];
      KK_FLOAT d0_re[V], d0_im[V], d1_re[V], d1_im[V], d2_re[V], d2_im[V];
      for (int lane = 0; lane < V; lane++) {
        ylm_re[lane] = phasem_re[lane] * plm[lane];
        ylm_im[lane] = phasem_im[lane] * plm[lane];
        // dyx = mphasem1 * plm, dyy = i * dyx, dyz = phasem * dplm
        const KK_FLOAT dyx_re = mph_re[lane] * plm[lane];
        const KK_FLOAT dyx_im = mph_im[lane] * plm[lane];
        const KK_FLOAT dyy_re = -dyx_im;
        const KK_FLOAT dyy_im = dyx_re;
        const KK_FLOAT dyz_re = phasem_re[lane] * dplm[lane];
        const KK_FLOAT dyz_im = phasem_im[lane] * dplm[lane];
        const KK_FLOAT rdy_re = rx[lane] * dyx_re + ry[lane] * dyy_re + rz[lane] * dyz_re;
        const KK_FLOAT rdy_im = rx[lane] * dyx_im + ry[lane] * dyy_im + rz[lane] * dyz_im;
        d0_re[lane] = dyx_re - rdy_re * rx[lane];
        d0_im[lane] = dyx_im - rdy_im * rx[lane];
        d1_re[lane] = dyy_re - rdy_re * ry[lane];
        d1_im[lane] = dyy_im - rdy_im * ry[lane];
        d2_re[lane] = dyz_re - rdy_re * rz[lane];
        d2_im[lane] = dyz_im - rdy_im * rz[lane];
      }

      const KK_FLOAT *fr_l = fr_b + l * nradmax * V;
      const KK_FLOAT *dfr_l = dfr_b + l * nradmax * V;
      const KK_FLOAT *w_l = w + 2 * (idx_sph * w_ss);
      for (int n = 0; n < nradmax; n++) {
        const KK_FLOAT wre = 2.0 * w_l[2 * (n * w_sn)];
        const KK_FLOAT wim = 2.0 * w_l[2 * (n * w_sn) + 1];
        for (int lane = 0; lane < V; lane++) {
          const KK_FLOAT R_over_r = fr_l[n * V + lane] * rinv[lane];
          const KK_FLOAT DR = dfr_l[n * V + lane];
          const KK_FLOAT ydr_re = ylm_re[lane] * DR;
          const KK_FLOAT ydr_im = ylm_im[lane] * DR;
          const KK_FLOAT gx_re = ydr_re * rx[lane] + d0_re[lane] * R_over_r;
          const KK_FLOAT gx_im = ydr_im * rx[lane] + d0_im[lane] * R_over_r;
          const KK_FLOAT gy_re = ydr_re * ry[lane] + d1_re[lane] * R_over_r;
          const KK_FLOAT gy_im = ydr_im * ry[lane] + d1_im[lane] * R_over_r;
          const KK_FLOAT gz_re = ydr_re * rz[lane] + d2_re[lane] * R_over_r;
          const KK_FLOAT gz_im = ydr_im * rz[lane] + d2_im[lane] * R_over_r;
          fx[lane] += wre * gx_re - wim * gx_im;
          fy[lane] += wre * gy_re - wim * gy_im;
          fz[lane] += wre * gz_re - wim * gz_im;
        }
      }

      for (int lane = 0; lane < V; lane++) {
        plm2[lane] = plm1[lane]; dplm2[lane] = dplm1[lane];
        plm1[lane] = plm[lane]; dplm1[lane] = dplm[lane];
      }
      idx_sph++;
    }
  }
}

}    // namespace



/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairPACEKokkos<DeviceType>::PairPACEKokkos(LAMMPS *lmp) : PairPACE(lmp)
{
  respa_enable = 0;

  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;

  host_fallback = 0;
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

template<class DeviceType>
PairPACEKokkos<DeviceType>::~PairPACEKokkos()
{
  if (copymode) return;

  // when init_style() fell back to the base class, its compute() is used and
  // eatom/vatom are plain arrays from Pair::ev_setup() that are freed in
  // ~Pair(); calling destroy_kokkos() on them would clear the pointers without
  // freeing them.  In every other case the KOKKOS compute() allocated them, so
  // they must be freed here -- otherwise ~Pair() frees a Kokkos allocation.
  if (!host_fallback) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->destroy_kokkos(k_vatom,vatom);
  }

  deallocate_views_of_views();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairPACEKokkos<DeviceType>::deallocate_views_of_views()
{
  // deallocate views of views in serial to prevent race conditions

  if (k_splines_gk.view_host().data()) {
    for (int i = 0; i < nelements; i++) {
      for (int j = 0; j < nelements; j++) {
        k_splines_gk.view_host()(i, j).deallocate();
        k_splines_rnl.view_host()(i, j).deallocate();
        k_splines_hc.view_host()(i, j).deallocate();
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairPACEKokkos<DeviceType>::grow(int natom, int maxneigh)
{
  auto basis_set = aceimpl->basis_set;

  // A_sph_re is allocated on the device path only, so it cannot be the
  // growth test: on a host backend its extent stays 0 and every call would
  // re-allocate (and re-zero, and re-fault) the whole set.  A_rank1 is
  // allocated by both paths.
  if ((int)A_rank1.extent(0) < natom) {

    // Dual layout: the host kernels use the interleaved-complex A/A_sph (re
    // and im in the same cache line, which is what one thread walking one
    // atom wants), the device kernels use the split re/im arrays (which is
    // what a warp scattering across atoms wants). Only one set is allocated.
    if constexpr (host_flag) {
      MemKK::realloc_kokkos(A_sph, "pace:A_sph", natom, nelements, idx_sph_max, nradmax + 1);
      MemKK::realloc_kokkos(A, "pace:A", natom, nelements, (lmax + 1) * (lmax + 1), nradmax + 1);
    } else {
      MemKK::realloc_kokkos(A_sph_re, "pace:A_sph_re", natom, nelements, idx_sph_max, nradmax + 1);
      MemKK::realloc_kokkos(A_sph_im, "pace:A_sph_im", natom, nelements, idx_sph_max, nradmax + 1);
    }
    MemKK::realloc_kokkos(A_rank1, "pace:A_rank1", natom, nelements, nradbase);


    MemKK::realloc_kokkos(e_atom, "pace:e_atom", natom);
    MemKK::realloc_kokkos(rhos, "pace:rhos", natom, basis_set->ndensitymax + 1); // +1 density for core repulsion
    MemKK::realloc_kokkos(dF_drho, "pace:dF_drho", natom, basis_set->ndensitymax + 1); // +1 density for core repulsion

    if constexpr (host_flag) {
      // one extra row on the host: the sink for (l,m) entries outside the
      // packed triangle, see d_idx_sph_cpu
      MemKK::realloc_kokkos(weights, "pace:weights", natom, nelements,
                            idx_sph_max + 1, nradmax + 1);
    } else {
      MemKK::realloc_kokkos(weights_re, "pace:weights_re", natom, nelements, idx_sph_max, nradmax + 1);
      MemKK::realloc_kokkos(weights_im, "pace:weights_im", natom, nelements, idx_sph_max, nradmax + 1);
    }
    MemKK::realloc_kokkos(weights_rank1, "pace:weights_rank1", natom, nelements, nradbase);

    // hard-core repulsion
    MemKK::realloc_kokkos(rho_core, "pace:rho_core", natom);
    MemKK::realloc_kokkos(dF_drho_core, "pace:dF_drho_core", natom);
    MemKK::realloc_kokkos(dF_dfcut, "pace:dF_dfcut", natom);
    MemKK::realloc_kokkos(d_d_min, "pace:r_min_pair", natom);
    MemKK::realloc_kokkos(d_jj_min, "pace:j_min_pair", natom);
    MemKK::realloc_kokkos(d_corerep, "pace:corerep", natom); // per-atom corerep

    // dB_flatten is a host-only array: the device path rebuilds the
    // leave-one-out products in ComputeWeights instead of storing them
    if constexpr (host_flag)
      // rank-compacted on the host, indexed through d_dB_off
      MemKK::realloc_kokkos(dB_flatten, "pace:dB_flatten", natom, dB_total_max, 1);
  }

  if constexpr (host_flag) {
    // the offset tables assume the trailing dimensions are packed
    if (((int)A.stride(2) != (int)A.extent(3)) ||
        ((int)A.stride(1) != (int)A.extent(2) * (int)A.extent(3)) ||
        ((int)weights.stride(2) != (int)weights.extent(3)) ||
        ((int)weights.stride(1) != (int)weights.extent(2) * (int)weights.extent(3)))
      error->all(FLERR, "Pair style pace/kk: unexpected padding in the "
                        "coefficient arrays on a CPU backend");
  }

  if (((int)fr.extent(0) < natom) || ((int)fr.extent(1) < maxneigh)) {

    // radial functions
    MemKK::realloc_kokkos(fr, "pace:fr", natom, maxneigh, (lmax + 1) * nradmax);
    MemKK::realloc_kokkos(dfr, "pace:dfr", natom, maxneigh, (lmax + 1) * nradmax);
    MemKK::realloc_kokkos(gr, "pace:gr", natom, maxneigh, nradbase);
    MemKK::realloc_kokkos(dgr, "pace:dgr", natom, maxneigh, nradbase);
    const int max_num_functions = MAX(nradbase, nradmax*(lmax + 1));
    MemKK::realloc_kokkos(d_values, "pace:d_values", natom, maxneigh, max_num_functions);
    MemKK::realloc_kokkos(d_derivatives, "pace:d_derivatives", natom, maxneigh, max_num_functions);

    // hard-core repulsion
    MemKK::realloc_kokkos(cr, "pace:cr", natom, maxneigh);
    MemKK::realloc_kokkos(dcr, "pace:dcr", natom, maxneigh);

    // short neigh list
    MemKK::realloc_kokkos(d_ncount, "pace:ncount", natom);
    MemKK::realloc_kokkos(d_mu, "pace:mu", natom, maxneigh);
    MemKK::realloc_kokkos(d_rhats, "pace:rhats", natom, maxneigh);
    MemKK::realloc_kokkos(d_rnorms, "pace:rnorms", natom, maxneigh);
    MemKK::realloc_kokkos(d_nearest, "pace:nearest", natom, maxneigh);

    MemKK::realloc_kokkos(f_ij, "pace:f_ij", natom, maxneigh);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairPACEKokkos<DeviceType>::copy_pertype()
{
  auto basis_set = aceimpl->basis_set;

  MemKK::realloc_kokkos(d_rho_core_cutoff, "pace:rho_core_cutoff", nelements);
  MemKK::realloc_kokkos(d_drho_core_cutoff, "pace:drho_core_cutoff", nelements);
  MemKK::realloc_kokkos(d_E0vals, "pace:E0vals", nelements);
  MemKK::realloc_kokkos(d_ndensity, "pace:ndensity", nelements);
  MemKK::realloc_kokkos(d_npoti, "pace:npoti", nelements);

  auto h_rho_core_cutoff = Kokkos::create_mirror_view(d_rho_core_cutoff);
  auto h_drho_core_cutoff = Kokkos::create_mirror_view(d_drho_core_cutoff);
  auto h_E0vals = Kokkos::create_mirror_view(d_E0vals);
  auto h_ndensity = Kokkos::create_mirror_view(d_ndensity);
  auto h_npoti = Kokkos::create_mirror_view(d_npoti);

  for (int n = 0; n < nelements; n++) {
    h_rho_core_cutoff[n] = basis_set->map_embedding_specifications.at(n).rho_core_cutoff;
    h_drho_core_cutoff[n] = basis_set->map_embedding_specifications.at(n).drho_core_cutoff;

    h_E0vals(n) = basis_set->E0vals(n);

    h_ndensity(n) = basis_set->map_embedding_specifications.at(n).ndensity;

    string npoti = basis_set->map_embedding_specifications.at(n).npoti;
    if (npoti == "FinnisSinclair")
      h_npoti(n) = FS;
    else if (npoti == "FinnisSinclairShiftedScaled")
      h_npoti(n) = FS_SHIFTEDSCALED;
  }

  Kokkos::deep_copy(d_rho_core_cutoff, h_rho_core_cutoff);
  Kokkos::deep_copy(d_drho_core_cutoff, h_drho_core_cutoff);
  Kokkos::deep_copy(d_E0vals, h_E0vals);
  Kokkos::deep_copy(d_ndensity, h_ndensity);
  Kokkos::deep_copy(d_npoti, h_npoti);

  MemKK::realloc_kokkos(d_wpre, "pace:wpre", nelements, basis_set->ndensitymax);
  MemKK::realloc_kokkos(d_mexp, "pace:mexp", nelements, basis_set->ndensitymax);

  auto h_wpre = Kokkos::create_mirror_view(d_wpre);
  auto h_mexp = Kokkos::create_mirror_view(d_mexp);

  for (int n = 0; n < nelements; n++) {
    const int ndensity = basis_set->map_embedding_specifications.at(n).ndensity;
    for (int p = 0; p < ndensity; p++) {
      h_wpre(n, p) = basis_set->map_embedding_specifications.at(n).FS_parameters.at(p * 2 + 0);
      h_mexp(n, p) = basis_set->map_embedding_specifications.at(n).FS_parameters.at(p * 2 + 1);
    }
  }

  Kokkos::deep_copy(d_wpre, h_wpre);
  Kokkos::deep_copy(d_mexp, h_mexp);

  // ZBL core-rep
  MemKK::realloc_kokkos(d_cut_in, "pace:d_cut_in", nelements, nelements);
  MemKK::realloc_kokkos(d_dcut_in, "pace:d_dcut_in", nelements, nelements);
  auto h_cut_in = Kokkos::create_mirror_view(d_cut_in);
  auto h_dcut_in = Kokkos::create_mirror_view(d_dcut_in);

  for (int mu_i = 0; mu_i < nelements; ++mu_i) {
    for (int mu_j = 0; mu_j < nelements; ++mu_j) {
      h_cut_in(mu_i,mu_j) = basis_set->map_bond_specifications.at({mu_i,mu_j}).rcut_in;
      h_dcut_in(mu_i,mu_j) = basis_set->map_bond_specifications.at({mu_i,mu_j}).dcut_in;
    }
  }
  Kokkos::deep_copy(d_cut_in, h_cut_in);
  Kokkos::deep_copy(d_dcut_in, h_dcut_in);

  is_zbl = basis_set->radial_functions->inner_cutoff_type == "zbl";
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairPACEKokkos<DeviceType>::copy_splines()
{
  auto basis_set = aceimpl->basis_set;

  deallocate_views_of_views();

  k_splines_gk = Kokkos::DualView<SplineInterpolatorKokkos**, DeviceType>("pace:splines_gk", nelements, nelements);
  k_splines_rnl = Kokkos::DualView<SplineInterpolatorKokkos**, DeviceType>("pace:splines_rnl", nelements, nelements);
  k_splines_hc = Kokkos::DualView<SplineInterpolatorKokkos**, DeviceType>("pace:splines_hc", nelements, nelements);

  ACERadialFunctions* radial_functions = dynamic_cast<ACERadialFunctions*>(basis_set->radial_functions);

  if (radial_functions == nullptr)
    error->all(FLERR,"Chosen radial basis style not supported by pair style pace/kk");

  for (int i = 0; i < nelements; i++) {
    for (int j = 0; j < nelements; j++) {
      k_splines_gk.view_host()(i, j) = radial_functions->splines_gk(i, j);
      k_splines_rnl.view_host()(i, j) = radial_functions->splines_rnl(i, j);
      k_splines_hc.view_host()(i, j) = radial_functions->splines_hc(i, j);
    }
  }

  k_splines_gk.modify_host();
  k_splines_rnl.modify_host();
  k_splines_hc.modify_host();

  k_splines_gk.sync_device();
  k_splines_rnl.sync_device();
  k_splines_hc.sync_device();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairPACEKokkos<DeviceType>::copy_tilde()
{
  auto basis_set = aceimpl->basis_set;

  // flatten loops, get per-element count and max

  idx_ms_combs_max = 0;
  int total_basis_size_max = 0;

  MemKK::realloc_kokkos(d_idx_ms_combs_count, "pace:idx_ms_combs_count", nelements);
  auto h_idx_ms_combs_count = Kokkos::create_mirror_view(d_idx_ms_combs_count);
  MemKK::realloc_kokkos(d_nms_rank1, "pace:nms_rank1", nelements);
  auto h_nms_rank1 = Kokkos::create_mirror_view(d_nms_rank1);

  for (int mu = 0; mu < nelements; mu++) {
    int idx_ms_combs = 0;
    const int total_basis_size_rank1 = basis_set->total_basis_size_rank1[mu];
    const int total_basis_size = basis_set->total_basis_size[mu];

    ACECTildeBasisFunction *basis = basis_set->basis[mu];

    // rank=1
    for (int func_rank1_ind = 0; func_rank1_ind < total_basis_size_rank1; ++func_rank1_ind)
      idx_ms_combs++;

    // rank > 1
    for (int idx_func = 0; idx_func < total_basis_size; ++idx_func) {
      ACECTildeBasisFunction *func = &basis[idx_func];

      // loop over {ms} combinations in sum
      for (int ms_ind = 0; ms_ind < func->num_ms_combs; ++ms_ind)
        idx_ms_combs++;
    }
    h_idx_ms_combs_count(mu) = idx_ms_combs;
    h_nms_rank1(mu) = total_basis_size_rank1;
    idx_ms_combs_max = MAX(idx_ms_combs_max, idx_ms_combs);
    total_basis_size_max = MAX(total_basis_size_max, total_basis_size_rank1 + total_basis_size);
  }

  Kokkos::deep_copy(d_idx_ms_combs_count, h_idx_ms_combs_count);
  Kokkos::deep_copy(d_nms_rank1, h_nms_rank1);

  MemKK::realloc_kokkos(d_rank, "pace:rank", nelements, total_basis_size_max);
  MemKK::realloc_kokkos(d_num_ms_combs, "pace:num_ms_combs", nelements, total_basis_size_max);
  MemKK::realloc_kokkos(d_idx_funcs, "pace:idx_func", nelements, idx_ms_combs_max);
  MemKK::realloc_kokkos(d_mus, "pace:mus", nelements, total_basis_size_max, basis_set->rankmax);
  MemKK::realloc_kokkos(d_ns, "pace:ns", nelements, total_basis_size_max, basis_set->rankmax);
  MemKK::realloc_kokkos(d_ls, "pace:ls", nelements, total_basis_size_max, basis_set->rankmax);
  MemKK::realloc_kokkos(d_ms_combs, "pace:ms_combs", nelements, idx_ms_combs_max, basis_set->rankmax);
  MemKK::realloc_kokkos(d_ctildes, "pace:ctildes", nelements, idx_ms_combs_max, basis_set->ndensitymax);

  auto h_rank = Kokkos::create_mirror_view(d_rank);
  auto h_num_ms_combs = Kokkos::create_mirror_view(d_num_ms_combs);
  auto h_idx_funcs = Kokkos::create_mirror_view(d_idx_funcs);
  auto h_mus = Kokkos::create_mirror_view(d_mus);
  auto h_ns = Kokkos::create_mirror_view(d_ns);
  auto h_ls = Kokkos::create_mirror_view(d_ls);
  auto h_ms_combs = Kokkos::create_mirror_view(d_ms_combs);
  auto h_ctildes = Kokkos::create_mirror_view(d_ctildes);

  // copy values on host

  for (int mu = 0; mu < nelements; mu++) {
    const int total_basis_size_rank1 = basis_set->total_basis_size_rank1[mu];
    const int total_basis_size = basis_set->total_basis_size[mu];

    ACECTildeBasisFunction *basis_rank1 = basis_set->basis_rank1[mu];
    ACECTildeBasisFunction *basis = basis_set->basis[mu];

    const int ndensity = basis_set->map_embedding_specifications.at(mu).ndensity;

    int idx_ms_combs = 0;

    // rank=1
    for (int idx_func = 0; idx_func < total_basis_size_rank1; ++idx_func) {
      ACECTildeBasisFunction *func = &basis_rank1[idx_func];
      h_rank(mu, idx_func) = 1;
      h_mus(mu, idx_func, 0) = func->mus[0];
      h_ns(mu, idx_func, 0) = func->ns[0];

      for (int p = 0; p < ndensity; ++p)
        h_ctildes(mu, idx_ms_combs, p) = func->ctildes[p];

      h_idx_funcs(mu, idx_ms_combs) = idx_func;
      idx_ms_combs++;
    }

    // rank > 1
    for (int idx_func = 0; idx_func < total_basis_size; ++idx_func) {
      ACECTildeBasisFunction *func = &basis[idx_func];
      // TODO: check if func->ctildes are zero, then skip

      const int idx_func_through = total_basis_size_rank1 + idx_func;

      const int rank = h_rank(mu, idx_func_through) = func->rank;
      h_num_ms_combs(mu, idx_func_through) = func->num_ms_combs;
      for (int t = 0; t < rank; t++) {
        h_mus(mu, idx_func_through, t) = func->mus[t];
        h_ns(mu, idx_func_through, t) = func->ns[t];
        h_ls(mu, idx_func_through, t) = func->ls[t];
      }

      // loop over {ms} combinations in sum
      for (int ms_ind = 0; ms_ind < func->num_ms_combs; ++ms_ind) {
        auto ms = &func->ms_combs[ms_ind * rank]; // current ms-combination (of length = rank)
        for (int t = 0; t < rank; t++)
          h_ms_combs(mu, idx_ms_combs, t) = ms[t];

        for (int p = 0; p < ndensity; ++p) {
          // real-part only multiplication
          h_ctildes(mu, idx_ms_combs, p) = func->ctildes[ms_ind * ndensity + p];
        }

        h_idx_funcs(mu, idx_ms_combs) = idx_func_through;
        idx_ms_combs++;
      }
    }
  }

  Kokkos::deep_copy(d_rank, h_rank);
  Kokkos::deep_copy(d_num_ms_combs, h_num_ms_combs);
  Kokkos::deep_copy(d_idx_funcs, h_idx_funcs);
  Kokkos::deep_copy(d_mus, h_mus);
  Kokkos::deep_copy(d_ns, h_ns);
  Kokkos::deep_copy(d_ls, h_ls);
  Kokkos::deep_copy(d_ms_combs, h_ms_combs);
  Kokkos::deep_copy(d_ctildes, h_ctildes);
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

template<class DeviceType>
void PairPACEKokkos<DeviceType>::init_style()
{
  // two cases have no KOKKOS kernels on a CPU backend.  On a GPU each of them
  // is an error (see compute()), but on a CPU the non-accelerated evaluator is
  // available and is what pace/kk used before it gained CPU kernels, so use it
  // rather than rejecting the potential.  host_fallback records the choice for
  // compute() and the destructor.

  if (host_flag && recursive) {
    if (comm->me == 0)
      error->warning(FLERR, "Pair style pace/kk has no KOKKOS implementation of the "
                            "recursive evaluator and falls back to the non-accelerated "
                            "one; use the 'product' keyword for the threaded KOKKOS "
                            "calculation");
    host_fallback = 1;
  } else if (host_flag && (aceimpl->basis_set->rankmax > MAX_RANK_CPU)) {
    // the CPU kernels hold the leave-one-out products of a basis function in
    // small fixed-size stack arrays, which caps the correlation order
    if (comm->me == 0)
      error->warning(FLERR, "Pair style pace/kk supports a maximum correlation order of "
                            "{} on CPU backends, but this potential uses {}; falling back "
                            "to the non-accelerated evaluator",
                     MAX_RANK_CPU, aceimpl->basis_set->rankmax);
    host_fallback = 1;
  }

  if (host_fallback) {
    PairPACE::init_style();
    return;
  }

  if (atom->tag_enable == 0) error->all(FLERR, "Pair style PACE requires atom IDs");
  if (force->newton_pair == 0) error->all(FLERR, "Pair style PACE requires newton pair on");

  // neighbor list request for KOKKOS

  neighflag = lmp->kokkos->neighflag;

  auto request = neighbor->add_request(this, NeighConst::REQ_FULL);
  request->set_kokkos_host(std::is_same_v<DeviceType,LMPHostType> &&
                           !std::is_same_v<DeviceType,LMPDeviceType>);
  request->set_kokkos_device(std::is_same_v<DeviceType,LMPDeviceType>);
  if (neighflag == FULL)
    error->all(FLERR,"Must use half neighbor list style with pair pace/kk");

  auto basis_set = aceimpl->basis_set;

  nelements = basis_set->nelements;
  lmax = basis_set->lmax;
  nradmax = basis_set->nradmax;
  nradbase = basis_set->nradbase;

  // spherical harmonics

  MemKK::realloc_kokkos(d_idx_sph, "pace:idx_sph", (lmax + 1) * (lmax + 1));
  MemKK::realloc_kokkos(d_idx_sph_cpu, "pace:idx_sph_cpu", (lmax + 1) * (lmax + 1));
  MemKK::realloc_kokkos(alm, "pace:alm", (lmax + 1) * (lmax + 1));
  MemKK::realloc_kokkos(blm, "pace:blm", (lmax + 1) * (lmax + 1));
  MemKK::realloc_kokkos(cl, "pace:cl", lmax + 1);
  MemKK::realloc_kokkos(dl, "pace:dl", lmax + 1);

  pre_compute_harmonics(lmax);
  copy_pertype();
  copy_splines();
  copy_tilde();

  if constexpr (host_flag) {
    build_cpu_offset_tables();
    use_batched_cpu = (((lmax + 1) * nradmax <= PACE_BATCH_NRL_MAX) &&
                       (idx_sph_max * nradmax <= PACE_BATCH_SPH_MAX) &&
                       (nradbase <= PACE_BATCH_NRB_MAX));
  }
}

/* ----------------------------------------------------------------------
   CPU only: precompute, for every (ms-combination, rank) term, the flat
   offset of the A entry it gathers and of the two weights entries it
   scatters into.  The inner loops of the density and weights passes then do
   one integer load each instead of four loads out of mus/ns/ls/ms_combs plus
   the index arithmetic to rebuild the subscript.
------------------------------------------------------------------------- */

template<class DeviceType>
void PairPACEKokkos<DeviceType>::build_cpu_offset_tables()
{
  const int rankmax = (int) d_ms_combs.extent(2);
  const int A_l = (lmax + 1) * (lmax + 1);
  const int A_n = nradmax + 1;
  const int w_l = idx_sph_max + 1;      // includes the trash row
  const int w_n = nradmax + 1;

  MemKK::realloc_kokkos(d_A_off, "pace:A_off", nelements, idx_ms_combs_max, rankmax);
  MemKK::realloc_kokkos(d_w_off, "pace:w_off", nelements, idx_ms_combs_max, rankmax);
  MemKK::realloc_kokkos(d_wm_off, "pace:wm_off", nelements, idx_ms_combs_max, rankmax);
  MemKK::realloc_kokkos(d_r1_off, "pace:r1_off", nelements, idx_ms_combs_max);
  MemKK::realloc_kokkos(d_dB_off, "pace:dB_off", nelements, idx_ms_combs_max + 1);

  auto h_A_off = Kokkos::create_mirror_view(d_A_off);
  auto h_w_off = Kokkos::create_mirror_view(d_w_off);
  auto h_wm_off = Kokkos::create_mirror_view(d_wm_off);
  auto h_r1_off = Kokkos::create_mirror_view(d_r1_off);
  auto h_dB_off = Kokkos::create_mirror_view(d_dB_off);

  auto h_idx_funcs = Kokkos::create_mirror_view(d_idx_funcs);
  auto h_rank = Kokkos::create_mirror_view(d_rank);
  auto h_mus = Kokkos::create_mirror_view(d_mus);
  auto h_ns = Kokkos::create_mirror_view(d_ns);
  auto h_ls = Kokkos::create_mirror_view(d_ls);
  auto h_ms = Kokkos::create_mirror_view(d_ms_combs);
  auto h_isph = Kokkos::create_mirror_view(d_idx_sph_cpu);
  auto h_count = Kokkos::create_mirror_view(d_idx_ms_combs_count);
  auto h_nms1 = Kokkos::create_mirror_view(d_nms_rank1);
  Kokkos::deep_copy(h_idx_funcs, d_idx_funcs);
  Kokkos::deep_copy(h_rank, d_rank);
  Kokkos::deep_copy(h_mus, d_mus);
  Kokkos::deep_copy(h_ns, d_ns);
  Kokkos::deep_copy(h_ls, d_ls);
  Kokkos::deep_copy(h_ms, d_ms_combs);
  Kokkos::deep_copy(h_isph, d_idx_sph_cpu);
  Kokkos::deep_copy(h_count, d_idx_ms_combs_count);
  Kokkos::deep_copy(h_nms1, d_nms_rank1);

  dB_total_max = 0;
  for (int mu_i = 0; mu_i < nelements; mu_i++) {
    const int nms = h_count(mu_i);
    const int nms1 = h_nms1(mu_i);

    // rank-compacted prefix offsets for the dB store
    int run = 0;
    for (int idx = 0; idx < nms; idx++) {
      h_dB_off(mu_i, idx) = run;
      if (idx >= nms1) run += h_rank(mu_i, h_idx_funcs(mu_i, idx));
    }
    h_dB_off(mu_i, nms) = run;
    dB_total_max = MAX(dB_total_max, run);

    for (int idx = 0; idx < nms1; idx++) {
      const int idx_func = h_idx_funcs(mu_i, idx);
      h_r1_off(mu_i, idx) = h_mus(mu_i, idx_func, 0) * nradbase
                            + h_ns(mu_i, idx_func, 0) - 1;
    }

    for (int idx = nms1; idx < nms; idx++) {
      const int idx_func = h_idx_funcs(mu_i, idx);
      const int rank = h_rank(mu_i, idx_func);
      for (int t = 0; t < rank; t++) {
        const int mu = h_mus(mu_i, idx_func, t);
        const int n = h_ns(mu_i, idx_func, t);
        const int l = h_ls(mu_i, idx_func, t);
        const int m = h_ms(mu_i, idx, t);
        const int lbase = l * (l + 1);
        h_A_off(mu_i, idx, t) = (mu * A_l + lbase + m) * A_n + n - 1;
        h_w_off(mu_i, idx, t) = (mu * w_l + h_isph(lbase + m)) * w_n + n - 1;
        // low bit carries the (-1)^m factor of the half-basis
        h_wm_off(mu_i, idx, t) =
            (((mu * w_l + h_isph(lbase - m)) * w_n + n - 1) << 1) | (m & 1);
      }
    }
  }

  Kokkos::deep_copy(d_A_off, h_A_off);
  Kokkos::deep_copy(d_w_off, h_w_off);
  Kokkos::deep_copy(d_wm_off, h_wm_off);
  Kokkos::deep_copy(d_r1_off, h_r1_off);
  Kokkos::deep_copy(d_dB_off, h_dB_off);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

template<class DeviceType>
double PairPACEKokkos<DeviceType>::init_one(int i, int j)
{
  double cutone = PairPACE::init_one(i,j);

  k_scale.view_host()(i,j) = k_scale.view_host()(j,i) = scale[i][j];
  k_scale.modify_host();

  k_cutsq.view_host()(i,j) = k_cutsq.view_host()(j,i) = cutone*cutone;
  k_cutsq.modify_host();

  return cutone;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

template<class DeviceType>
void PairPACEKokkos<DeviceType>::coeff(int narg, char **arg)
{
  PairPACE::coeff(narg,arg);

  // Set up element lists

  auto h_map = Kokkos::create_mirror_view(d_map);

  for (int i = 1; i <= atom->ntypes; i++)
    h_map(i) = map[i];

  Kokkos::deep_copy(d_map,h_map);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairPACEKokkos<DeviceType>::allocate()
{
  PairPACE::allocate();

  int n = atom->ntypes + 1;
  MemKK::realloc_kokkos(d_map, "pace:map", n);

  MemKK::realloc_kokkos(k_cutsq, "pace:cutsq", n, n);
  d_cutsq = k_cutsq.template view<DeviceType>();

  MemKK::realloc_kokkos(k_scale, "pace:scale", n, n);
  d_scale = k_scale.template view<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
struct FindMaxNumNeighs {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  NeighListKokkos<DeviceType> k_list;

  FindMaxNumNeighs(NeighListKokkos<DeviceType>* nl): k_list(*nl) {}
  ~FindMaxNumNeighs() {k_list.copymode = 1;}

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator() (const int& ii, int& maxneigh) const {
    const int i = k_list.d_ilist[ii];
    const int num_neighs = k_list.d_numneigh[i];
    if (maxneigh < num_neighs) maxneigh = num_neighs;
  }
};

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairPACEKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
  if (host_fallback) {
    atomKK->sync(Host,X_MASK|TYPE_MASK);
    PairPACE::compute(eflag_in,vflag_in);
    atomKK->modified(Host,F_MASK);
    return;
  }

  eflag = eflag_in;
  vflag = vflag_in;

  if (neighflag == FULL) no_virial_fdotr_compute = 1;

  ev_init(eflag,vflag,0);

  // reallocate per-atom arrays if necessary

  if (eflag_atom) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->create_kokkos(k_eatom,eatom,maxeatom,"pair:eatom");
    d_eatom = k_eatom.view<DeviceType>();
  }
  if (vflag_atom) {
    memoryKK->destroy_kokkos(k_vatom,vatom);
    memoryKK->create_kokkos(k_vatom,vatom,maxvatom,"pair:vatom");
    d_vatom = k_vatom.view<DeviceType>();
  }
  if (flag_corerep_factor && atom->nlocal > nmax_corerep) {
    memory->destroy(corerep_factor);
    nmax_corerep = atom->nlocal;
    memory->create(corerep_factor, nmax_corerep, "pace/atom:corerep");
    //zeroify array
    memset(corerep_factor, 0, nmax_corerep * sizeof(*corerep_factor));
  }

  copymode = 1;
  if (!force->newton_pair)
    error->all(FLERR,"PairPACEKokkos requires 'newton on'");

  if (recursive)
    error->all(FLERR,"Must use 'product' algorithm with pair pace/kk");

  atomKK->sync(execution_space,X_MASK|F_MASK|TYPE_MASK);
  x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  k_scale.template sync<DeviceType>();
  k_cutsq.template sync<DeviceType>();

  NeighListKokkos<DeviceType>* k_list = static_cast<NeighListKokkos<DeviceType>*>(list);
  d_numneigh = k_list->d_numneigh;
  d_neighbors = k_list->d_neighbors;
  d_ilist = k_list->d_ilist;
  inum = list->inum;

  need_dup = lmp->kokkos->need_dup<DeviceType>();
  if (need_dup) {
    dup_f     = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(f);
    dup_vatom = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(d_vatom);
  } else {
    ndup_f     = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(f);
    ndup_vatom = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(d_vatom);
  }

  maxneigh = 0;
  Kokkos::parallel_reduce("pace::find_maxneigh", inum, FindMaxNumNeighs<DeviceType>(k_list), Kokkos::Max<int>(maxneigh));

  int vector_length_default = 1;

  chunk_size = MIN(chunksize,inum); // "chunksize" variable is set by user
  chunk_offset = 0;

  grow(chunk_size, maxneigh);

  EV_FLOAT ev;

  while (chunk_offset < inum) { // chunk up loop to prevent running out of memory

    // weights_re/weights_im/weights_rank1 are zeroed via first-touch in
    // ComputeFS (one thread per atom, run just before ComputeWeights).
    if constexpr (host_flag) {
      Kokkos::deep_copy(A_sph, complex::zero());
      // bulk zeroing beats the per-atom first-touch loop on a CPU backend
      Kokkos::deep_copy(weights, complex::zero());
      Kokkos::deep_copy(weights_rank1, 0.0);
    } else {
      Kokkos::deep_copy(A_sph_re, 0.0);
      Kokkos::deep_copy(A_sph_im, 0.0);
    }
    Kokkos::deep_copy(A_rank1, 0.0);
    Kokkos::deep_copy(rhos, 0.0);
    Kokkos::deep_copy(rho_core, 0.0);
    Kokkos::deep_copy(d_d_min, PairPACE::aceimpl->basis_set->cutoffmax);
    Kokkos::deep_copy(d_jj_min, -1);
    Kokkos::deep_copy(d_corerep, 0.0);

    EV_FLOAT ev_tmp;

    if (chunk_size > inum - chunk_offset)
      chunk_size = inum - chunk_offset;

    //Neigh
    {
      int vector_length = vector_length_default;
      int team_size = team_size_compute_neigh;
      check_team_size_for<TagPairPACEComputeNeigh>(chunk_size,team_size,vector_length);
      int scratch_size = scratch_size_helper<int>(team_size * maxneigh);
      typename Kokkos::TeamPolicy<DeviceType, TagPairPACEComputeNeigh> policy_neigh(chunk_size,team_size,vector_length);
      policy_neigh = policy_neigh.set_scratch_size(0, Kokkos::PerTeam(scratch_size));
      Kokkos::parallel_for("ComputeNeigh",policy_neigh,*this);
    }

    // ComputeRadial runs as its own kernel only on the host; on the device the
    // radial evaluation is fused into ComputeAi (one less launch, no round trip)
    if constexpr (host_flag)
      Kokkos::parallel_for("ComputeRadial", host_atom_policy<TagPairPACEComputeRadialCPU>(), *this);

    //ComputeAi (radial evaluation fused in on the device)
    if constexpr (host_flag) {
      // one atom per thread, neighbors looped inside: no atomics, and the
      // atom's A_sph block stays in cache for the whole neighbor loop
      Kokkos::parallel_for("ComputeAi", host_atom_policy<TagPairPACEComputeAiCPU>(), *this);
    } else {
      int vector_length = vector_length_default;
      int team_size = team_size_compute_ai;
      check_team_size_for<TagPairPACEComputeAi>(((chunk_size+team_size-1)/team_size)*maxneigh,team_size,vector_length);
      typename Kokkos::TeamPolicy<DeviceType, TagPairPACEComputeAi> policy_ai(((chunk_size+team_size-1)/team_size)*maxneigh,team_size,vector_length);
      Kokkos::parallel_for("ComputeAi",policy_ai,*this);
    }

    // ConjugateAi runs on the host only: it expands the packed half-basis into
    // the full A array the CPU kernels gather from. The device path reads
    // A_sph directly through read_A() and needs no expansion pass.
    if constexpr (host_flag) {
      typename Kokkos::RangePolicy<DeviceType,TagPairPACEConjugateAi> policy_conj_ai(0,chunk_size);
      Kokkos::parallel_for("ConjugateAi",policy_conj_ai,*this);
    }

    //ComputeRho
    if constexpr (host_flag) {
      const int nbatch = (chunk_size + PACE_VLEN - 1) / PACE_VLEN;
      Kokkos::parallel_for("ComputeRhoFSWeights",
          Kokkos::RangePolicy<DeviceType, Kokkos::Schedule<Kokkos::Dynamic>,
                              TagPairPACEComputeRhoBatchCPU>(0, nbatch), *this);
    } else {
      typename Kokkos::RangePolicy<DeviceType,TagPairPACEComputeRho> policy_rho(0,chunk_size*idx_ms_combs_max);
      Kokkos::parallel_for("ComputeRho",policy_rho,*this);
    }

    //ComputeFS
    if constexpr (!host_flag) {
      typename Kokkos::RangePolicy<DeviceType,TagPairPACEComputeFS> policy_fs(0,chunk_size);
      Kokkos::parallel_for("ComputeFS",policy_fs,*this);
    }

    //ComputeWeights (fused into ComputeRho on the host)
    if constexpr (!host_flag) {
      typename Kokkos::RangePolicy<DeviceType,TagPairPACEComputeWeights> policy_weights(0,chunk_size * idx_ms_combs_max);
      Kokkos::parallel_for("ComputeWeights",policy_weights,*this);
    }

    //ComputeDerivative
    if constexpr (host_flag) {
      Kokkos::parallel_for("ComputeDerivative", host_atom_policy<TagPairPACEComputeDerivativeCPU>(), *this);
    } else {
      int vector_length = vector_length_default;
      int team_size = team_size_compute_derivative;
      check_team_size_for<TagPairPACEComputeDerivative>(((chunk_size+team_size-1)/team_size)*maxneigh,team_size,vector_length);
      typename Kokkos::TeamPolicy<DeviceType,TagPairPACEComputeDerivative> policy_derivative(((chunk_size+team_size-1)/team_size)*maxneigh,team_size,vector_length);
      Kokkos::parallel_for("ComputeDerivative",policy_derivative,*this);
    }

    //ComputeForce
    {
      if (evflag) {
        if (neighflag == HALF) {
          typename Kokkos::RangePolicy<DeviceType,TagPairPACEComputeForce<HALF,1> > policy_force(0,chunk_size);
          Kokkos::parallel_reduce(policy_force, *this, ev_tmp);
        } else if (neighflag == HALFTHREAD) {
          typename Kokkos::RangePolicy<DeviceType,TagPairPACEComputeForce<HALFTHREAD,1> > policy_force(0,chunk_size);
          Kokkos::parallel_reduce("ComputeForce",policy_force, *this, ev_tmp);
        }
      } else {
        if (neighflag == HALF) {
          typename Kokkos::RangePolicy<DeviceType,TagPairPACEComputeForce<HALF,0> > policy_force(0,chunk_size);
          Kokkos::parallel_for(policy_force, *this);
        } else if (neighflag == HALFTHREAD) {
          typename Kokkos::RangePolicy<DeviceType,TagPairPACEComputeForce<HALFTHREAD,0> > policy_force(0,chunk_size);
          Kokkos::parallel_for("ComputeForce",policy_force, *this);
        }
      }
    }
    ev += ev_tmp;

    if (flag_corerep_factor) {
      h_corerep = Kokkos::create_mirror_view(d_corerep);
      Kokkos::deep_copy(h_corerep,d_corerep);
      memcpy(corerep_factor+chunk_offset, (void *) h_corerep.data(), sizeof(double)*chunk_size);
    }

    chunk_offset += chunk_size;
  } // end while

  if (need_dup)
    Kokkos::Experimental::contribute(f, dup_f);

  if (eflag_global) eng_vdwl += ev.evdwl;
  if (vflag_global) {
    virial[0] += ev.v[0];
    virial[1] += ev.v[1];
    virial[2] += ev.v[2];
    virial[3] += ev.v[3];
    virial[4] += ev.v[4];
    virial[5] += ev.v[5];
  }

  if (vflag_fdotr) pair_virial_fdotr_compute(this);

  if (eflag_atom) {
    k_eatom.template modify<DeviceType>();
    k_eatom.sync_host();
  }

  if (vflag_atom) {
    if (need_dup)
      Kokkos::Experimental::contribute(d_vatom, dup_vatom);
    k_vatom.template modify<DeviceType>();
    k_vatom.sync_host();
  }

  atomKK->modified(execution_space,F_MASK);

  copymode = 0;

  // free duplicated memory
  if (need_dup) {
    dup_f     = {};
    dup_vatom = {};
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairPACEKokkos<DeviceType>::operator() (TagPairPACEComputeNeigh,const typename Kokkos::TeamPolicy<DeviceType, TagPairPACEComputeNeigh>::member_type& team) const
{
  const int ii = team.league_rank();
  const int i = d_ilist[ii + chunk_offset];
  const int itype = type[i];
  const KK_FLOAT xtmp = x(i,0);
  const KK_FLOAT ytmp = x(i,1);
  const KK_FLOAT ztmp = x(i,2);
  const int jnum = d_numneigh[i];
  const int mu_i = d_map(type(i));

  // get a pointer to scratch memory
  // This is used to cache whether or not an atom is within the cutoff
  // If it is, inside is assigned to 1, otherwise -1
  const int team_rank = team.team_rank();
  const int scratch_shift = team_rank * maxneigh; // offset into pointer for entire team
  int* inside = (int*)team.team_shmem().get_shmem(team.team_size() * maxneigh * sizeof(int), 0) + scratch_shift;

  // loop over list of all neighbors within force cutoff
  // distsq[] = distance sq to each
  // rlist[] = distance vector to each
  // nearest[] = atom indices of neighbors

  int ncount = 0;
  Kokkos::parallel_reduce(Kokkos::TeamThreadRange(team,jnum),
      [&] (const int jj, int& count) {
    int j = d_neighbors(i,jj);
    j &= NEIGHMASK;

    const int jtype = type(j);

    const KK_FLOAT delx = xtmp - x(j,0);
    const KK_FLOAT dely = ytmp - x(j,1);
    const KK_FLOAT delz = ztmp - x(j,2);
    const KK_FLOAT rsq = delx*delx + dely*dely + delz*delz;

    inside[jj] = -1;
    if (rsq < d_cutsq(itype,jtype)) {
     inside[jj] = 1;
     count++;
    }
  },ncount);

  d_ncount(ii) = ncount;

  Kokkos::parallel_scan(Kokkos::TeamThreadRange(team,jnum),
      [&] (const int jj, int& offset, bool final) {

    if (inside[jj] < 0) return;

    if (final) {
      int j = d_neighbors(i,jj);
      j &= NEIGHMASK;
      const KK_FLOAT delx = xtmp - x(j,0);
      const KK_FLOAT dely = ytmp - x(j,1);
      const KK_FLOAT delz = ztmp - x(j,2);
      const KK_FLOAT rsq = delx*delx + dely*dely + delz*delz;
      const KK_FLOAT r = sqrt(rsq);
      const KK_FLOAT rinv = 1.0/r;
      const int mu_j = d_map(type(j));
      d_mu(ii,offset) = mu_j;
      d_rnorms(ii,offset) = r;
      d_rhats(ii,offset,0) = -delx*rinv;
      d_rhats(ii,offset,1) = -dely*rinv;
      d_rhats(ii,offset,2) = -delz*rinv;
      d_nearest(ii,offset) = j;
    }
    offset++;
  });

  if (is_zbl) {
    //adapted from https://www.osti.gov/servlets/purl/1429450
    if (ncount > 0) {
      using minloc_value_type=Kokkos::MinLoc<KK_FLOAT,int>::value_type;
      minloc_value_type djjmin;
      djjmin.val=1e20;
      djjmin.loc=-1;
      Kokkos::MinLoc<KK_FLOAT,int> reducer_scalar(djjmin);
      // loop over ncount (actual neighbours withing cutoff) rather than jnum (total number of neigh in cutoff+skin)
      Kokkos::parallel_reduce(Kokkos::TeamThreadRange(team, ncount),
               [&](const int offset, minloc_value_type &min_d_dist) {
                 int j = d_nearest(ii,offset);
                 j &= NEIGHMASK;
                 auto r = d_rnorms(ii,offset);
                 const int mu_j = d_map(type(j));
                 const KK_FLOAT d = r - (d_cut_in(mu_i, mu_j) - d_dcut_in(mu_i, mu_j));
                 if (d < min_d_dist.val) {
                   min_d_dist.val = d;
                   min_d_dist.loc = offset;
                 }
       }, reducer_scalar);
      d_d_min(ii) = djjmin.val;
      d_jj_min(ii) = djjmin.loc;// d_jj_min should be NOT in 0..jnum range, but in 0..d_ncount(<=jnum)
    } else {
      d_d_min(ii) = 1e20;
      d_jj_min(ii) = -1;
    }
  }
}

/* ----------------------------------------------------------------------
   CPU backend: one atom per thread.  Besides dropping the two integer
   divisions per work item, this skips the maxneigh - ncount empty items
   the flattened team policy has to launch for every atom.
------------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairPACEKokkos<DeviceType>::operator() (TagPairPACEComputeRadialCPU, const int& ii) const
{
  // Host-only kernel. The body is guarded so the explicit instantiation of
  // PairPACEKokkos<LMPDeviceType> does not compile it for the device, where
  // it would pull in host-only helpers and the LayoutRight assumption.
  // "if constexpr (host_flag)" alone is not enough: PairPACEKokkos<LMPHostType>
  // (pace/kk/host) is also explicitly instantiated in GPU builds, and there
  // host_flag is always true, so the device compiler pass still type-checks
  // this body for that instantiation -- hence the extra LMP_KK_DEVICE_COMPILE
  // preprocessor guard, which is false for both instantiations in that pass.
#ifndef LMP_KK_DEVICE_COMPILE
  if constexpr (host_flag) {
    const int i = d_ilist[ii + chunk_offset];
    const int mu_i = d_map(type(i));
    const int ncount = d_ncount(ii);

    // The lookup tables are far larger than L2 and the bin is data dependent,
    // so this kernel is bound by L3/DRAM gather latency (measured: ~52 lines
    // per pair).  Both tables share the same binning, so the next neighbor's
    // rows can be prefetched while the current one is evaluated.
  #if defined(__GNUC__)
    const int mu0 = (ncount > 0) ? (int) d_mu(ii, 0) : 0;
    if (ncount > 0) {
      auto &sp_gk = k_splines_gk.template view<DeviceType>()(mu_i, mu0);
      auto &sp_rnl = k_splines_rnl.template view<DeviceType>()(mu_i, mu0);
      const int nl0 = (int) (d_rnorms(ii, 0) * sp_gk.rscalelookup);
      const KK_FLOAT *g0 = &sp_gk.lookupTable(nl0, 0, 0);
      const KK_FLOAT *r0 = &sp_rnl.lookupTable(nl0, 0, 0);
      for (int q = 0; q < sp_gk.num_of_functions * 4; q += 8) __builtin_prefetch(g0 + q, 0, 1);
      for (int q = 0; q < sp_rnl.num_of_functions * 4; q += 8) __builtin_prefetch(r0 + q, 0, 1);
    }
  #endif

    for (int jj = 0; jj < ncount; jj++) {
  #if defined(__GNUC__)
      if (jj + 1 < ncount) {
        const int mu_n = (int) d_mu(ii, jj + 1);
        auto &sp_gk = k_splines_gk.template view<DeviceType>()(mu_i, mu_n);
        auto &sp_rnl = k_splines_rnl.template view<DeviceType>()(mu_i, mu_n);
        const int nl_n = (int) (d_rnorms(ii, jj + 1) * sp_gk.rscalelookup);
        const KK_FLOAT *gn = &sp_gk.lookupTable(nl_n, 0, 0);
        const KK_FLOAT *rn = &sp_rnl.lookupTable(nl_n, 0, 0);
        for (int q = 0; q < sp_gk.num_of_functions * 4; q += 8) __builtin_prefetch(gn + q, 0, 1);
        for (int q = 0; q < sp_rnl.num_of_functions * 4; q += 8) __builtin_prefetch(rn + q, 0, 1);
      }
  #endif
      evaluate_splines(ii, jj, d_rnorms(ii, jj), nradbase, nradmax, mu_i, d_mu(ii, jj));
    }
  }
#endif
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairPACEKokkos<DeviceType>::operator() (TagPairPACEComputeAi, const typename Kokkos::TeamPolicy<DeviceType, TagPairPACEComputeAi>::member_type& team) const
{
  // Extract the atom number
  int ii = team.team_rank() + team.team_size() * (team.league_rank() %
           ((chunk_size+team.team_size()-1)/team.team_size()));
  if (ii >= chunk_size) return;

  // Extract the neighbor number
  const int jj = team.league_rank() / ((chunk_size+team.team_size()-1)/team.team_size());
  const int ncount = d_ncount(ii);
  if (jj >= ncount) return;

  // several teams accumulate into the same atom, so the adds must be atomic
  compute_ai_one<true>(ii, jj);
}

/* ----------------------------------------------------------------------
   CPU backend: one thread owns atom ii and walks its whole neighbor list,
   so nothing else writes A_sph/A_rank1/rho_core for this atom and the
   accumulation needs no atomics
------------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairPACEKokkos<DeviceType>::operator() (TagPairPACEComputeAiCPU, const int& ii) const
{
  // Host-only kernel. The body is guarded so the explicit instantiation of
  // PairPACEKokkos<LMPDeviceType> does not compile it for the device, where
  // it would pull in host-only helpers and the LayoutRight assumption.
  // "if constexpr (host_flag)" alone is not enough: PairPACEKokkos<LMPHostType>
  // (pace/kk/host) is also explicitly instantiated in GPU builds, and there
  // host_flag is always true, so the device compiler pass still type-checks
  // this body for that instantiation -- hence the extra LMP_KK_DEVICE_COMPILE
  // preprocessor guard, which is false for both instantiations in that pass.
#ifndef LMP_KK_DEVICE_COMPILE
  if constexpr (host_flag) {
    // measured: batching this kernel across neighbors loses -- its arithmetic
    // is too small to amortize the staging and the cross-lane reduction of the
    // A_sph accumulators, unlike ComputeDerivative where each lane owns its
    // own output
    const int ncount = d_ncount(ii);
    for (int jj = 0; jj < ncount; jj++)
      compute_ai_one<false>(ii, jj);
  }
#endif
}

/* ---------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Host only: expand the packed half-basis A_sph into the full (l,m) A array
   used by the CPU basis-function kernels, filling the -m entries by the
   conjugate symmetry A(l,-m) = (-1)^m conj(A(l,m)). The device kernels skip
   this pass entirely and apply the same symmetry on the fly in read_A().
------------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairPACEKokkos<DeviceType>::operator() (TagPairPACEConjugateAi, const int& ii) const
{
  // Host-only kernel. The body is guarded so the explicit instantiation of
  // PairPACEKokkos<LMPDeviceType> does not compile it for the device, where
  // it would pull in host-only helpers and the LayoutRight assumption.
  if constexpr (host_flag) {
    for (int mu_j = 0; mu_j < nelements; mu_j++) {

      // transpose
      int idx_sph = 0;
      for (int m = 0; m <= lmax; m++) {
        for (int l = m; l <= lmax; l++) {
          const int idx = l * (l + 1) + m;
          for (int n = 0; n < nradmax; n++)
            A(ii, mu_j, idx, n) = A_sph(ii, mu_j, idx_sph, n);
          idx_sph++;
        }
      }

      // complex conjugate A's (for NEGATIVE (-m) terms) for rank > 1
      for (int l = 0; l <= lmax; l++) {
        for (int m = 1; m <= l; m++) {
          const int idx = l * (l + 1) + m;   // (l, m)
          const int idxm = l * (l + 1) - m;  // (l, -m)
          const int idx_sph = d_idx_sph(idx);
          const int factor = m % 2 == 0 ? 1 : -1;
          for (int n = 0; n < nradmax; n++)
            A(ii, mu_j, idxm, n) = A_sph(ii, mu_j, idx_sph, n).conj() * (KK_FLOAT)factor;
        }
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<bool NEED_ATOMICS>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairPACEKokkos<DeviceType>::compute_ai_one(const int ii, const int jj) const
{
  const int mu_j = d_mu(ii, jj);

  if constexpr (NEED_ATOMICS) {
    // Device: the radial evaluation is fused in here rather than run as its
    // own ComputeRadial kernel -- it saves a launch and the round trip of the
    // radial arrays, and each thread owns its (ii,jj) slot so there is no
    // cross-thread dependency. The host path keeps ComputeRadialCPU, which
    // walks a whole atom and can prefetch the next neighbor's spline rows.
    const int i = d_ilist[ii + chunk_offset];
    const KK_FLOAT r_norm = d_rnorms(ii, jj);
    const int mu_i = d_map(type(i));
    evaluate_splines(ii, jj, r_norm, nradbase, nradmax, mu_i, mu_j);
  }

  // Hoist the (atom, neighbor) and (atom, element) parts of the subscripts out
  // of the innermost loops.  Using View::stride() keeps this correct for any
  // layout, so it is a strength reduction rather than a layout assumption.
  const KK_FLOAT * const fr_ij = &fr(ii, jj, 0);
  const KK_FLOAT * const gr_ij = &gr(ii, jj, 0);
  KK_FLOAT * const A_rank1_i = &A_rank1(ii, mu_j, 0);
  const int fr_sl = 1;              // fr is flat (n*(lmax+1)+l): l-stride
  const int fr_sn = lmax + 1;       // ... and n-stride
  const int gr_sn = (int) gr.stride(2);
  const int ar1_sn = (int) A_rank1.stride(2);

  // A_sph is the host-only interleaved array; the device accumulates into
  // A_sph_re/A_sph_im instead and grow() never allocates A_sph there, so its
  // address must not be taken on that path. With a null base and zero stride
  // the derived a_sph below stays a well-defined (never dereferenced) null.
  complex *A_sph_i = nullptr;
  int sph_ss = 0, sph_sn = 0;
  if constexpr (!NEED_ATOMICS) {
    A_sph_i = &A_sph(ii, mu_j, 0, 0);
    sph_ss = (int) A_sph.stride(2);
    sph_sn = (int) A_sph.stride(3);
  }

  // rank = 1
  for (int n = 0; n < nradbase; n++) {
    if constexpr (NEED_ATOMICS)
      Kokkos::atomic_add(&A_rank1_i[n * ar1_sn], gr_ij[n * gr_sn] * Y00);
    else
      A_rank1_i[n * ar1_sn] += gr_ij[n * gr_sn] * Y00;
  }

  // rank > 1

  // Compute plm and ylm

  // requires rx^2 + ry^2 + rz^2 = 1 , NO CHECKING IS PERFORMED !!!!!!!!!
  // requires -1 <= rz <= 1 , NO CHECKING IS PERFORMED !!!!!!!!!
  // prefactors include 1/sqrt(2) factor compared to reference

  complex ylm, phase;
  complex phasem, mphasem1;
  complex dyx, dyy, dyz;
  complex rdy;

  const KK_FLOAT rx = d_rhats(ii, jj, 0);
  const KK_FLOAT ry = d_rhats(ii, jj, 1);
  const KK_FLOAT rz = d_rhats(ii, jj, 2);

  phase.re = rx;
  phase.im = ry;

  KK_FLOAT plm_idx,plm_idx1,plm_idx2;

  plm_idx = plm_idx1 = plm_idx2 = 0.0;

  int idx_sph = 0;

  // m = 0
  for (int l = 0; l <= lmax; l++) {
    // const int idx = l * (l + 1);

    if (l == 0) {
      // l=0, m=0
      // plm[0] = Y00/sq1o4pi; //= sq1o4pi;
      plm_idx = Y00; //= 1;
    } else if (l == 1) {
      // l=1, m=0
      plm_idx = Y00 * sq3 * rz;
    } else {
      // l>=2, m=0
      plm_idx = alm(idx_sph) * (rz * plm_idx1 + blm(idx_sph) * plm_idx2);
    }

    ylm.re = plm_idx;
    ylm.im = 0.0;

    complex * const a_sph = A_sph_i + idx_sph * sph_ss;
    const KK_FLOAT * const r_nl = fr_ij + l * fr_sl;
    for (int n = 0; n < nradmax; n++) {
      if constexpr (NEED_ATOMICS) {
        // device: split re/im arrays, several teams share the atom
        Kokkos::atomic_add(&A_sph_re(ii, mu_j, idx_sph, n), fr(ii, jj, n * (lmax + 1) + l) * ylm.re);
        Kokkos::atomic_add(&A_sph_im(ii, mu_j, idx_sph, n), fr(ii, jj, n * (lmax + 1) + l) * ylm.im);
      } else {
        // host: interleaved A_sph through hoisted pointers, no atomics
        a_sph[n * sph_sn].re += r_nl[n * fr_sn] * ylm.re;
        a_sph[n * sph_sn].im += r_nl[n * fr_sn] * ylm.im;
      }
    }

    plm_idx2 = plm_idx1;
    plm_idx1 = plm_idx;

    idx_sph++;
  }

  plm_idx = plm_idx1 = plm_idx2 = 0.0;

  // m = 1
  for (int l = 1; l <= lmax; l++) {
    // const int idx = l * (l + 1) + 1; // (l, 1)

    if (l == 1) {
      // l=1, m=1
      plm_idx = -sq3o2 * Y00;
    } else if (l == 2) {
      const KK_FLOAT t = dl(l) * plm_idx1;
      plm_idx = t * rz;
    } else {
      plm_idx = alm(idx_sph) * (rz * plm_idx1 + blm(idx_sph) * plm_idx2);
    }

    ylm = phase * plm_idx;

    complex * const a_sph = A_sph_i + idx_sph * sph_ss;
    const KK_FLOAT * const r_nl = fr_ij + l * fr_sl;
    for (int n = 0; n < nradmax; n++) {
      if constexpr (NEED_ATOMICS) {
        // device: split re/im arrays, several teams share the atom
        Kokkos::atomic_add(&A_sph_re(ii, mu_j, idx_sph, n), fr(ii, jj, n * (lmax + 1) + l) * ylm.re);
        Kokkos::atomic_add(&A_sph_im(ii, mu_j, idx_sph, n), fr(ii, jj, n * (lmax + 1) + l) * ylm.im);
      } else {
        // host: interleaved A_sph through hoisted pointers, no atomics
        a_sph[n * sph_sn].re += r_nl[n * fr_sn] * ylm.re;
        a_sph[n * sph_sn].im += r_nl[n * fr_sn] * ylm.im;
      }
    }

    plm_idx2 = plm_idx1;
    plm_idx1 = plm_idx;

    idx_sph++;
  }

  plm_idx = plm_idx1 = plm_idx2 = 0.0;

  KK_FLOAT plm_mm1_mm1 = -sq3o2 * Y00; // (1, 1)

  // m > 1
  phasem = phase;
  for (int m = 2; m <= lmax; m++) {

    mphasem1.re = phasem.re * KK_FLOAT(m);
    mphasem1.im = phasem.im * KK_FLOAT(m);
    phasem = phasem * phase;

    for (int l = m; l <= lmax; l++) {
      // const int idx = l * (l + 1) + m;

      if (l == m) {
        plm_idx = cl(l) * plm_mm1_mm1; // (m+1, m)
        plm_mm1_mm1 = plm_idx;
      } else if (l == (m + 1)) {
        const KK_FLOAT t = dl(l) * plm_mm1_mm1; // (m - 1, m - 1)
        plm_idx = t * rz; // (m, m)
      } else {
        plm_idx = alm(idx_sph) * (rz * plm_idx1 + blm(idx_sph) * plm_idx2);
      }

      ylm.re = phasem.re * plm_idx;
      ylm.im = phasem.im * plm_idx;

      complex * const a_sph = A_sph_i + idx_sph * sph_ss;
      const KK_FLOAT * const r_nl = fr_ij + l * fr_sl;
      for (int n = 0; n < nradmax; n++) {
        if constexpr (NEED_ATOMICS) {
          // device: split re/im arrays, several teams share the atom
          Kokkos::atomic_add(&A_sph_re(ii, mu_j, idx_sph, n), fr(ii, jj, n * (lmax + 1) + l) * ylm.re);
          Kokkos::atomic_add(&A_sph_im(ii, mu_j, idx_sph, n), fr(ii, jj, n * (lmax + 1) + l) * ylm.im);
        } else {
          // host: interleaved A_sph through hoisted pointers, no atomics
          a_sph[n * sph_sn].re += r_nl[n * fr_sn] * ylm.re;
          a_sph[n * sph_sn].im += r_nl[n * fr_sn] * ylm.im;
        }
      }

      plm_idx2 = plm_idx1;
      plm_idx1 = plm_idx;

      idx_sph++;
    }
  }

  // hard-core repulsion
  if constexpr (NEED_ATOMICS)
    Kokkos::atomic_add(&rho_core(ii), cr(ii, jj));
  else
    rho_core(ii) += cr(ii, jj);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
typename PairPACEKokkos<DeviceType>::complex
PairPACEKokkos<DeviceType>::read_A(const int ii, const int mu, const int l, const int m, const int n) const
{
  complex A_t;
  if (m >= 0) {
    const int idx_sph = d_idx_sph(l * (l + 1) + m);
    A_t.re = A_sph_re(ii, mu, idx_sph, n);
    A_t.im = A_sph_im(ii, mu, idx_sph, n);
  } else {
    const int p = -m;
    const int idx_sph = d_idx_sph(l * (l + 1) + p);
    const KK_FLOAT factor = (p % 2 == 0) ? 1.0 : -1.0;
    A_t.re =  A_sph_re(ii, mu, idx_sph, n) * factor;
    A_t.im = -A_sph_im(ii, mu, idx_sph, n) * factor;
  }
  return A_t;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairPACEKokkos<DeviceType>::operator() (TagPairPACEComputeRho, const int& iter) const
{
  const int idx_ms_combs = iter / chunk_size;
  const int ii = iter % chunk_size;

  const int i = d_ilist[ii + chunk_offset];
  const int mu_i = d_map(type(i));

  if (idx_ms_combs >= d_idx_ms_combs_count(mu_i)) return;

  const int ndensity = d_ndensity(mu_i);

  const int idx_func = d_idx_funcs(mu_i, idx_ms_combs);
  const int rank = d_rank(mu_i, idx_func);

  // Basis function B with iterative product and density rho(p) calculation
  if (rank == 1) {
    const int mu = d_mus(mu_i, idx_func, 0);
    const int n = d_ns(mu_i, idx_func, 0);
    KK_FLOAT A_cur = A_rank1(ii, mu, n - 1);
    for (int p = 0; p < ndensity; ++p) {
      //for rank=1 (r=0) only 1 ms-combination exists (ms_ind=0), so index of func.ctildes is 0..ndensity-1
      Kokkos::atomic_add(&rhos(ii, p), d_ctildes(mu_i, idx_ms_combs, p) * A_cur);
    }
  } else { // rank > 1
    // B = product of A over the ms-combination, accumulated in a register (no
    // global product-chain scratch). The leave-one-out products needed for the
    // weights are recomputed in ComputeWeights rather than stored in dB_flatten.
    complex B = complex::one();
    for (int t = 0; t < rank; t++) {
      //TODO: optimize ns[t]-1 -> ns[t] during functions construction
      const int mu = d_mus(mu_i, idx_func, t);
      const int n = d_ns(mu_i, idx_func, t);
      const int l = d_ls(mu_i, idx_func, t);
      const int m = d_ms_combs(mu_i, idx_ms_combs, t); // current ms-combination (of length = rank)
      B = B * read_A(ii, mu, l, m, n - 1);
    }

    for (int p = 0; p < ndensity; ++p) {
      // real-part only multiplication
      Kokkos::atomic_add(&rhos(ii, p), B.real_part_product(d_ctildes(mu_i, idx_ms_combs, p)));
    }
  }
}

/* ----------------------------------------------------------------------
   CPU backend: one thread owns atom ii and loops over all of its basis
   functions.  A_list and A_forward_prod never outlive a single
   (atom, ms-combination) iteration, so on the CPU they are stack arrays of
   at most rankmax+1 entries instead of chunk-sized global arrays.  That
   removes the dominant store stream of the whole evaluator and, because
   only this thread touches rhos(ii,:), the density accumulation needs no
   atomics either.
------------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairPACEKokkos<DeviceType>::operator() (TagPairPACEComputeRhoCPU, const int& ii) const
{
  // Host-only kernel. The body is guarded so the explicit instantiation of
  // PairPACEKokkos<LMPDeviceType> does not compile it for the device, where
  // it would pull in host-only helpers and the LayoutRight assumption.
  // "if constexpr (host_flag)" alone is not enough: PairPACEKokkos<LMPHostType>
  // (pace/kk/host) is also explicitly instantiated in GPU builds, and there
  // host_flag is always true, so the device compiler pass still type-checks
  // this body for that instantiation -- hence the extra LMP_KK_DEVICE_COMPILE
  // preprocessor guard, which is false for both instantiations in that pass.
#ifndef LMP_KK_DEVICE_COMPILE
  if constexpr (host_flag) {
    // ndensity is 1 or 2 for every published ACE potential; fixing it at
    // compile time lets the short coefficient loops unroll
    const int nd = d_ndensity(d_map(type(d_ilist[ii + chunk_offset])));
    if (nd == 2) rho_fs_weights_cpu<2>(ii);
    else if (nd == 1) rho_fs_weights_cpu<1>(ii);
    else rho_fs_weights_cpu<0>(ii);
  }
#endif
}

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairPACEKokkos<DeviceType>::operator() (TagPairPACEComputeRhoBatchCPU, const int& ib) const
{
  // Host-only kernel. The body is guarded so the explicit instantiation of
  // PairPACEKokkos<LMPDeviceType> does not compile it for the device, where
  // it would pull in host-only helpers and the LayoutRight assumption.
  // "if constexpr (host_flag)" alone is not enough: PairPACEKokkos<LMPHostType>
  // (pace/kk/host) is also explicitly instantiated in GPU builds, and there
  // host_flag is always true, so the device compiler pass still type-checks
  // this body for that instantiation -- hence the extra LMP_KK_DEVICE_COMPILE
  // preprocessor guard, which is false for both instantiations in that pass.
#ifndef LMP_KK_DEVICE_COMPILE
  if constexpr (host_flag) {
    const int ii0 = ib * PACE_VLEN;
    const int nb = (chunk_size - ii0 < PACE_VLEN) ? (chunk_size - ii0) : PACE_VLEN;

    // the shared basis tables require every lane to have the same element
    const int mu_i = d_map(type(d_ilist[ii0 + chunk_offset]));
    bool uniform = true;
    for (int k = 1; k < nb; k++)
      if ((int) d_map(type(d_ilist[ii0 + k + chunk_offset])) != mu_i) { uniform = false; break; }

    if (!uniform) {
      for (int k = 0; k < nb; k++) {
        const int nd = d_ndensity(d_map(type(d_ilist[ii0 + k + chunk_offset])));
        if (nd == 2) rho_fs_weights_cpu<2>(ii0 + k);
        else if (nd == 1) rho_fs_weights_cpu<1>(ii0 + k);
        else rho_fs_weights_cpu<0>(ii0 + k);
      }
      return;
    }

    const int nd = d_ndensity(mu_i);
    if (nd == 2) rho_fs_weights_batch_cpu<2>(ii0, nb, mu_i);
    else if (nd == 1) rho_fs_weights_batch_cpu<1>(ii0, nb, mu_i);
    else rho_fs_weights_batch_cpu<0>(ii0, nb, mu_i);
  }
#endif
}

/* ----------------------------------------------------------------------
   the fused density -> F(rho) -> weights pass for a batch of atoms of the
   same element.  Every basis-table entry (offsets, coefficients, ranks) is
   read once per batch instead of once per atom; only the per-atom A blocks,
   dB stores and weights blocks are touched per lane.  On this evaluator's
   validated profile the per-atom re-streaming of those shared tables from
   L3 was the largest remaining traffic term.
------------------------------------------------------------------------- */

template<class DeviceType>
template<int NDENSITY>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairPACEKokkos<DeviceType>::rho_fs_weights_batch_cpu(const int ii0, const int nb,
                                                          const int mu_i) const
{
  const int ndensity = NDENSITY ? NDENSITY : d_ndensity(mu_i);
  const int nms = d_idx_ms_combs_count(mu_i);
  const int nms1 = d_nms_rank1(mu_i);

  BasisPtrs bl[PACE_VLEN];
  for (int k = 0; k < nb; k++) set_basis_ptrs(bl[k], ii0 + k, mu_i);
  const BasisPtrs &b0 = bl[0];

  // ---- densities ----
  for (int idx = 0; idx < nms1; idx++) {
    const int ro = b0.r1_off[idx];
    const KK_FLOAT *ct = b0.ctildes + idx * b0.ndensitymax;
    const int ndl = NDENSITY ? NDENSITY : ndensity;
    for (int k = 0; k < nb; k++) {
      const KK_FLOAT A_cur = bl[k].A_rank1[ro];
      for (int p = 0; p < ndl; ++p) bl[k].rho[p] += ct[p] * A_cur;
    }
  }

  for (int idx = nms1; idx < nms; idx++) {
    const int rank = b0.rank[b0.idx_funcs[idx]];
    const int *aoff = b0.A_off + idx * b0.rankmax;
    const int dbo = b0.dB_off[idx];
    const KK_FLOAT *ct = b0.ctildes + idx * b0.ndensitymax;
    const int ndl = NDENSITY ? NDENSITY : ndensity;

    if (rank == 3) {
      const int o0 = aoff[0], o1 = aoff[1], o2 = aoff[2];
      for (int k = 0; k < nb; k++) {
        const complex a0 = bl[k].A[o0], a1 = bl[k].A[o1], a2 = bl[k].A[o2];
        const complex a01 = a0 * a1;
        complex *dB = bl[k].dB + dbo;
        dB[0] = a1 * a2;
        dB[1] = a0 * a2;
        dB[2] = a01;
        const complex B = a01 * a2;
        for (int p = 0; p < ndl; ++p) bl[k].rho[p] += B.real_part_product(ct[p]);
      }
    } else {
      for (int k = 0; k < nb; k++)
        compute_rho_one_cpu<NDENSITY>(bl[k], ndensity, idx);
    }
  }

  // ---- embedding function per lane ----
  for (int k = 0; k < nb; k++) compute_fs_one(ii0 + k);

  // ---- weights ----
  for (int idx = 0; idx < nms1; idx++) {
    const int ro = b0.r1_off[idx];
    const KK_FLOAT *ct = b0.ctildes + idx * b0.ndensitymax;
    const int ndl = NDENSITY ? NDENSITY : ndensity;
    for (int k = 0; k < nb; k++) {
      KK_FLOAT theta = 0.0;
      for (int p = 0; p < ndl; ++p) theta += bl[k].dF[p] * ct[p];
      bl[k].w_rank1[ro] += theta;
    }
  }

  for (int idx = nms1; idx < nms; idx++) {
    const int rank = b0.rank[b0.idx_funcs[idx]];
    const int *woff = b0.w_off + idx * b0.rankmax;
    const int *wmoff = b0.wm_off + idx * b0.rankmax;
    const int dbo = b0.dB_off[idx];
    const KK_FLOAT *ct = b0.ctildes + idx * b0.ndensitymax;
    const int ndl = NDENSITY ? NDENSITY : ndensity;

    if (rank == 3) {
      for (int k = 0; k < nb; k++) {
        KK_FLOAT theta = 0.0;
        for (int p = 0; p < ndl; ++p) theta += bl[k].dF[p] * ct[p];
        theta *= 0.5;
        const complex *dB = bl[k].dB + dbo;
        complex *w = bl[k].w;
        for (int t = 0; t < 3; ++t) {
          const complex value = theta * dB[t];
          w[woff[t]].re += value.re;
          w[woff[t]].im += value.im;
          const int packed = wmoff[t];
          const KK_FLOAT factor = (packed & 1) ? KK_FLOAT(-1.0) : KK_FLOAT(1.0);
          const complex valuem = theta * dB[t].conj() * factor;
          w[packed >> 1].re += valuem.re;
          w[packed >> 1].im += valuem.im;
        }
      }
    } else {
      for (int k = 0; k < nb; k++)
        compute_weights_one_cpu<NDENSITY>(bl[k], ndensity, idx);
    }
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NDENSITY>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairPACEKokkos<DeviceType>::rho_fs_weights_cpu(const int ii) const
{
  const int i = d_ilist[ii + chunk_offset];
  const int mu_i = d_map(type(i));
  const int ndensity = NDENSITY ? NDENSITY : d_ndensity(mu_i);
  const int nms = d_idx_ms_combs_count(mu_i);
  // rank-1 basis functions occupy the first entries, so the two cases are
  // separate loops rather than a branch inside one
  const int nms1 = d_nms_rank1(mu_i);

  // The atom index and the element index are loop invariant here, but the
  // multi-dimensional View subscripts below recompute their contribution to
  // the offset on every one of the millions of accesses this kernel makes.
  // Hoisting them into base pointers and explicit strides removes about half
  // of the kernel's instructions.  This is a host-only code path and the host
  // layout is LayoutRight, so the last index is the contiguous one.
  static_assert(std::is_same_v<typename t_ace_4c::array_layout, Kokkos::LayoutRight>,
                "pace/kk CPU kernels assume LayoutRight on the host");
  BasisPtrs b;
  set_basis_ptrs(b, ii, mu_i);

  // densities
  for (int idx = 0; idx < nms1; idx++)
    rho_one_rank1_cpu<NDENSITY>(b, ndensity, idx);
  for (int idx = nms1; idx < nms; idx++)
    compute_rho_one_cpu<NDENSITY>(b, ndensity, idx);

  // embedding function F(rho) and its derivatives for this atom
  compute_fs_one(ii);

  // weights, reusing the dB products just written for this atom while they
  // are still cache resident
  b.dF = &dF_drho(ii, 0);
  for (int idx = 0; idx < nms1; idx++)
    weights_one_rank1_cpu<NDENSITY>(b, ndensity, idx);
  for (int idx = nms1; idx < nms; idx++)
    compute_weights_one_cpu<NDENSITY>(b, ndensity, idx);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairPACEKokkos<DeviceType>::set_basis_ptrs(BasisPtrs &b, const int ii, const int mu_i) const
{
  b.A       = &A(ii, 0, 0, 0);
  b.A_rank1 = &A_rank1(ii, 0, 0);
  b.dB      = &dB_flatten(ii, 0, 0);
  b.w       = &weights(ii, 0, 0, 0);
  b.w_rank1 = &weights_rank1(ii, 0, 0);
  b.rho     = &rhos(ii, 0);
  b.dF      = &dF_drho(ii, 0);

  b.mus     = &d_mus(mu_i, 0, 0);
  b.ns      = &d_ns(mu_i, 0, 0);
  b.ls      = &d_ls(mu_i, 0, 0);
  b.ms      = &d_ms_combs(mu_i, 0, 0);
  b.ctildes = &d_ctildes(mu_i, 0, 0);
  b.idx_funcs = &d_idx_funcs(mu_i, 0);
  b.idx_sph = &d_idx_sph_cpu(0);
  b.A_off   = &d_A_off(mu_i, 0, 0);
  b.w_off   = &d_w_off(mu_i, 0, 0);
  b.wm_off  = &d_wm_off(mu_i, 0, 0);
  b.r1_off  = &d_r1_off(mu_i, 0);
  b.dB_off  = &d_dB_off(mu_i, 0);
  b.rank    = &d_rank(mu_i, 0);

  b.A_l   = (int) A.extent(2);            // (lmax+1)^2
  b.A_n   = (int) A.extent(3);            // nradmax+1
  b.w_l   = (int) weights.extent(2);      // idx_sph_max
  b.w_n   = (int) weights.extent(3);      // nradmax+1
  b.rankmax = (int) d_ms_combs.extent(2);
  b.ndensitymax = (int) d_ctildes.extent(2);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NDENSITY>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairPACEKokkos<DeviceType>::rho_one_rank1_cpu(const BasisPtrs &b,
                                                   const int ndensity, const int idx) const
{
  const KK_FLOAT A_cur = b.A_rank1[b.r1_off[idx]];
  const KK_FLOAT *ct = b.ctildes + idx * b.ndensitymax;
  const int nd = NDENSITY ? NDENSITY : ndensity;
  for (int p = 0; p < nd; ++p)
    b.rho[p] += ct[p] * A_cur;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NDENSITY>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairPACEKokkos<DeviceType>::weights_one_rank1_cpu(const BasisPtrs &b,
                                                       const int ndensity, const int idx) const
{
  const KK_FLOAT *ct = b.ctildes + idx * b.ndensitymax;
  KK_FLOAT theta = 0.0;
  const int nd = NDENSITY ? NDENSITY : ndensity;
  for (int p = 0; p < nd; ++p)
    theta += b.dF[p] * ct[p];
  b.w_rank1[b.r1_off[idx]] += theta;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NDENSITY>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairPACEKokkos<DeviceType>::compute_rho_one_cpu(const BasisPtrs &b,
                                                     const int ndensity,
                                                     const int idx_ms_combs) const
{
  const int rank = b.rank[b.idx_funcs[idx_ms_combs]];
  const int r = rank - 1;

  const int *aoff = b.A_off + idx_ms_combs * b.rankmax;
  complex *dB     = b.dB    + b.dB_off[idx_ms_combs];

  complex B;
  if (rank == 3) {
    // the dominant correlation order in typical bases; fixed trip count so
    // the product chains unroll fully, and the terms are grouped by rank in
    // the tables so this branch is long-run predictable
    const complex a0 = b.A[aoff[0]], a1 = b.A[aoff[1]], a2 = b.A[aoff[2]];
    const complex a01 = a0 * a1;
    dB[0] = a1 * a2;
    dB[1] = a0 * a2;
    dB[2] = a01;
    B = a01 * a2;
  } else {
    // general rank > 1: forward and backward products over the ms-combination
    complex A_list_l[MAX_RANK_CPU];
    complex A_fwd_l[MAX_RANK_CPU + 1];

    A_fwd_l[0] = complex::one();
    for (int t = 0; t < rank; t++) {
      A_list_l[t] = b.A[aoff[t]];
      A_fwd_l[t + 1] = A_fwd_l[t] * A_list_l[t];
    }

    complex A_backward_prod = complex::one();
    for (int t = r; t >= 1; t--) {
      dB[t] = A_fwd_l[t] * A_backward_prod;
      A_backward_prod = A_backward_prod * A_list_l[t];
    }
    dB[0] = A_fwd_l[0] * A_backward_prod;

    B = A_fwd_l[rank];
  }
  const KK_FLOAT *ct = b.ctildes + idx_ms_combs * b.ndensitymax;
  const int nd = NDENSITY ? NDENSITY : ndensity;
  for (int p = 0; p < nd; ++p)
    b.rho[p] += B.real_part_product(ct[p]);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairPACEKokkos<DeviceType>::operator() (TagPairPACEComputeFS, const int& ii) const
{
  compute_fs_one(ii);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairPACEKokkos<DeviceType>::compute_fs_one(const int ii) const
{
  const int i = d_ilist[ii + chunk_offset];
  const int mu_i = d_map(type(i));

  // First-touch zeroing of this atom's weight accumulators. ComputeFS runs one
  // thread per atom immediately before ComputeWeights (which accumulates into
  // these via atomic_add), so the slice is cleared here instead of with
  // separate per-chunk deep_copy passes over the full (large) weight arrays.
  // Device: first-touch zeroing of this atom's weight accumulators, so the
  // large weight arrays need no separate per-chunk deep_copy pass. On the
  // host this loses: a bulk deep_copy (issued in compute()) is much faster
  // than per-atom scalar stores, so the host path skips this.
  if constexpr (!host_flag) {
    for (int mu = 0; mu < nelements; mu++) {
      for (int n = 0; n < nradbase; n++)
        weights_rank1(ii, mu, n) = 0.0;
      for (int idx = 0; idx < idx_sph_max; idx++)
        for (int n = 0; n <= nradmax; n++) {
          weights_re(ii, mu, idx, n) = 0.0;
          weights_im(ii, mu, idx, n) = 0.0;
        }
    }
  }

  const KK_FLOAT rho_cut = d_rho_core_cutoff(mu_i);
  const KK_FLOAT drho_cut = d_drho_core_cutoff(mu_i);
  const int ndensity = d_ndensity(mu_i);

  KK_FLOAT evdwl, fcut, dfcut;
  KK_FLOAT evdwl_cut;
  evdwl = fcut = dfcut = 0.0;

  FS_values_and_derivatives(ii, evdwl, mu_i);

  if (is_zbl) {
    if (d_jj_min(ii) != -1) {
      const int mu_jmin = d_mu(ii,d_jj_min(ii));
      KK_FLOAT dcutin = d_dcut_in(mu_i, mu_jmin);
      KK_FLOAT transition_coordinate =  dcutin  - d_d_min(ii); // == cutin - r_min
      cutoff_func_poly(transition_coordinate, dcutin, dcutin, fcut, dfcut);
      dfcut = -dfcut; // invert, because rho_core = cutin - r_min
    } else {
      // no neighbours
      fcut = 1;
      dfcut = 0;
    }
    evdwl_cut = evdwl * fcut + rho_core(ii) * (1 - fcut); // evdwl * fcut + rho_core_uncut  - rho_core_uncut* fcut
    dF_drho_core(ii) = 1 - fcut;
    dF_dfcut(ii) = evdwl * dfcut - rho_core(ii) * dfcut;
  } else {
    inner_cutoff(rho_core(ii), rho_cut, drho_cut, fcut, dfcut);
    dF_drho_core(ii) = evdwl * dfcut + 1;
    evdwl_cut = evdwl * fcut + rho_core(ii);
  }
  for (int p = 0; p < ndensity; ++p)
    dF_drho(ii, p) *= fcut;

  // tally energy contribution
  if (eflag) {
    // E0 shift
    evdwl_cut += d_E0vals(mu_i);
    e_atom(ii) = evdwl_cut;
  }

  if (flag_corerep_factor)
    d_corerep(ii) = 1-fcut;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairPACEKokkos<DeviceType>::operator() (TagPairPACEComputeWeights, const int& iter) const
{
  const int idx_ms_combs = iter / chunk_size;
  const int ii = iter % chunk_size;

  const int i = d_ilist[ii + chunk_offset];
  const int mu_i = d_map(type(i));

  if (idx_ms_combs >= d_idx_ms_combs_count(mu_i)) return;

  const int ndensity = d_ndensity(mu_i);

  const int idx_func = d_idx_funcs(mu_i, idx_ms_combs);
  const int rank = d_rank(mu_i, idx_func);

  // Weights and theta calculation

  if (rank == 1) {
    const int mu = d_mus(mu_i, idx_func, 0);
    const int n = d_ns(mu_i, idx_func, 0);
    KK_FLOAT theta = 0.0;
    for (int p = 0; p < ndensity; ++p) {
      // for rank=1 (r=0) only 1 ms-combination exists (ms_ind=0), so index of func.ctildes is 0..ndensity-1
      theta += dF_drho(ii, p) * d_ctildes(mu_i, idx_ms_combs, p);
    }
    Kokkos::atomic_add(&weights_rank1(ii, mu, n - 1), theta);
  } else { // rank > 1
    KK_FLOAT theta = 0.0;
    for (int p = 0; p < ndensity; ++p)
      theta += dF_drho(ii, p) * d_ctildes(mu_i, idx_ms_combs, p);

    theta *= 0.5; // 0.5 factor due to possible KK_FLOAT counting ???
    for (int t = 0; t < rank; ++t) {
      const int m_t = d_ms_combs(mu_i, idx_ms_combs, t);
      const int factor = (m_t % 2 == 0 ? 1 : -1);
      // dB = product of all factors except t (leave-one-out), recomputed here
      // from A_sph instead of reading a stored dB_flatten array.
      complex dB = complex::one();
      for (int s = 0; s < rank; ++s) {
        if (s == t) continue;
        const int mu_s = d_mus(mu_i, idx_func, s);
        const int n_s = d_ns(mu_i, idx_func, s);
        const int l_s = d_ls(mu_i, idx_func, s);
        const int m_s = d_ms_combs(mu_i, idx_ms_combs, s);
        dB = dB * read_A(ii, mu_s, l_s, m_s, n_s - 1);
      }
      const int mu_t = d_mus(mu_i, idx_func, t);
      const int n_t = d_ns(mu_i, idx_func, t);
      const int l_t = d_ls(mu_i, idx_func, t);
      const int idx = l_t * (l_t + 1) + m_t; // (l, m)
      const int idx_sph = d_idx_sph(idx);
      if (idx_sph >= 0) {
        const complex value = theta * dB;
        Kokkos::atomic_add(&(weights_re(ii, mu_t, idx_sph, n_t - 1)), value.re);
        Kokkos::atomic_add(&(weights_im(ii, mu_t, idx_sph, n_t - 1)), value.im);
      }
      // update -m_t (that could also be positive), because the basis is half_basis
      const int idxm = l_t * (l_t + 1) - m_t; // (l, -m)
      const int idxm_sph = d_idx_sph(idxm);
      if (idxm_sph >= 0) {
        const complex valuem = theta * dB.conj() * (KK_FLOAT)factor;
        Kokkos::atomic_add(&(weights_re(ii, mu_t, idxm_sph, n_t - 1)), valuem.re);
        Kokkos::atomic_add(&(weights_im(ii, mu_t, idxm_sph, n_t - 1)), valuem.im);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   CPU backend: one thread owns atom ii, so the weights scatter is a plain
   accumulation into that atom's own block instead of four atomics per
   (ms-combination, rank) step
------------------------------------------------------------------------- */

template<class DeviceType>
template<int NDENSITY>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairPACEKokkos<DeviceType>::compute_weights_one_cpu(const BasisPtrs &b,
                                                         const int ndensity,
                                                         const int idx_ms_combs) const
{
  const int rank = b.rank[b.idx_funcs[idx_ms_combs]];

  const int *woff  = b.w_off  + idx_ms_combs * b.rankmax;
  const int *wmoff = b.wm_off + idx_ms_combs * b.rankmax;
  const complex *dB = b.dB + b.dB_off[idx_ms_combs];
  const KK_FLOAT *ct = b.ctildes + idx_ms_combs * b.ndensitymax;

  KK_FLOAT theta = 0.0;
  const int nd = NDENSITY ? NDENSITY : ndensity;
  for (int p = 0; p < nd; ++p)
    theta += b.dF[p] * ct[p];

  theta *= 0.5; // 0.5 factor due to possible double counting
  for (int t = 0; t < rank; ++t) {
    // both (l,m) and (l,-m) are updated unconditionally; entries outside the
    // packed triangle land in the trash row
    const complex value = theta * dB[t];
    complex &w = b.w[woff[t]];
    w.re += value.re;
    w.im += value.im;

    const int packed = wmoff[t];
    const KK_FLOAT factor = (packed & 1) ? KK_FLOAT(-1.0) : KK_FLOAT(1.0);
    const complex valuem = theta * dB[t].conj() * factor;
    complex &wm = b.w[packed >> 1];
    wm.re += valuem.re;
    wm.im += valuem.im;
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairPACEKokkos<DeviceType>::compute_derivative_radial(const int ii, const int jj,
    const int mu_j, const int idx_sph, const int l, const complex &ylm,
    const complex (&dylm)[3], const KK_FLOAT rinv, const KK_FLOAT (&r_hat)[3],
    const KK_FLOAT wscale, KK_ACC_FLOAT (&f_ji)[3]) const
{
  for (int n = 0; n < nradmax; n++) {

    // Read and test the (idx_sph, n) weight first: skipping the radial reads
    // and complex products for zero weights avoids needless memory traffic.
    complex w = complex(weights_re(ii, mu_j, idx_sph, n), weights_im(ii, mu_j, idx_sph, n));
    if (w.re == 0.0 && w.im == 0.0) continue;
    // wscale folds in the factor-of-2 that accounts for the -m cases (m > 0)
    w.re *= wscale;
    w.im *= wscale;

    const KK_FLOAT R_over_r = fr(ii, jj, n * (lmax + 1) + l) * rinv;
    const KK_FLOAT DR = dfr(ii, jj, n * (lmax + 1) + l);
    const complex Y_DR = ylm * DR;

    complex grad_phi_nlm[3];
    grad_phi_nlm[0] = Y_DR * r_hat[0] + dylm[0] * R_over_r;
    grad_phi_nlm[1] = Y_DR * r_hat[1] + dylm[1] * R_over_r;
    grad_phi_nlm[2] = Y_DR * r_hat[2] + dylm[2] * R_over_r;
    // real-part multiplication only
    f_ji[0] += w.real_part_product(grad_phi_nlm[0]);
    f_ji[1] += w.real_part_product(grad_phi_nlm[1]);
    f_ji[2] += w.real_part_product(grad_phi_nlm[2]);
  }
}

/* ---------------------------------------------------------------------- */
template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairPACEKokkos<DeviceType>::operator() (TagPairPACEComputeDerivative, const typename Kokkos::TeamPolicy<DeviceType, TagPairPACEComputeDerivative>::member_type& team) const
{
  // Extract the atom number
  int ii = team.team_rank() + team.team_size() * (team.league_rank() %
           ((chunk_size+team.team_size()-1)/team.team_size()));
  if (ii >= chunk_size) return;

  // Extract the neighbor number
  const int jj = team.league_rank() / ((chunk_size+team.team_size()-1)/team.team_size());
  const int ncount = d_ncount(ii);
  if (jj >= ncount) return;

  compute_derivative_one(ii, jj);
}

/* ----------------------------------------------------------------------
   CPU backend: one atom per thread, so the atom's weights block is loaded
   once and reused across its whole neighbor list
------------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairPACEKokkos<DeviceType>::operator() (TagPairPACEComputeDerivativeCPU, const int& ii) const
{
  // Host-only kernel. The body is guarded so the explicit instantiation of
  // PairPACEKokkos<LMPDeviceType> does not compile it for the device, where
  // it would pull in host-only helpers and the LayoutRight assumption.
  // "if constexpr (host_flag)" alone is not enough: PairPACEKokkos<LMPHostType>
  // (pace/kk/host) is also explicitly instantiated in GPU builds, and there
  // host_flag is always true, so the device compiler pass still type-checks
  // this body for that instantiation -- hence the extra LMP_KK_DEVICE_COMPILE
  // preprocessor guard, which is false for both instantiations in that pass.
#ifndef LMP_KK_DEVICE_COMPILE
  if constexpr (host_flag) {
    const int ncount = d_ncount(ii);

    if (!use_batched_cpu) {
      for (int jj = 0; jj < ncount; jj++)
        compute_derivative_one(ii, jj);
      return;
    }

    constexpr int V = PACE_VLEN;
    const int nrl = (lmax + 1) * nradmax;

    const int i = d_ilist[ii + chunk_offset];
    const int itype = type(i);
    const KK_FLOAT scale = d_scale(itype, itype);

    const int fr_sl = 1;              // fr is flat (n*(lmax+1)+l): l-stride
    const int fr_sn = lmax + 1;       // ... and n-stride
    const int dgr_sn = (int) dgr.stride(2);
    const int w_ss = (int) weights.stride(2);
    const int w_sn = (int) weights.stride(3);
    const int wr1_sn = (int) weights_rank1.stride(2);

    const KK_FLOAT *alm_p = alm.data();
    const KK_FLOAT *blm_p = blm.data();
    const KK_FLOAT *cl_p = cl.data();
    const KK_FLOAT *dl_p = dl.data();

    const KK_FLOAT dfcore = dF_drho_core(ii);
    const int jj_min = is_zbl ? (int) d_jj_min(ii) : -1;
    const KK_FLOAT dfc = is_zbl ? dF_dfcut(ii) : 0.0;

    double rxb[V], ryb[V], rzb[V], rinvb[V];
    double fr_b[PACE_BATCH_NRL_MAX * V];
    double dfr_b[PACE_BATCH_NRL_MAX * V];
    double dgr_b[PACE_BATCH_NRB_MAX * V];
    double fxb[V], fyb[V], fzb[V];
    int jidx[V];

    for (int mu = 0; mu < nelements; mu++) {
      const KK_FLOAT *w_mu = reinterpret_cast<const KK_FLOAT *>(&weights(ii, mu, 0, 0));
      const KK_FLOAT *wr1_mu = &weights_rank1(ii, mu, 0);

      int nb = 0;
      for (int jj = 0; jj < ncount; jj++) {
        if ((int) d_mu(ii, jj) != mu) continue;

        rxb[nb] = d_rhats(ii, jj, 0);
        ryb[nb] = d_rhats(ii, jj, 1);
        rzb[nb] = d_rhats(ii, jj, 2);
        rinvb[nb] = 1.0 / d_rnorms(ii, jj);
        jidx[nb] = jj;
        const KK_FLOAT *frp = &fr(ii, jj, 0);
        const KK_FLOAT *dfrp = &dfr(ii, jj, 0);
        for (int l = 0; l <= lmax; l++)
          for (int n = 0; n < nradmax; n++) {
            fr_b[(l * nradmax + n) * V + nb] = frp[l * fr_sl + n * fr_sn];
            dfr_b[(l * nradmax + n) * V + nb] = dfrp[l * fr_sl + n * fr_sn];
          }
        const KK_FLOAT *dgrp = &dgr(ii, jj, 0);
        for (int n = 0; n < nradbase; n++)
          dgr_b[n * V + nb] = dgrp[n * dgr_sn];

        if (++nb < V) continue;

        pace_batched_derivative(lmax, nradmax, nradbase, alm_p, blm_p, cl_p, dl_p,
                                rxb, ryb, rzb, rinvb, fr_b, dfr_b, dgr_b,
                                w_mu, w_ss, w_sn, wr1_mu, wr1_sn, fxb, fyb, fzb);
        for (int lane = 0; lane < V; lane++) {
          const int jjl = jidx[lane];
          const KK_FLOAT fpair = dfcore * dcr(ii, jjl);
          f_ij(ii, jjl, 0) = scale * fxb[lane] + fpair * rxb[lane];
          f_ij(ii, jjl, 1) = scale * fyb[lane] + fpair * ryb[lane];
          f_ij(ii, jjl, 2) = scale * fzb[lane] + fpair * rzb[lane];
          if (jjl == jj_min) {
            f_ij(ii, jjl, 0) += dfc * rxb[lane];
            f_ij(ii, jjl, 1) += dfc * ryb[lane];
            f_ij(ii, jjl, 2) += dfc * rzb[lane];
          }
        }
        nb = 0;
      }

      if (nb > 0) {
        // tail lanes copy lane-0 data so every value stays finite; only the
        // first nb lanes are consumed
        for (int lane = nb; lane < V; lane++) {
          rxb[lane] = rxb[0];
          ryb[lane] = ryb[0];
          rzb[lane] = rzb[0];
          rinvb[lane] = rinvb[0];
        }
        for (int f = 0; f < nrl; f++)
          for (int lane = nb; lane < V; lane++) {
            fr_b[f * V + lane] = fr_b[f * V];
            dfr_b[f * V + lane] = dfr_b[f * V];
          }
        for (int n = 0; n < nradbase; n++)
          for (int lane = nb; lane < V; lane++) dgr_b[n * V + lane] = dgr_b[n * V];

        pace_batched_derivative(lmax, nradmax, nradbase, alm_p, blm_p, cl_p, dl_p,
                                rxb, ryb, rzb, rinvb, fr_b, dfr_b, dgr_b,
                                w_mu, w_ss, w_sn, wr1_mu, wr1_sn, fxb, fyb, fzb);
        for (int lane = 0; lane < nb; lane++) {
          const int jjl = jidx[lane];
          const KK_FLOAT fpair = dfcore * dcr(ii, jjl);
          f_ij(ii, jjl, 0) = scale * fxb[lane] + fpair * rxb[lane];
          f_ij(ii, jjl, 1) = scale * fyb[lane] + fpair * ryb[lane];
          f_ij(ii, jjl, 2) = scale * fzb[lane] + fpair * rzb[lane];
          if (jjl == jj_min) {
            f_ij(ii, jjl, 0) += dfc * rxb[lane];
            f_ij(ii, jjl, 1) += dfc * ryb[lane];
            f_ij(ii, jjl, 2) += dfc * rzb[lane];
          }
        }
      }
    }
  }
#endif
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairPACEKokkos<DeviceType>::compute_derivative_one(const int ii, const int jj) const
{
  const int i = d_ilist[ii + chunk_offset];
  const int itype = type(i);
  const KK_FLOAT scale = d_scale(itype,itype);

  const int mu_j = d_mu(ii, jj);

  // same strength reduction as in compute_ai_one: the (atom, neighbor) and
  // (atom, element) parts of the subscripts are loop invariant here
  const KK_FLOAT * const fr_ij = &fr(ii, jj, 0);
  const KK_FLOAT * const dfr_ij = &dfr(ii, jj, 0);
  const KK_FLOAT * const dgr_ij = &dgr(ii, jj, 0);
  const complex * const w_i = &weights(ii, mu_j, 0, 0);
  const KK_FLOAT * const wr1_i = &weights_rank1(ii, mu_j, 0);
  const int fr_sl = 1;              // fr is flat (n*(lmax+1)+l): l-stride
  const int fr_sn = lmax + 1;       // ... and n-stride
  const int dgr_sn = (int) dgr.stride(2);
  const int w_ss = (int) weights.stride(2);
  const int w_sn = (int) weights.stride(3);
  const int wr1_sn = (int) weights_rank1.stride(2);

  KK_FLOAT r_hat[3];
  r_hat[0] = d_rhats(ii, jj, 0);
  r_hat[1] = d_rhats(ii, jj, 1);
  r_hat[2] = d_rhats(ii, jj, 2);
  const KK_FLOAT r = d_rnorms(ii, jj);
  const KK_FLOAT rinv = 1.0/r;

  KK_ACC_FLOAT f_ji[3];
  f_ji[0] = f_ji[1] = f_ji[2] = 0;

  // for rank = 1
  for (int n = 0; n < nradbase; ++n) {
    KK_FLOAT DGR = dgr_ij[n * dgr_sn] * Y00;
    DGR *= wr1_i[n * wr1_sn];
    f_ji[0] += DGR * r_hat[0];
    f_ji[1] += DGR * r_hat[1];
    f_ji[2] += DGR * r_hat[2];
  }

  // for rank > 1

  // compute plm, dplm, ylm and dylm
  // requires rx^2 + ry^2 + rz^2 = 1 , NO CHECKING IS PERFORMED !!!!!!!!!
  // requires -1 <= rz <= 1 , NO CHECKING IS PERFORMED !!!!!!!!!
  // prefactors include 1/sqrt(2) factor compared to reference

  complex ylm,dylm[3];
  complex phase;
  complex phasem, mphasem1;
  complex dyx, dyy, dyz;
  complex rdy;

  const KK_FLOAT rx = d_rhats(ii, jj, 0);
  const KK_FLOAT ry = d_rhats(ii, jj, 1);
  const KK_FLOAT rz = d_rhats(ii, jj, 2);

  phase.re = rx;
  phase.im = ry;

  KK_FLOAT plm_idx,plm_idx1,plm_idx2;
  KK_FLOAT dplm_idx,dplm_idx1,dplm_idx2;

  plm_idx = plm_idx1 = plm_idx2 = 0.0;
  dplm_idx = dplm_idx1 = dplm_idx2 = 0.0;

  int idx_sph = 0;

  // m = 0
  for (int l = 0; l <= lmax; l++) {
    // const int idx = l * (l + 1);

    if (l == 0) {
      // l=0, m=0
      // plm[0] = Y00/sq1o4pi; //= sq1o4pi;
      plm_idx = Y00; //= 1;
      dplm_idx = 0.0;
    } else if (l == 1) {
      // l=1, m=0
      plm_idx = Y00 * sq3 * rz;
      dplm_idx = Y00 * sq3;
    } else {
      // l>=2, m=0
      plm_idx = alm(idx_sph) * (rz * plm_idx1 + blm(idx_sph) * plm_idx2);
      dplm_idx = alm(idx_sph) * (plm_idx1 + rz * dplm_idx1 + blm(idx_sph) * dplm_idx2);
    }

    ylm.re = plm_idx;
    ylm.im = 0.0;

    dyz.re = dplm_idx;
    rdy.re = dyz.re * rz;

    dylm[0].re = -rdy.re * rx;
    dylm[0].im = 0.0;
    dylm[1].re = -rdy.re * ry;
    dylm[1].im = 0.0;
    dylm[2].re = dyz.re - rdy.re * rz;
    dylm[2].im = 0;

    if constexpr (host_flag) {
      const KK_FLOAT * const r_nl = fr_ij + l * fr_sl;
      const KK_FLOAT * const dr_nl = dfr_ij + l * fr_sl;
      const complex * const w_nl = w_i + idx_sph * w_ss;
      for (int n = 0; n < nradmax; n++) {
        const KK_FLOAT R_over_r = r_nl[n * fr_sn] * rinv;
        const KK_FLOAT DR = dr_nl[n * fr_sn];
        const complex Y_DR = ylm * DR;
        complex w = w_nl[n * w_sn];
        complex grad_phi_nlm[3];
        grad_phi_nlm[0] = Y_DR * r_hat[0] + dylm[0] * R_over_r;
        grad_phi_nlm[1] = Y_DR * r_hat[1] + dylm[1] * R_over_r;
        grad_phi_nlm[2] = Y_DR * r_hat[2] + dylm[2] * R_over_r;
        f_ji[0] += w.real_part_product(grad_phi_nlm[0]);
        f_ji[1] += w.real_part_product(grad_phi_nlm[1]);
        f_ji[2] += w.real_part_product(grad_phi_nlm[2]);
      }
    } else {
      compute_derivative_radial(ii, jj, mu_j, idx_sph, l, ylm, dylm, rinv, r_hat, 1.0, f_ji);
    }

    plm_idx2 = plm_idx1;
    dplm_idx2 = dplm_idx1;

    plm_idx1 = plm_idx;
    dplm_idx1 = dplm_idx;

    idx_sph++;
  }

  plm_idx = plm_idx1 = plm_idx2 = 0.0;
  dplm_idx = dplm_idx1 = dplm_idx2 = 0.0;

  // m = 1
  for (int l = 1; l <= lmax; l++) {
    // const int idx = l * (l + 1) + 1; // (l, 1)

    if (l == 1) {
      // l=1, m=1
      plm_idx = -sq3o2 * Y00;
      dplm_idx = 0.0;
    } else if (l == 2) {
      const KK_FLOAT t = dl(l) * plm_idx1;
      plm_idx = t * rz;
      dplm_idx = t;
    } else {
      plm_idx = alm(idx_sph) * (rz * plm_idx1 + blm(idx_sph) * plm_idx2);
      dplm_idx = alm(idx_sph) * (plm_idx1 + rz * dplm_idx1 + blm(idx_sph) * dplm_idx2);
    }

    ylm = phase * plm_idx;

    dyx.re = plm_idx;
    dyx.im = 0.0;
    dyy.re = 0.0;
    dyy.im = plm_idx;
    dyz.re = phase.re * dplm_idx;
    dyz.im = phase.im * dplm_idx;

    rdy.re = rx * dyx.re + +rz * dyz.re;
    rdy.im = ry * dyy.im + rz * dyz.im;

    dylm[0].re = dyx.re - rdy.re * rx;
    dylm[0].im = -rdy.im * rx;
    dylm[1].re = -rdy.re * ry;
    dylm[1].im = dyy.im - rdy.im * ry;
    dylm[2].re = dyz.re - rdy.re * rz;
    dylm[2].im = dyz.im - rdy.im * rz;

    if constexpr (host_flag) {
      const KK_FLOAT * const r_nl = fr_ij + l * fr_sl;
      const KK_FLOAT * const dr_nl = dfr_ij + l * fr_sl;
      const complex * const w_nl = w_i + idx_sph * w_ss;
      for (int n = 0; n < nradmax; n++) {
        const KK_FLOAT R_over_r = r_nl[n * fr_sn] * rinv;
        const KK_FLOAT DR = dr_nl[n * fr_sn];
        const complex Y_DR = ylm * DR;
        complex w = w_nl[n * w_sn];
          w.re *= 2.0;
          w.im *= 2.0;
        complex grad_phi_nlm[3];
        grad_phi_nlm[0] = Y_DR * r_hat[0] + dylm[0] * R_over_r;
        grad_phi_nlm[1] = Y_DR * r_hat[1] + dylm[1] * R_over_r;
        grad_phi_nlm[2] = Y_DR * r_hat[2] + dylm[2] * R_over_r;
        f_ji[0] += w.real_part_product(grad_phi_nlm[0]);
        f_ji[1] += w.real_part_product(grad_phi_nlm[1]);
        f_ji[2] += w.real_part_product(grad_phi_nlm[2]);
      }
    } else {
      compute_derivative_radial(ii, jj, mu_j, idx_sph, l, ylm, dylm, rinv, r_hat, 2.0, f_ji);
    }

    plm_idx2 = plm_idx1;
    dplm_idx2 = dplm_idx1;

    plm_idx1 = plm_idx;
    dplm_idx1 = dplm_idx;

    idx_sph++;
  }

  plm_idx = plm_idx1 = plm_idx2 = 0.0;
  dplm_idx = dplm_idx1 = dplm_idx2 = 0.0;

  KK_FLOAT plm_mm1_mm1 = -sq3o2 * Y00; // (1, 1)

  // m > 1
  phasem = phase;
  for (int m = 2; m <= lmax; m++) {

    mphasem1.re = phasem.re * KK_FLOAT(m);
    mphasem1.im = phasem.im * KK_FLOAT(m);
    phasem = phasem * phase;

    for (int l = m; l <= lmax; l++) {
      // const int idx = l * (l + 1) + m;

      if (l == m) {
        plm_idx = cl(l) * plm_mm1_mm1; // (m+1, m)
        dplm_idx = 0.0;
        plm_mm1_mm1 = plm_idx;
      } else if (l == (m + 1)) {
        const KK_FLOAT t = dl(l) * plm_mm1_mm1; // (m - 1, m - 1)
        plm_idx = t * rz; // (m, m)
        dplm_idx = t;
      } else {
        plm_idx = alm(idx_sph) * (rz * plm_idx1 + blm(idx_sph) * plm_idx2);
        dplm_idx = alm(idx_sph) * (plm_idx1 + rz * dplm_idx1 + blm(idx_sph) * dplm_idx2);
      }

      ylm.re = phasem.re * plm_idx;
      ylm.im = phasem.im * plm_idx;

      dyx = mphasem1 * plm_idx;
      dyy.re = -dyx.im;
      dyy.im = dyx.re;
      dyz = phasem * dplm_idx;

      rdy.re = rx * dyx.re + ry * dyy.re + rz * dyz.re;
      rdy.im = rx * dyx.im + ry * dyy.im + rz * dyz.im;

      dylm[0].re = dyx.re - rdy.re * rx;
      dylm[0].im = dyx.im - rdy.im * rx;
      dylm[1].re = dyy.re - rdy.re * ry;
      dylm[1].im = dyy.im - rdy.im * ry;
      dylm[2].re = dyz.re - rdy.re * rz;
      dylm[2].im = dyz.im - rdy.im * rz;

      if constexpr (host_flag) {
        const KK_FLOAT * const r_nl = fr_ij + l * fr_sl;
        const KK_FLOAT * const dr_nl = dfr_ij + l * fr_sl;
        const complex * const w_nl = w_i + idx_sph * w_ss;
        for (int n = 0; n < nradmax; n++) {
          const KK_FLOAT R_over_r = r_nl[n * fr_sn] * rinv;
          const KK_FLOAT DR = dr_nl[n * fr_sn];
          const complex Y_DR = ylm * DR;
          complex w = w_nl[n * w_sn];
            w.re *= 2.0;
            w.im *= 2.0;
          complex grad_phi_nlm[3];
          grad_phi_nlm[0] = Y_DR * r_hat[0] + dylm[0] * R_over_r;
          grad_phi_nlm[1] = Y_DR * r_hat[1] + dylm[1] * R_over_r;
          grad_phi_nlm[2] = Y_DR * r_hat[2] + dylm[2] * R_over_r;
          f_ji[0] += w.real_part_product(grad_phi_nlm[0]);
          f_ji[1] += w.real_part_product(grad_phi_nlm[1]);
          f_ji[2] += w.real_part_product(grad_phi_nlm[2]);
        }
      } else {
        compute_derivative_radial(ii, jj, mu_j, idx_sph, l, ylm, dylm, rinv, r_hat, 2.0, f_ji);
      }

      plm_idx2 = plm_idx1;
      dplm_idx2 = dplm_idx1;

      plm_idx1 = plm_idx;
      dplm_idx1 = dplm_idx;

      idx_sph++;
    }
  }

  // hard-core repulsion
  const KK_FLOAT fpair = dF_drho_core(ii) * dcr(ii,jj);
  f_ij(ii, jj, 0) = scale * f_ji[0] + fpair * r_hat[0];
  f_ij(ii, jj, 1) = scale * f_ji[1] + fpair * r_hat[1];
  f_ij(ii, jj, 2) = scale * f_ji[2] + fpair * r_hat[2];

  if (is_zbl) {
    if (jj==d_jj_min(ii)) {
      // DCRU = 1.0
      f_ij(ii, jj, 0) += dF_dfcut(ii) * r_hat[0];
      f_ij(ii, jj, 1) += dF_dfcut(ii) * r_hat[1];
      f_ij(ii, jj, 2) += dF_dfcut(ii) * r_hat[2];
    }
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG, int EVFLAG>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairPACEKokkos<DeviceType>::operator() (TagPairPACEComputeForce<NEIGHFLAG,EVFLAG>, const int& ii, EV_FLOAT& ev) const
{
  // The f array is duplicated for OpenMP, atomic for GPU, and neither for Serial
  const auto v_f = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_f),decltype(ndup_f)>::get(dup_f,ndup_f);
  const auto a_f = v_f.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  const int i = d_ilist[ii + chunk_offset];
  const int itype = type(i);
  const KK_FLOAT scale = d_scale(itype,itype);

  const int ncount = d_ncount(ii);

  KK_ACC_FLOAT fitmp[3] = {0.0,0.0,0.0};
  for (int jj = 0; jj < ncount; jj++) {
    int j = d_nearest(ii,jj);

    KK_FLOAT r_hat[3];
    r_hat[0] = d_rhats(ii, jj, 0);
    r_hat[1] = d_rhats(ii, jj, 1);
    r_hat[2] = d_rhats(ii, jj, 2);
    const KK_FLOAT r = d_rnorms(ii, jj);
    const KK_FLOAT delx = -r_hat[0]*r;
    const KK_FLOAT dely = -r_hat[1]*r;
    const KK_FLOAT delz = -r_hat[2]*r;

    const KK_FLOAT fpairx = f_ij(ii, jj, 0);
    const KK_FLOAT fpairy = f_ij(ii, jj, 1);
    const KK_FLOAT fpairz = f_ij(ii, jj, 2);

    fitmp[0] += fpairx;
    fitmp[1] += fpairy;
    fitmp[2] += fpairz;
    a_f(j,0) -= fpairx;
    a_f(j,1) -= fpairy;
    a_f(j,2) -= fpairz;

    // tally per-atom virial contribution
    if (EVFLAG && vflag_either)
      v_tally_xyz<NEIGHFLAG>(ev, i, j, fpairx, fpairy, fpairz, delx, dely, delz);
  }

  a_f(i,0) += fitmp[0];
  a_f(i,1) += fitmp[1];
  a_f(i,2) += fitmp[2];

  // tally energy contribution
  if (EVFLAG && eflag_either) {
    const KK_FLOAT evdwl = scale*e_atom(ii);
    //ev_tally_full(i, 2.0 * evdwl, 0.0, 0.0, 0.0, 0.0, 0.0);
    if (eflag_global) ev.evdwl += evdwl;
    if (eflag_atom) d_eatom[i] += evdwl;
  }
}

template<class DeviceType>
template<int NEIGHFLAG, int EVFLAG>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairPACEKokkos<DeviceType>::operator() (TagPairPACEComputeForce<NEIGHFLAG,EVFLAG>,const int& ii) const {
  EV_FLOAT ev;
  this->template operator()<NEIGHFLAG,EVFLAG>(TagPairPACEComputeForce<NEIGHFLAG,EVFLAG>(), ii, ev);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairPACEKokkos<DeviceType>::v_tally_xyz(EV_FLOAT &ev, const int &i, const int &j,
      const KK_FLOAT &fx, const KK_FLOAT &fy, const KK_FLOAT &fz,
      const KK_FLOAT &delx, const KK_FLOAT &dely, const KK_FLOAT &delz) const
{
  // The vatom array is duplicated for OpenMP, atomic for GPU, and neither for Serial

  auto v_vatom = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_vatom),decltype(ndup_vatom)>::get(dup_vatom,ndup_vatom);
  auto a_vatom = v_vatom.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  const KK_FLOAT v0 = delx*fx;
  const KK_FLOAT v1 = dely*fy;
  const KK_FLOAT v2 = delz*fz;
  const KK_FLOAT v3 = delx*fy;
  const KK_FLOAT v4 = delx*fz;
  const KK_FLOAT v5 = dely*fz;

  if (vflag_global) {
    ev.v[0] += v0;
    ev.v[1] += v1;
    ev.v[2] += v2;
    ev.v[3] += v3;
    ev.v[4] += v4;
    ev.v[5] += v5;
  }

  if (vflag_atom) {
    a_vatom(i,0) += 0.5*v0;
    a_vatom(i,1) += 0.5*v1;
    a_vatom(i,2) += 0.5*v2;
    a_vatom(i,3) += 0.5*v3;
    a_vatom(i,4) += 0.5*v4;
    a_vatom(i,5) += 0.5*v5;
    a_vatom(j,0) += 0.5*v0;
    a_vatom(j,1) += 0.5*v1;
    a_vatom(j,2) += 0.5*v2;
    a_vatom(j,3) += 0.5*v3;
    a_vatom(j,4) += 0.5*v4;
    a_vatom(j,5) += 0.5*v5;
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairPACEKokkos<DeviceType>::pre_compute_harmonics(int lmax)
{
  auto h_idx_sph = Kokkos::create_mirror_view(d_idx_sph);
  auto h_alm = Kokkos::create_mirror_view(alm);
  auto h_blm = Kokkos::create_mirror_view(blm);
  auto h_cl = Kokkos::create_mirror_view(cl);
  auto h_dl = Kokkos::create_mirror_view(dl);

  Kokkos::deep_copy(h_idx_sph,-1);

  int idx_sph = 0;
  for (int m = 0; m <= lmax; m++) {
    const double msq = m * m;
    for (int l = m; l <= lmax; l++) {
      const int idx = l * (l + 1) + m; // (l, m)
      h_idx_sph(idx) = idx_sph;

      double a = 0.0;
      double b = 0.0;

      if (l > 1 && l != m) {
        const double lsq = l * l;
        const double ld = 2 * l;
        const double l1 = (4 * lsq - 1);
        const double l2 = lsq - ld + 1;

        a = sqrt((double(l1)) / (double(lsq - msq)));
        b = -sqrt((double(l2 - msq)) / (double(4 * l2 - 1)));
      }
      h_alm(idx_sph) = a;
      h_blm(idx_sph) = b;
      idx_sph++;
    }
  }
  idx_sph_max = idx_sph;

  for (int l = 1; l <= lmax; l++) {
    h_cl(l) = -sqrt(1.0 + 0.5 / (double(l)));
    h_dl(l) = sqrt(double(2 * (l - 1) + 3));
  }

  Kokkos::deep_copy(d_idx_sph, h_idx_sph);

  // Host copy of the same table, as int, with the "no such (l,m)" sentinel
  // pointing at a trash slot appended to weights instead of -1.  That turns
  // the two data-dependent guards in the weights scatter, which alternate
  // with the sign of m and are essentially unpredictable, into unconditional
  // stores.
  auto h_idx_sph_cpu = Kokkos::create_mirror_view(d_idx_sph_cpu);
  for (int idx = 0; idx < (lmax + 1) * (lmax + 1); idx++) {
    const int v = (int) h_idx_sph(idx);
    h_idx_sph_cpu(idx) = (v >= 0) ? v : idx_sph_max;
  }
  Kokkos::deep_copy(d_idx_sph_cpu, h_idx_sph_cpu);
  Kokkos::deep_copy(alm, h_alm);
  Kokkos::deep_copy(blm, h_blm);
  Kokkos::deep_copy(cl, h_cl);
  Kokkos::deep_copy(dl, h_dl);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairPACEKokkos<DeviceType>::cutoff_func_poly(const KK_FLOAT r, const KK_FLOAT r_in, const KK_FLOAT delta_in, KK_FLOAT &fc, KK_FLOAT &dfc) const
{
  if (r <= r_in-delta_in) {
    fc = 1;
    dfc = 0;
  } else if (r >= r_in ) {
    fc = 0;
    dfc = 0;
  } else {
    KK_FLOAT x = 1 - 2 * (1 + (r - r_in) / delta_in);
    // explicit integer powers (avoid pow(): ~hundreds of cycles on GPU)
    const KK_FLOAT x2 = x * x;
    const KK_FLOAT x3 = x2 * x;
    const KK_FLOAT x4 = x2 * x2;
    const KK_FLOAT x5 = x4 * x;
    fc = 0.5 + 7.5 / 2. * (x / 4. - x3 / 6. + x5 / 20.);
    dfc = -7.5 / delta_in * (0.25 - x2 / 2.0 + x4 / 4.);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairPACEKokkos<DeviceType>::Fexp(const KK_FLOAT x, const KK_FLOAT m, KK_FLOAT &F, KK_FLOAT &DF) const
{
  const KK_FLOAT w = 1.e6;
  const KK_FLOAT eps = 1e-10;

  const KK_FLOAT lambda = pow(1.0 / w, m - 1.0);
  if (abs(x) > eps) {
    KK_FLOAT g;
    const KK_FLOAT a = abs(x);
    const KK_FLOAT am = pow(a, m);
    const KK_FLOAT wa = w * a;
    const KK_FLOAT w3x3 = wa * wa * wa; // cube (avoid pow())
    const KK_FLOAT sign_factor = (signbit(x) ? -1 : 1);
    if (w3x3 > 30.0)
        g = 0.0;
    else
        g = exp(-w3x3);

    const KK_FLOAT omg = 1.0 - g;
    F = sign_factor * (omg * am + lambda * g * a);
    const KK_FLOAT dg = -3.0 * w * w * w * a * a * g;
    DF = m * pow(a, m - 1.0) * omg - am * dg + lambda * dg * a + lambda * g;
  } else {
    F = lambda * x;
    DF = lambda;
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairPACEKokkos<DeviceType>::FexpShiftedScaled(const KK_FLOAT rho, const KK_FLOAT mexp, KK_FLOAT &F, KK_FLOAT &DF) const
{
  const KK_FLOAT eps = 1e-10;

  if (abs(mexp - 1.0) < eps) {
    F = rho;
    DF = 1;
  } else {
    const KK_FLOAT a = abs(rho);
    const KK_FLOAT exprho = exp(-a);
    const KK_FLOAT nx = 1. / mexp;
    const KK_FLOAT xoff = pow(nx, (nx / (1.0 - nx))) * exprho;
    const KK_FLOAT yoff = pow(nx, (1 / (1.0 - nx))) * exprho;
    const KK_FLOAT sign_factor = (signbit(rho) ? -1 : 1);
    F = sign_factor * (pow(xoff + a, mexp) - yoff);
    DF = yoff + mexp * (-xoff + 1.0) * pow(xoff + a, mexp - 1.);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairPACEKokkos<DeviceType>::inner_cutoff(const KK_FLOAT rho_core, const KK_FLOAT rho_cut, const KK_FLOAT drho_cut,
                                     KK_FLOAT &fcut, KK_FLOAT &dfcut) const
{
  KK_FLOAT rho_low = rho_cut - drho_cut;
  if (rho_core >= rho_cut) {
    fcut = 0;
    dfcut = 0;
  } else if (rho_core <= rho_low) {
    fcut = 1;
    dfcut = 0;
  } else {
    cutoff_func_poly(rho_core, rho_cut, drho_cut, fcut, dfcut);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairPACEKokkos<DeviceType>::FS_values_and_derivatives(const int ii, KK_FLOAT &evdwl, const int mu_i) const
{
  KK_FLOAT F, DF = 0;
  int npoti = d_npoti(mu_i);
  int ndensity = d_ndensity(mu_i);
  for (int p = 0; p < ndensity; p++) {
    const KK_FLOAT wpre = d_wpre(mu_i, p);
    const KK_FLOAT mexp = d_mexp(mu_i, p);

    if (npoti == FS)
      Fexp(rhos(ii, p), mexp, F, DF);
    else if (npoti == FS_SHIFTEDSCALED)
      FexpShiftedScaled(rhos(ii, p), mexp, F, DF);

    evdwl += F * wpre; // * weight (wpre)
    dF_drho(ii, p) = DF * wpre; // * weight (wpre)
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairPACEKokkos<DeviceType>::evaluate_splines(const int ii, const int jj, KK_FLOAT r,
                                                  int /*nradbase_c*/, int /*nradial_c*/,
                                                  int mu_i, int mu_j) const
{
  auto &spline_gk = k_splines_gk.template view<DeviceType>()(mu_i, mu_j);
  auto &spline_rnl = k_splines_rnl.template view<DeviceType>()(mu_i, mu_j);
  auto &spline_hc = k_splines_hc.template view<DeviceType>()(mu_i, mu_j);

  spline_gk.calcSplines(ii, jj, r, gr, dgr);

  // fr/dfr use the spline's flat (n*(lmax+1)+l) function order, so the rnl
  // spline writes them directly -- no separate d_values buffer + copy pass.
  // (Transposing to n-innermost for the host gather loops was measured and
  // was a net loss: the extra copy pass costs more than the unit stride buys.)
  spline_rnl.calcSplines(ii, jj, r, fr, dfr);

  // the hard-core repulsion is always taken from its (single-function) spline
  spline_hc.calcSplines(ii, jj, r, d_values, d_derivatives);
  cr(ii, jj) = d_values(ii, jj, 0);
  dcr(ii, jj) = d_derivatives(ii, jj, 0);
}

/* ---------------------------------------------------------------------- */
template<class DeviceType>
void PairPACEKokkos<DeviceType>::SplineInterpolatorKokkos::operator=(const SplineInterpolator &spline) {
    cutoff = spline.cutoff;
    deltaSplineBins = spline.deltaSplineBins;
    ntot = spline.ntot;
    nlut = spline.nlut;
    invrscalelookup = spline.invrscalelookup;
    rscalelookup = spline.rscalelookup;
    num_of_functions = spline.num_of_functions;

    lookupTable = t_ace_3d4_lr("lookupTable", ntot+1, num_of_functions);
    auto h_lookupTable = Kokkos::create_mirror_view(lookupTable);
    for (int i = 0; i < ntot+1; i++)
        for (int j = 0; j < num_of_functions; j++)
            for (int k = 0; k < 4; k++)
                h_lookupTable(i, j, k) = spline.lookupTable(i, j, k);
    Kokkos::deep_copy(lookupTable, h_lookupTable);
}
/* ---------------------------------------------------------------------- */
template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairPACEKokkos<DeviceType>::SplineInterpolatorKokkos::calcSplines(const int ii, const int jj, const KK_FLOAT r, const t_ace_3d &d_values, const t_ace_3d &d_derivatives) const
{
  KK_FLOAT wl;
  KK_FLOAT c[4];
  KK_FLOAT x = r * rscalelookup;
  int nl = static_cast<int>(floor(x));

  if (nl <= 0)
    Kokkos::abort("Encountered very small distance. Stopping.");

  if (nl < nlut) {
    wl = x - KK_FLOAT(nl);
    for (int func_id = 0; func_id < num_of_functions; func_id++) {
      for (int idx = 0; idx < 4; idx++)
        c[idx] = lookupTable(nl, func_id, idx);
      // Horner: three FMAs for the value, two for the derivative
      d_values(ii, jj, func_id) = c[0] + wl * (c[1] + wl * (c[2] + wl * c[3]));
      d_derivatives(ii, jj, func_id) =
          (c[1] + wl * (2.0 * c[2] + wl * (3.0 * c[3]))) * rscalelookup;
    }
  } else { // fill with zeroes
    for (int func_id = 0; func_id < num_of_functions; func_id++) {
      d_values(ii, jj, func_id) = 0.0;
      d_derivatives(ii, jj, func_id) = 0.0;
    }
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<class TagStyle>
void PairPACEKokkos<DeviceType>::check_team_size_for(int inum, int &team_size, int vector_length) {
  int team_size_max;

  team_size_max = Kokkos::TeamPolicy<DeviceType,TagStyle>(inum,Kokkos::AUTO).team_size_max(*this,Kokkos::ParallelForTag());

  if (team_size*vector_length > team_size_max)
    team_size = team_size_max/vector_length;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<class TagStyle>
void PairPACEKokkos<DeviceType>::check_team_size_reduce(int inum, int &team_size, int vector_length) {
  int team_size_max;

  team_size_max = Kokkos::TeamPolicy<DeviceType,TagStyle>(inum,Kokkos::AUTO).team_size_max(*this,Kokkos::ParallelReduceTag());

  if (team_size*vector_length > team_size_max)
    team_size = team_size_max/vector_length;
}

template<class DeviceType>
template<typename scratch_type>
int PairPACEKokkos<DeviceType>::scratch_size_helper(int values_per_team) {
  typedef Kokkos::View<scratch_type*, Kokkos::DefaultExecutionSpace::scratch_memory_space, Kokkos::MemoryTraits<Kokkos::Unmanaged> > ScratchViewType;

  return ScratchViewType::shmem_size(values_per_team);
}

/* ----------------------------------------------------------------------
   memory usage of arrays
------------------------------------------------------------------------- */

template<class DeviceType>
double PairPACEKokkos<DeviceType>::memory_usage()
{
  double bytes = 0;

  bytes += MemKK::memory_usage(A_rank1);
  // dual layout: only one of these sets is allocated, so both are counted and
  // the unallocated one contributes nothing
  bytes += MemKK::memory_usage(A_sph_re);
  bytes += MemKK::memory_usage(A_sph_im);
  bytes += MemKK::memory_usage(A_sph);
  bytes += MemKK::memory_usage(A);
  bytes += MemKK::memory_usage(dB_flatten);
  bytes += MemKK::memory_usage(e_atom);
  bytes += MemKK::memory_usage(rhos);
  bytes += MemKK::memory_usage(dF_drho);
  bytes += MemKK::memory_usage(weights_re);
  bytes += MemKK::memory_usage(weights_im);
  bytes += MemKK::memory_usage(weights);
  bytes += MemKK::memory_usage(weights_rank1);
  bytes += MemKK::memory_usage(rho_core);
  bytes += MemKK::memory_usage(dF_drho_core);
  bytes += MemKK::memory_usage(dF_dfcut);
  bytes += MemKK::memory_usage(d_corerep);
  bytes += MemKK::memory_usage(fr);
  bytes += MemKK::memory_usage(dfr);
  bytes += MemKK::memory_usage(gr);
  bytes += MemKK::memory_usage(dgr);
  bytes += MemKK::memory_usage(d_values);
  bytes += MemKK::memory_usage(d_derivatives);
  bytes += MemKK::memory_usage(cr);
  bytes += MemKK::memory_usage(dcr);
  bytes += MemKK::memory_usage(d_ncount);
  bytes += MemKK::memory_usage(d_mu);
  bytes += MemKK::memory_usage(d_rhats);
  bytes += MemKK::memory_usage(d_rnorms);
  bytes += MemKK::memory_usage(d_d_min);
  bytes += MemKK::memory_usage(d_jj_min);
  bytes += MemKK::memory_usage(d_nearest);
  bytes += MemKK::memory_usage(f_ij);
  bytes += MemKK::memory_usage(d_rho_core_cutoff);
  bytes += MemKK::memory_usage(d_drho_core_cutoff);
  bytes += MemKK::memory_usage(d_E0vals);
  bytes += MemKK::memory_usage(d_ndensity);
  bytes += MemKK::memory_usage(d_npoti);
  bytes += MemKK::memory_usage(d_wpre);
  bytes += MemKK::memory_usage(d_mexp);
  bytes += MemKK::memory_usage(d_idx_ms_combs_count);
  bytes += MemKK::memory_usage(d_rank);
  bytes += MemKK::memory_usage(d_num_ms_combs);
  bytes += MemKK::memory_usage(d_idx_funcs);
  bytes += MemKK::memory_usage(d_mus);
  bytes += MemKK::memory_usage(d_ns);
  bytes += MemKK::memory_usage(d_ls);
  bytes += MemKK::memory_usage(d_ms_combs);
  bytes += MemKK::memory_usage(d_ctildes);
  bytes += MemKK::memory_usage(d_idx_sph);
  bytes += MemKK::memory_usage(d_idx_sph_cpu);
  bytes += MemKK::memory_usage(alm);
  bytes += MemKK::memory_usage(blm);
  bytes += MemKK::memory_usage(cl);
  bytes += MemKK::memory_usage(dl);

  if (k_splines_gk.view_host().data()) {
    for (int i = 0; i < nelements; i++) {
      for (int j = 0; j < nelements; j++) {
        bytes += k_splines_gk.view_host()(i, j).memory_usage();
        bytes += k_splines_rnl.view_host()(i, j).memory_usage();
        bytes += k_splines_hc.view_host()(i, j).memory_usage();
      }
    }
  }

  return bytes;
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class PairPACEKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairPACEKokkos<LMPHostType>;
#endif
}
