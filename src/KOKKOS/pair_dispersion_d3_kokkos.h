/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(dispersion/d3/kk,PairDispersionD3Kokkos<LMPDeviceType>);
PairStyle(dispersion/d3/kk/device,PairDispersionD3Kokkos<LMPDeviceType>);
PairStyle(dispersion/d3/kk/host,PairDispersionD3Kokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_DISPERSION_D3_KOKKOS_H
#define LMP_PAIR_DISPERSION_D3_KOKKOS_H

#include "kokkos_base.h"
#include "math_special_kokkos.h"
#include "pair_dispersion_d3.h"
#include "pair_kokkos.h"

#include <limits>

namespace LAMMPS_NS {

using DispersionD3::AUTOANG;
using DispersionD3::AUTOANG6;
using DispersionD3::AUTOEV;
using DispersionD3::K1;
using DispersionD3::K3;

/* ---------------------------------------------------------------------- */
//    Functor to initialize cn and dc6 arrays
/* ---------------------------------------------------------------------- */

template <class DeviceType>
struct PairDispersionD3InitializeFunctor {
  typedef ArrayTypes<DeviceType> AT;

  typename AT::t_kkfloat_1d d_cn;
  typename AT::t_kkfloat_1d d_dc6;

  // NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(const int &i) const {
    d_cn(i) = 0.0;
    d_dc6(i) = 0.0;
  }
};

template <class DeviceType>
struct PairDispersionD3PackForwardCommFunctor {
  typedef ArrayTypes<DeviceType> AT;

  typename AT::t_kkfloat_1d d_arr;
  typename AT::t_int_1d d_sendlist;
  typename AT::t_double_1d_um v_buf;

  // NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(const int &i) const {
    int j = d_sendlist(i);
    v_buf(i) = d_arr(j);
  }
};

template <class DeviceType>
struct PairDispersionD3UnpackForwardCommFunctor {
  typedef ArrayTypes<DeviceType> AT;

  typename AT::t_kkfloat_1d d_arr;
  int first;
  typename AT::t_double_1d_um v_buf;

  // NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(const int &i) const {
    d_arr(i + first) = v_buf(i);
  }
};

template <class DeviceType>
struct PairDispersionD3PackReverseCommFunctor {
  typedef ArrayTypes<DeviceType> AT;

  typename AT::t_kkfloat_1d d_arr;
  int first;
  typename AT::t_double_1d_um v_buf;

  // NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(const int &i) const {
    v_buf(i) = d_arr(i + first);
  }
};

template <class DeviceType>
struct PairDispersionD3UnpackReverseCommFunctor {
  typedef ArrayTypes<DeviceType> AT;

  typename AT::t_kkfloat_1d d_arr;
  typename AT::t_int_1d d_recvlist;
  typename AT::t_double_1d_um v_buf;

  // NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(const int &i) const {
    int j = d_recvlist(i);
    d_arr(j) += v_buf(i);
  }
};

template <class DeviceType, int NEIGHFLAG, int NEWTON_PAIR>
struct PairDispersionD3CoordinationNumberKernel {
  typedef ArrayTypes<DeviceType> AT;
  using DUP = NeedDup_v<NEIGHFLAG,DeviceType>;
  using ScatterAccess = std::conditional_t<
      std::is_same_v<DUP,Kokkos::Experimental::ScatterDuplicated>,
      Kokkos::Experimental::ScatterNonAtomic,
      Kokkos::Experimental::ScatterAtomic>;

  typename AT::t_kkfloat_1d_3_lr_randomread x;
  typename AT::t_int_1d_randomread type;
  typename AT::t_kkfloat_1d d_rcov;
  typename AT::t_kkfloat_1d d_cn;
  typename AT::t_int_1d d_ilist;
  typename AT::t_int_1d d_numneigh;
  typename AT::t_neighbors_2d d_neighbors;
  int nlocal;
  KK_FLOAT cn_thr;

  KKScatterView<KK_FLOAT*, typename DAT::t_kkfloat_1d::array_layout,
                typename KKDevice<DeviceType>::value, KKScatterSum, DUP> dup_cn;

  PairDispersionD3CoordinationNumberKernel(
      const typename AT::t_kkfloat_1d_3_lr_randomread &x_in,
      const typename AT::t_int_1d_randomread &type_in,
      const typename AT::t_kkfloat_1d &d_rcov_in,
      const typename AT::t_kkfloat_1d &d_cn_in,
      const typename AT::t_int_1d &d_ilist_in,
      const typename AT::t_int_1d &d_numneigh_in,
      const typename AT::t_neighbors_2d &d_neighbors_in,
      int nlocal_in, KK_FLOAT cn_thr_in)
    : x(x_in), type(type_in), d_rcov(d_rcov_in), d_cn(d_cn_in),
      d_ilist(d_ilist_in), d_numneigh(d_numneigh_in),
      d_neighbors(d_neighbors_in), nlocal(nlocal_in), cn_thr(cn_thr_in)
  {
    dup_cn = Kokkos::Experimental::create_scatter_view<KKScatterSum, DUP>(d_cn);
  }

  void contribute() {
    Kokkos::Experimental::contribute(d_cn, dup_cn);
  }

  // NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(const int &ii) const {
    auto a_cn = dup_cn.template access<ScatterAccess>();

    int i = d_ilist(ii);
    int itype = type(i);
    int jnum = d_numneigh(i);

    const KK_FLOAT xtmp = x(i, 0);
    const KK_FLOAT ytmp = x(i, 1);
    const KK_FLOAT ztmp = x(i, 2);

    KK_FLOAT cn_i = 0.0;

    for (int jj = 0; jj < jnum; jj++) {
      int j = d_neighbors(i,jj);
      j &= NEIGHMASK;
      int jtype = type(j);

      const KK_FLOAT delx = xtmp - x(j, 0);
      const KK_FLOAT dely = ytmp - x(j, 1);
      const KK_FLOAT delz = ztmp - x(j, 2);
      const KK_FLOAT rsq = delx * delx + dely * dely + delz * delz;

      if (rsq > cn_thr) continue;

      const KK_FLOAT rr = Kokkos::sqrt(rsq);
      const KK_FLOAT rcov_ij = (d_rcov(itype) + d_rcov(jtype)) * AUTOANG;
      const KK_FLOAT cn_ij = 1.0 / (1.0 + Kokkos::exp(-K1 * ((rcov_ij / rr) - 1.0)));

      cn_i += cn_ij;
      if (NEIGHFLAG != FULL && (NEWTON_PAIR || j < nlocal)) a_cn(j) += cn_ij;
    }

    a_cn(i) += cn_i;
  }
};


/* ----------------------------------------------------------------------
   Kernel A: compute energy/force  and dC6
------------------------------------------------------------------------- */

template <class DeviceType, int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
struct PairDispersionD3KernelA {
  typedef ArrayTypes<DeviceType> AT;
  using value_type = EV_FLOAT;
  using DUP = NeedDup_v<NEIGHFLAG,DeviceType>;
  using ScatterAccess = std::conditional_t<
      std::is_same_v<DUP,Kokkos::Experimental::ScatterDuplicated>,
      Kokkos::Experimental::ScatterNonAtomic,
      Kokkos::Experimental::ScatterAtomic>;

  typedef typename KKDevice<DeviceType>::value KKDeviceType;
  template<typename DataType, typename Layout>
  using DupScatterView = KKScatterView<DataType, Layout, KKDeviceType, KKScatterSum, KKScatterDuplicated>;
  template<typename DataType, typename Layout>
  using NonDupScatterView = KKScatterView<DataType, Layout, KKDeviceType, KKScatterSum, KKScatterNonDuplicated>;

  typename AT::t_kkfloat_1d_3_lr_randomread x;
  typename AT::t_int_1d_randomread type;
  typename AT::t_kkfloat_2d d_cutsq;
  typename AT::t_kkfloat_1d d_cn;
  typename AT::t_kkfloat_1d d_dc6;
  typename AT::t_kkfloat_1d d_r2r4;
  typename AT::t_kkfloat_2d d_r0ab;
  Kokkos::View<KK_FLOAT*****, LMPDeviceLayout, DeviceType> d_c6ab;
  typename AT::t_int_1d d_mxci;

  typename AT::t_int_1d d_numneigh;
  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_int_1d d_ilist;

  DupScatterView<KK_ACC_FLOAT*[3], typename DAT::t_kkacc_1d_3::array_layout> dup_f;
  NonDupScatterView<KK_ACC_FLOAT*[3], typename DAT::t_kkacc_1d_3::array_layout> ndup_f;
  DupScatterView<KK_ACC_FLOAT*, typename DAT::t_kkacc_1d::array_layout> dup_eatom;
  NonDupScatterView<KK_ACC_FLOAT*, typename DAT::t_kkacc_1d::array_layout> ndup_eatom;
  DupScatterView<KK_ACC_FLOAT*[6], typename DAT::t_kkacc_1d_6::array_layout> dup_vatom;
  NonDupScatterView<KK_ACC_FLOAT*[6], typename DAT::t_kkacc_1d_6::array_layout> ndup_vatom;
  DupScatterView<KK_FLOAT*, typename DAT::t_kkfloat_1d::array_layout> dup_dc6;
  NonDupScatterView<KK_FLOAT*, typename DAT::t_kkfloat_1d::array_layout> ndup_dc6;

  KK_FLOAT special_lj[4];
  int nlocal;
  int eflag;
  int vflag_either;
  int eflag_global;
  int eflag_atom;
  int vflag_global;
  int vflag_atom;

  int dampingCode;
  KK_FLOAT s6, s8, rs6, rs8, a1, a2, alpha;

  PairDispersionD3KernelA(
      const typename AT::t_kkfloat_1d_3_lr_randomread &x_in,
      const typename AT::t_int_1d_randomread &type_in,
      const typename AT::t_kkfloat_2d &d_cutsq_in,
      const typename AT::t_kkfloat_1d &d_cn_in,
      const typename AT::t_kkfloat_1d &d_dc6_in,
      const typename AT::t_kkfloat_1d &d_r2r4_in,
      const typename AT::t_kkfloat_2d &d_r0ab_in,
      const Kokkos::View<KK_FLOAT*****, LMPDeviceLayout, DeviceType> &d_c6ab_in,
      const typename AT::t_int_1d &d_mxci_in,
      const typename AT::t_int_1d &d_numneigh_in,
      const typename AT::t_neighbors_2d &d_neighbors_in,
      const typename AT::t_int_1d &d_ilist_in,
      const DupScatterView<KK_ACC_FLOAT*[3], typename DAT::t_kkacc_1d_3::array_layout> &dup_f_in,
      const NonDupScatterView<KK_ACC_FLOAT*[3], typename DAT::t_kkacc_1d_3::array_layout> &ndup_f_in,
      const DupScatterView<KK_ACC_FLOAT*, typename DAT::t_kkacc_1d::array_layout> &dup_eatom_in,
      const NonDupScatterView<KK_ACC_FLOAT*, typename DAT::t_kkacc_1d::array_layout> &ndup_eatom_in,
      const DupScatterView<KK_ACC_FLOAT*[6], typename DAT::t_kkacc_1d_6::array_layout> &dup_vatom_in,
      const NonDupScatterView<KK_ACC_FLOAT*[6], typename DAT::t_kkacc_1d_6::array_layout> &ndup_vatom_in,
      const DupScatterView<KK_FLOAT*, typename DAT::t_kkfloat_1d::array_layout> &dup_dc6_in,
      const NonDupScatterView<KK_FLOAT*, typename DAT::t_kkfloat_1d::array_layout> &ndup_dc6_in,
      const KK_FLOAT *special_lj_in,
      int nlocal_in, int eflag_in, int vflag_either_in,
      int eflag_global_in, int eflag_atom_in, int vflag_global_in, int vflag_atom_in,
      int dampingCode_in, KK_FLOAT s6_in, KK_FLOAT s8_in,
      KK_FLOAT rs6_in, KK_FLOAT rs8_in, KK_FLOAT a1_in, KK_FLOAT a2_in,
      KK_FLOAT alpha_in)
    : x(x_in), type(type_in), d_cutsq(d_cutsq_in), d_cn(d_cn_in), d_dc6(d_dc6_in),
      d_r2r4(d_r2r4_in), d_r0ab(d_r0ab_in),
      d_c6ab(d_c6ab_in), d_mxci(d_mxci_in),
      d_numneigh(d_numneigh_in), d_neighbors(d_neighbors_in), d_ilist(d_ilist_in),
      dup_f(dup_f_in), ndup_f(ndup_f_in),
      dup_eatom(dup_eatom_in), ndup_eatom(ndup_eatom_in),
      dup_vatom(dup_vatom_in), ndup_vatom(ndup_vatom_in),
      dup_dc6(dup_dc6_in), ndup_dc6(ndup_dc6_in),
      nlocal(nlocal_in), eflag(eflag_in), vflag_either(vflag_either_in),
      eflag_global(eflag_global_in), eflag_atom(eflag_atom_in),
      vflag_global(vflag_global_in), vflag_atom(vflag_atom_in),
      dampingCode(dampingCode_in), s6(s6_in), s8(s8_in),
      rs6(rs6_in), rs8(rs8_in), a1(a1_in), a2(a2_in),
      alpha(alpha_in)
  {
    special_lj[0] = special_lj_in[0];
    special_lj[1] = special_lj_in[1];
    special_lj[2] = special_lj_in[2];
    special_lj[3] = special_lj_in[3];
  }

  // Extract special bond mask
  // NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  static int sbmask(const int &j) {
    return j >> SBBITS & 3;
  }

  /* ----------------------------------------------------------------------
   Get derivative of C6 on device
  ------------------------------------------------------------------------- */

  // NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void get_dC6_kokkos(KK_FLOAT &c6, KK_FLOAT &dc6i, KK_FLOAT &dc6j,
      const int &iat, const int &jat, const KK_FLOAT &cni, const KK_FLOAT &cnj) const
  {
    KK_FLOAT c6_ref, cni_ref, cnj_ref;
    KK_FLOAT c6mem, r_save, r;
    KK_FLOAT expterm, term;
    KK_FLOAT num, den, d_num_i, d_num_j, d_den_i, d_den_j;

    c6mem = -1.0e20;
    r_save = 1.0e20;
    num = 0.0;
    den = 0.0;
    d_num_i = 0.0;
    d_num_j = 0.0;
    d_den_i = 0.0;
    d_den_j = 0.0;

    int maxci = d_mxci(iat);
    int maxcj = d_mxci(jat);
    for (int ci = 0; ci <= maxci; ci++) {
      for (int cj = 0; cj <= maxcj; cj++) {
        c6_ref = d_c6ab(iat, jat, ci, cj, 0);
        c6_ref *= AUTOEV * AUTOANG6;

        if (c6_ref > 0) {
          cni_ref = d_c6ab(iat, jat, ci, cj, 1);
          cnj_ref = d_c6ab(iat, jat, ci, cj, 2);

          r = (cni - cni_ref) * (cni - cni_ref) + (cnj - cnj_ref) * (cnj - cnj_ref);

          if (r < r_save) {
            r_save = r;
            c6mem = c6_ref;
          }

          expterm = Kokkos::exp(static_cast<KK_FLOAT>(K3) * r);

          num += c6_ref * expterm;
          den += expterm;

          expterm = expterm * static_cast<KK_FLOAT>(2.0 * K3);

          term = expterm * (cni - cni_ref);
          d_num_i += c6_ref * term;
          d_den_i += term;

          term = expterm * (cnj - cnj_ref);
          d_num_j += c6_ref * term;
          d_den_j += term;
        }
      }
    }

    // the reference threshold of 1.0e-99 underflows to zero in a single
    // precision build, so use the smallest normalized value of KK_FLOAT

    if (den > std::numeric_limits<KK_FLOAT>::min()) {
      c6 = num / den;
      dc6i = ((d_num_i * den) - (d_den_i * num)) / (den * den);
      dc6j = ((d_num_j * den) - (d_den_j * num)) / (den * den);
    } else {
      c6 = c6mem;
      dc6i = 0;
      dc6j = 0;
    }
  }

  // NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void ev_tally(value_type &ev, const int &i, const int &j,
        const KK_FLOAT &epair, const KK_FLOAT &fpair, const KK_FLOAT &delx,
                  const KK_FLOAT &dely, const KK_FLOAT &delz) const
  {
    const int EFLAG = eflag;
    const int VFLAG = vflag_either;

    if (EFLAG) {
      const KK_FLOAT epairhalf = 0.5 * epair;
      if (eflag_global) {
        if (NEIGHFLAG!=FULL && (NEWTON_PAIR || j < nlocal)) {
          ev.evdwl += epair;
        } else {
          ev.evdwl += epairhalf;
        }
      }
      if (eflag_atom) {
        auto v_eatom = ScatterViewHelper<DUP,decltype(dup_eatom),decltype(ndup_eatom)>::get(dup_eatom,ndup_eatom);
        auto a_eatom = v_eatom.template access<ScatterAccess>();
        if (NEIGHFLAG!=FULL) {
          if (NEWTON_PAIR || i < nlocal) a_eatom[i] += epairhalf;
          if (NEWTON_PAIR || j < nlocal) a_eatom[j] += epairhalf;
        } else {
          a_eatom[i] += epairhalf;
        }
      }
    }

    if (VFLAG) {
      const KK_FLOAT v0 = delx*delx*fpair;
      const KK_FLOAT v1 = dely*dely*fpair;
      const KK_FLOAT v2 = delz*delz*fpair;
      const KK_FLOAT v3 = delx*dely*fpair;
      const KK_FLOAT v4 = delx*delz*fpair;
      const KK_FLOAT v5 = dely*delz*fpair;

      if (vflag_global) {
        if (NEIGHFLAG!=FULL) {
          if (NEWTON_PAIR || i < nlocal) {
            ev.v[0] += 0.5*v0;
            ev.v[1] += 0.5*v1;
            ev.v[2] += 0.5*v2;
            ev.v[3] += 0.5*v3;
            ev.v[4] += 0.5*v4;
            ev.v[5] += 0.5*v5;
          }
          if (NEWTON_PAIR || j < nlocal) {
            ev.v[0] += 0.5*v0;
            ev.v[1] += 0.5*v1;
            ev.v[2] += 0.5*v2;
            ev.v[3] += 0.5*v3;
            ev.v[4] += 0.5*v4;
            ev.v[5] += 0.5*v5;
          }
        } else {
          ev.v[0] += 0.5*v0;
          ev.v[1] += 0.5*v1;
          ev.v[2] += 0.5*v2;
          ev.v[3] += 0.5*v3;
          ev.v[4] += 0.5*v4;
          ev.v[5] += 0.5*v5;
        }
      }

      if (vflag_atom) {
        auto v_vatom = ScatterViewHelper<DUP,decltype(dup_vatom),decltype(ndup_vatom)>::get(dup_vatom,ndup_vatom);
        auto a_vatom = v_vatom.template access<ScatterAccess>();
        if (NEIGHFLAG!=FULL) {
          if (NEWTON_PAIR || i < nlocal) {
            a_vatom(i,0) += 0.5*v0;
            a_vatom(i,1) += 0.5*v1;
            a_vatom(i,2) += 0.5*v2;
            a_vatom(i,3) += 0.5*v3;
            a_vatom(i,4) += 0.5*v4;
            a_vatom(i,5) += 0.5*v5;
          }
          if (NEWTON_PAIR || j < nlocal) {
            a_vatom(j,0) += 0.5*v0;
            a_vatom(j,1) += 0.5*v1;
            a_vatom(j,2) += 0.5*v2;
            a_vatom(j,3) += 0.5*v3;
            a_vatom(j,4) += 0.5*v4;
            a_vatom(j,5) += 0.5*v5;
          }
        } else {
          a_vatom(i,0) += 0.5*v0;
          a_vatom(i,1) += 0.5*v1;
          a_vatom(i,2) += 0.5*v2;
          a_vatom(i,3) += 0.5*v3;
          a_vatom(i,4) += 0.5*v4;
          a_vatom(i,5) += 0.5*v5;
        }
      }
    }
  }

  // NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(const int &ii, value_type &ev) const {
    auto v_f = ScatterViewHelper<DUP,decltype(dup_f),decltype(ndup_f)>::get(dup_f,ndup_f);
    auto a_f = v_f.template access<ScatterAccess>();

    auto v_dc6 = ScatterViewHelper<DUP,decltype(dup_dc6),decltype(ndup_dc6)>::get(dup_dc6,ndup_dc6);
    auto a_dc6 = v_dc6.template access<ScatterAccess>();

    const int i = d_ilist(ii);
    const int itype = type(i);
    int jnum = d_numneigh(i);

    const KK_FLOAT xtmp = x(i, 0);
    const KK_FLOAT ytmp = x(i, 1);
    const KK_FLOAT ztmp = x(i, 2);

    KK_FLOAT fxtmp = 0.0;
    KK_FLOAT fytmp = 0.0;
    KK_FLOAT fztmp = 0.0;
    KK_FLOAT dc6_i = 0.0;

    for (int jj = 0; jj < jnum; jj++) {
      int j = d_neighbors(i, jj);
      KK_FLOAT factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      const KK_FLOAT delx = xtmp - x(j, 0);
      const KK_FLOAT dely = ytmp - x(j, 1);
      const KK_FLOAT delz = ztmp - x(j, 2);

      const KK_FLOAT rsq = delx * delx + dely * dely + delz * delz;

      const int jtype = type(j);

      if (rsq < d_cutsq(itype, jtype)) {

        const KK_FLOAT r2inv = 1.0f / rsq;
        const KK_FLOAT r6inv = r2inv * r2inv * r2inv;
        const KK_FLOAT r8inv = r6inv * r2inv;
        const KK_FLOAT r10inv = r8inv * r2inv;

        KK_FLOAT c6 = 0.0;
        KK_FLOAT dc6i = 0.0;
        KK_FLOAT dc6j = 0.0;
        get_dC6_kokkos(c6, dc6i, dc6j, itype, jtype, d_cn(i), d_cn(j));

        const KK_FLOAT C6 = c6;
        const KK_FLOAT C8 = 3.0 * C6 * d_r2r4(itype) * d_r2r4(jtype) * AUTOANG * AUTOANG;

        const KK_FLOAT alpha6 = alpha;
        const KK_FLOAT alpha8 = alpha + 2;

        KK_FLOAT t6, t8, damp6, damp8, e6, e8;
        KK_FLOAT tmp6, tmp8, fpair1, fpair2, fpair;
        KK_FLOAT evdwl = 0.0;
        t6 = t8 = e6 = e8 = fpair = fpair1 = fpair2 = 0.0;

        switch (dampingCode) {

          case 1: {    // original

            const KK_FLOAT ip6 = rs6 * d_r0ab(itype, jtype);
            const KK_FLOAT ip8 = rs8 * d_r0ab(itype, jtype);

            const KK_FLOAT half_alpha6 = 0.5 * alpha6;
            const KK_FLOAT half_alpha8 = 0.5 * alpha8;

            t6 = MathSpecialKokkos::powauto(ip6, alpha6) * MathSpecialKokkos::powauto(rsq, -half_alpha6);
            t8 = MathSpecialKokkos::powauto(ip8, alpha8) * MathSpecialKokkos::powauto(rsq, -half_alpha8);

            damp6 = 1.0f / (1.0f + 6.0f * t6);
            damp8 = 1.0f / (1.0f + 6.0f * t8);

            e6 = C6 * damp6 * r6inv;
            e8 = C8 * damp8 * r8inv;

            tmp6 = 6 * s6 * C6 * r8inv * damp6;
            tmp8 = 8 * s8 * C8 * r10inv * damp8;

            fpair1 = -tmp6 - tmp8;
            fpair2 = tmp6 * alpha6 * t6 * damp6 + (3.0f / 4) * tmp8 * alpha8 * t8 * damp8;

            fpair = fpair1 + fpair2;
            fpair *= factor_lj;

          } break;

          case 2: {    // zerom

            const KK_FLOAT r0 = d_r0ab(itype, jtype);
            const KK_FLOAT r = Kokkos::sqrt(rsq);

            t6 = MathSpecialKokkos::powauto((r / (rs6 * r0)) + rs8 * r0, -alpha6);
            damp6 = 1.0f / (1.0f + 6.0f * t6);
            t8 = MathSpecialKokkos::powauto((r / r0) + rs8 * r0, -alpha8);
            damp8 = 1.0f / (1.0f + 6.0f * t8);

            e6 = C6 * damp6 * r6inv;
            e8 = C8 * damp8 * r8inv;

            tmp6 = 6 * s6 * C6 * r8inv * damp6;
            tmp8 = 8 * s8 * C8 * r10inv * damp8;

            fpair1 = -tmp6 - tmp8;

            const KK_FLOAT fp26 = tmp6 * alpha6 * t6 * damp6 * r / (r + rs6 * rs8 * r0 * r0);
            const KK_FLOAT fp28 = tmp8 * alpha8 * t8 * damp8 * r / (r + rs8 * r0 * r0);

            fpair2 = fp26 + (3.0f / 4) * fp28;

            fpair = fpair1 + fpair2;
            fpair *= factor_lj;
          } break;

          case 3:      // bj
          case 4: {    // bjm, same functional form as bj, different parameters

            const KK_FLOAT r0 = Kokkos::sqrt(C8 / C6);

            const KK_FLOAT r4 = rsq * rsq;
            KK_FLOAT r6 = rsq * rsq * rsq;
            KK_FLOAT r8 = r6 * rsq;

            const KK_FLOAT d = a1 * r0 + a2;
            const KK_FLOAT d2 = d * d;
            const KK_FLOAT d4 = d2 * d2;

            t6 = r6 + MathSpecialKokkos::cube(d2);
            t8 = r8 + MathSpecialKokkos::square(d4);

            e6 = C6 / t6;
            e8 = C8 / t8;

            tmp6 = 6.0 * s6 * C6 * r4 / (t6 * t6);
            tmp8 = 8.0 * s8 * C8 * r6 / (t8 * t8);

            fpair = -(tmp6 + tmp8);
            fpair *= factor_lj;
          } break;

          // no default case: dampingCode is validated in init_style()
        }

        if (EVFLAG) evdwl = -(s6 * e6 + s8 * e8) * factor_lj;

        const KK_FLOAT rest = (s6 * e6 + s8 * e8) / C6;

        dc6_i += rest * dc6i;
        if (NEIGHFLAG != FULL && (NEWTON_PAIR || j < nlocal)) { a_dc6(j) += rest * dc6j; }

        fxtmp += delx * fpair;
        fytmp += dely * fpair;
        fztmp += delz * fpair;

        if (NEIGHFLAG != FULL && (NEWTON_PAIR || j < nlocal)) {
          a_f(j, 0) -= delx * fpair;
          a_f(j, 1) -= dely * fpair;
          a_f(j, 2) -= delz * fpair;
        }

        if (EVFLAG) {
          ev_tally(ev, i, j, evdwl, fpair, delx, dely, delz);
        }
      }
    }

    a_f(i, 0) += fxtmp;
    a_f(i, 1) += fytmp;
    a_f(i, 2) += fztmp;
    a_dc6(i) += dc6_i;
  }

  // NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(const int &ii) const {
    value_type ev;
    operator()(ii, ev);
  }
};

/* ----------------------------------------------------------------------
   Kernel B: compute force contribution from dC6
------------------------------------------------------------------------- */

template <class DeviceType, int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
struct PairDispersionD3KernelB {
  typedef ArrayTypes<DeviceType> AT;
  using value_type = EV_FLOAT;
  using DUP = NeedDup_v<NEIGHFLAG,DeviceType>;
  using ScatterAccess = std::conditional_t<
      std::is_same_v<DUP,Kokkos::Experimental::ScatterDuplicated>,
      Kokkos::Experimental::ScatterNonAtomic,
      Kokkos::Experimental::ScatterAtomic>;

  typedef typename KKDevice<DeviceType>::value KKDeviceType;
  template<typename DataType, typename Layout>
  using DupScatterView = KKScatterView<DataType, Layout, KKDeviceType, KKScatterSum, KKScatterDuplicated>;
  template<typename DataType, typename Layout>
  using NonDupScatterView = KKScatterView<DataType, Layout, KKDeviceType, KKScatterSum, KKScatterNonDuplicated>;

  typename AT::t_kkfloat_1d_3_lr_randomread x;
  typename AT::t_int_1d_randomread type;
  typename AT::t_kkfloat_2d d_cutsq;
  typename AT::t_kkfloat_1d d_dc6;
  typename AT::t_kkfloat_1d d_rcov;
  typename AT::t_int_1d d_numneigh;
  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_int_1d d_ilist;

  DupScatterView<KK_ACC_FLOAT*[3], typename DAT::t_kkacc_1d_3::array_layout> dup_f;
  NonDupScatterView<KK_ACC_FLOAT*[3], typename DAT::t_kkacc_1d_3::array_layout> ndup_f;
  DupScatterView<KK_ACC_FLOAT*, typename DAT::t_kkacc_1d::array_layout> dup_eatom;
  NonDupScatterView<KK_ACC_FLOAT*, typename DAT::t_kkacc_1d::array_layout> ndup_eatom;
  DupScatterView<KK_ACC_FLOAT*[6], typename DAT::t_kkacc_1d_6::array_layout> dup_vatom;
  NonDupScatterView<KK_ACC_FLOAT*[6], typename DAT::t_kkacc_1d_6::array_layout> ndup_vatom;

  KK_FLOAT special_lj[4];
  int nlocal;
  int eflag;
  int vflag_either;
  int eflag_global;
  int eflag_atom;
  int vflag_global;
  int vflag_atom;
  KK_FLOAT cn_thr;

  PairDispersionD3KernelB(
      const typename AT::t_kkfloat_1d_3_lr_randomread &x_in,
      const typename AT::t_int_1d_randomread &type_in,
      const typename AT::t_kkfloat_2d &d_cutsq_in,
      const typename AT::t_kkfloat_1d &d_dc6_in,
      const typename AT::t_kkfloat_1d &d_rcov_in,
      const typename AT::t_int_1d &d_numneigh_in,
      const typename AT::t_neighbors_2d &d_neighbors_in,
      const typename AT::t_int_1d &d_ilist_in,
      const DupScatterView<KK_ACC_FLOAT*[3], typename DAT::t_kkacc_1d_3::array_layout> &dup_f_in,
      const NonDupScatterView<KK_ACC_FLOAT*[3], typename DAT::t_kkacc_1d_3::array_layout> &ndup_f_in,
      const DupScatterView<KK_ACC_FLOAT*, typename DAT::t_kkacc_1d::array_layout> &dup_eatom_in,
      const NonDupScatterView<KK_ACC_FLOAT*, typename DAT::t_kkacc_1d::array_layout> &ndup_eatom_in,
      const DupScatterView<KK_ACC_FLOAT*[6], typename DAT::t_kkacc_1d_6::array_layout> &dup_vatom_in,
      const NonDupScatterView<KK_ACC_FLOAT*[6], typename DAT::t_kkacc_1d_6::array_layout> &ndup_vatom_in,
      const KK_FLOAT *special_lj_in,
      int nlocal_in, int eflag_in, int vflag_either_in,
      int eflag_global_in, int eflag_atom_in, int vflag_global_in, int vflag_atom_in,
      KK_FLOAT cn_thr_in)
    : x(x_in), type(type_in), d_cutsq(d_cutsq_in), d_dc6(d_dc6_in),
      d_rcov(d_rcov_in), d_numneigh(d_numneigh_in), d_neighbors(d_neighbors_in),
      d_ilist(d_ilist_in), dup_f(dup_f_in), ndup_f(ndup_f_in),
      dup_eatom(dup_eatom_in), ndup_eatom(ndup_eatom_in),
      dup_vatom(dup_vatom_in), ndup_vatom(ndup_vatom_in),
      nlocal(nlocal_in), eflag(eflag_in), vflag_either(vflag_either_in),
      eflag_global(eflag_global_in), eflag_atom(eflag_atom_in),
      vflag_global(vflag_global_in), vflag_atom(vflag_atom_in),
      cn_thr(cn_thr_in)
  {
    special_lj[0] = special_lj_in[0];
    special_lj[1] = special_lj_in[1];
    special_lj[2] = special_lj_in[2];
    special_lj[3] = special_lj_in[3];
  }

  // Extract special bond mask
  // NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  static int sbmask(const int &j) {
    return j >> SBBITS & 3;
  }


  // NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void ev_tally(value_type &ev, const int &i, const int &j,
        const KK_FLOAT &epair, const KK_FLOAT &fpair, const KK_FLOAT &delx,
                  const KK_FLOAT &dely, const KK_FLOAT &delz) const
  {
    const int EFLAG = eflag;
    const int VFLAG = vflag_either;

    if (EFLAG) {
      const KK_FLOAT epairhalf = 0.5 * epair;
      if (eflag_global) {
        if (NEIGHFLAG!=FULL && (NEWTON_PAIR || j < nlocal)) {
          ev.evdwl += epair;
        } else {
          ev.evdwl += epairhalf;
        }
      }
      if (eflag_atom) {
        auto v_eatom = ScatterViewHelper<DUP,decltype(dup_eatom),decltype(ndup_eatom)>::get(dup_eatom,ndup_eatom);
        auto a_eatom = v_eatom.template access<ScatterAccess>();
        if (NEIGHFLAG!=FULL) {
          if (NEWTON_PAIR || i < nlocal) a_eatom[i] += epairhalf;
          if (NEWTON_PAIR || j < nlocal) a_eatom[j] += epairhalf;
        } else {
          a_eatom[i] += epairhalf;
        }
      }
    }

    if (VFLAG) {
      const KK_FLOAT v0 = delx*delx*fpair;
      const KK_FLOAT v1 = dely*dely*fpair;
      const KK_FLOAT v2 = delz*delz*fpair;
      const KK_FLOAT v3 = delx*dely*fpair;
      const KK_FLOAT v4 = delx*delz*fpair;
      const KK_FLOAT v5 = dely*delz*fpair;

      if (vflag_global) {
        if (NEIGHFLAG!=FULL) {
          if (NEWTON_PAIR || i < nlocal) {
            ev.v[0] += 0.5*v0;
            ev.v[1] += 0.5*v1;
            ev.v[2] += 0.5*v2;
            ev.v[3] += 0.5*v3;
            ev.v[4] += 0.5*v4;
            ev.v[5] += 0.5*v5;
          }
          if (NEWTON_PAIR || j < nlocal) {
            ev.v[0] += 0.5*v0;
            ev.v[1] += 0.5*v1;
            ev.v[2] += 0.5*v2;
            ev.v[3] += 0.5*v3;
            ev.v[4] += 0.5*v4;
            ev.v[5] += 0.5*v5;
          }
        } else {
          ev.v[0] += 0.5*v0;
          ev.v[1] += 0.5*v1;
          ev.v[2] += 0.5*v2;
          ev.v[3] += 0.5*v3;
          ev.v[4] += 0.5*v4;
          ev.v[5] += 0.5*v5;
        }
      }

      if (vflag_atom) {
        auto v_vatom = ScatterViewHelper<DUP,decltype(dup_vatom),decltype(ndup_vatom)>::get(dup_vatom,ndup_vatom);
        auto a_vatom = v_vatom.template access<ScatterAccess>();
        if (NEIGHFLAG!=FULL) {
          if (NEWTON_PAIR || i < nlocal) {
            a_vatom(i,0) += 0.5*v0;
            a_vatom(i,1) += 0.5*v1;
            a_vatom(i,2) += 0.5*v2;
            a_vatom(i,3) += 0.5*v3;
            a_vatom(i,4) += 0.5*v4;
            a_vatom(i,5) += 0.5*v5;
          }
          if (NEWTON_PAIR || j < nlocal) {
            a_vatom(j,0) += 0.5*v0;
            a_vatom(j,1) += 0.5*v1;
            a_vatom(j,2) += 0.5*v2;
            a_vatom(j,3) += 0.5*v3;
            a_vatom(j,4) += 0.5*v4;
            a_vatom(j,5) += 0.5*v5;
          }
        } else {
          a_vatom(i,0) += 0.5*v0;
          a_vatom(i,1) += 0.5*v1;
          a_vatom(i,2) += 0.5*v2;
          a_vatom(i,3) += 0.5*v3;
          a_vatom(i,4) += 0.5*v4;
          a_vatom(i,5) += 0.5*v5;
        }
      }
    }
  }

  // NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(const int &ii, value_type &ev) const {
    auto v_f = ScatterViewHelper<DUP,decltype(dup_f),decltype(ndup_f)>::get(dup_f,ndup_f);
    auto a_f = v_f.template access<ScatterAccess>();

    const int i = d_ilist(ii);
    const int itype = type(i);
    int jnum = d_numneigh(i);

    const KK_FLOAT xtmp = x(i, 0);
    const KK_FLOAT ytmp = x(i, 1);
    const KK_FLOAT ztmp = x(i, 2);

    KK_FLOAT fxtmp = 0.0;
    KK_FLOAT fytmp = 0.0;
    KK_FLOAT fztmp = 0.0;

    for (int jj = 0; jj < jnum; jj++) {
      int j = d_neighbors(i, jj);
      KK_FLOAT factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      const KK_FLOAT delx = xtmp - x(j, 0);
      const KK_FLOAT dely = ytmp - x(j, 1);
      const KK_FLOAT delz = ztmp - x(j, 2);
      const KK_FLOAT rsq = delx * delx + dely * dely + delz * delz;

      const int jtype = type(j);

      if (rsq < d_cutsq(itype, jtype)) {

        const KK_FLOAT r = Kokkos::sqrt(rsq);
        KK_FLOAT dcn;

        if (rsq < cn_thr) {
          const KK_FLOAT rcovij = (d_rcov(itype) + d_rcov(jtype)) * AUTOANG;
          const KK_FLOAT expterm = Kokkos::exp(-K1 * (rcovij / r - 1.0));
          dcn = -K1 * rcovij * expterm / (rsq * (expterm + 1.0) * (expterm + 1.0));
        } else {
          dcn = 0.0;
        }

        KK_FLOAT fpair = dcn * (d_dc6(i) + d_dc6(j)) / r;
        fpair *= factor_lj;

        fxtmp += delx * fpair;
        fytmp += dely * fpair;
        fztmp += delz * fpair;

        if (NEIGHFLAG != FULL && (NEWTON_PAIR || j < nlocal)) {
          a_f(j, 0) -= delx * fpair;
          a_f(j, 1) -= dely * fpair;
          a_f(j, 2) -= delz * fpair;
        }

        if (EVFLAG) {
          const KK_FLOAT epair = 0.0;
          ev_tally(ev, i, j, epair, fpair, delx, dely, delz);
        }
      }
    }

    a_f(i, 0) += fxtmp;
    a_f(i, 1) += fytmp;
    a_f(i, 2) += fztmp;
  }

  // NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(const int &ii) const {
    value_type ev;
    operator()(ii, ev);
  }
};

template<class DeviceType>
class PairDispersionD3Kokkos : public PairDispersionD3, public KokkosBase {
 public:
  enum {EnabledNeighFlags=FULL|HALFTHREAD|HALF};
  enum {COUL_FLAG=0};
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef ArrayTypes<LMPHostType> HAT;
  typedef EV_FLOAT value_type;
  PairDispersionD3Kokkos(class LAMMPS *);
  ~PairDispersionD3Kokkos() override;

  void calc_coordination_number() override;
  void compute(int, int) override;
  double init_one(int, int) override;
  void init_style() override;
  void allocate() override;
  void coeff(int, char **) override;

  int pack_forward_comm_kokkos(int, DAT::tdual_int_1d, DAT::tdual_double_1d&,
                       int, int *) override;
  void unpack_forward_comm_kokkos(int, int, DAT::tdual_double_1d&) override;
  int pack_reverse_comm_kokkos(int, int, DAT::tdual_double_1d&) override;
  void unpack_reverse_comm_kokkos(int, DAT::tdual_int_1d, DAT::tdual_double_1d&) override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;

 protected:
  typename AT::t_kkfloat_1d_3_lr_randomread x; // atom positions
  typename AT::t_kkacc_1d_3 f;                 // atom forces
  typename AT::t_int_1d_randomread type;       // atom types

  DAT::ttransform_kkacc_1d k_eatom;    // per-atom energy (dual view)
  DAT::ttransform_kkacc_1d_6 k_vatom;  // per-atom virial (dual view)
  typename AT::t_kkacc_1d d_eatom;     // device view of per-atom energy
  typename AT::t_kkacc_1d_6 d_vatom;   // device view of per-atom virial

  KK_FLOAT special_lj[4]; // special-bond scaling
  int inum;               // number of neighbor list atoms
  bool need_dup;          // whether duplicated scatter is required

  typedef typename KKDevice<DeviceType>::value KKDeviceType;
  template<typename DataType, typename Layout>
  using DupScatterView = KKScatterView<DataType, Layout, KKDeviceType, KKScatterSum, KKScatterDuplicated>;

  template<typename DataType, typename Layout>
  using NonDupScatterView = KKScatterView<DataType, Layout, KKDeviceType, KKScatterSum, KKScatterNonDuplicated>;

  DupScatterView<KK_ACC_FLOAT*[3], typename DAT::t_kkacc_1d_3::array_layout> dup_f; // duplicated force
  DupScatterView<KK_ACC_FLOAT*, typename DAT::t_kkacc_1d::array_layout> dup_eatom; // duplicated energy
  DupScatterView<KK_ACC_FLOAT*[6], typename DAT::t_kkacc_1d_6::array_layout> dup_vatom; // duplicated virial
  DupScatterView<KK_FLOAT*, typename DAT::t_kkfloat_1d::array_layout> dup_dc6;  // duplicated dC6
  NonDupScatterView<KK_ACC_FLOAT*[3], typename DAT::t_kkacc_1d_3::array_layout> ndup_f; // non-dup force
  NonDupScatterView<KK_ACC_FLOAT*, typename DAT::t_kkacc_1d::array_layout> ndup_eatom; // non-dup energy
  NonDupScatterView<KK_ACC_FLOAT*[6], typename DAT::t_kkacc_1d_6::array_layout> ndup_vatom; // non-dup virial
  NonDupScatterView<KK_FLOAT*, typename DAT::t_kkfloat_1d::array_layout> ndup_dc6; // non-dup dC6


  DAT::tdual_kkfloat_1d k_cn;   // coordination numbers (dual view)
  DAT::tdual_kkfloat_1d k_dc6;  // dC6 values (dual view)
  typename AT::t_kkfloat_1d d_cn;   // device CN
  typename AT::t_kkfloat_1d d_dc6;  // device dC6
  HAT::t_kkfloat_1d h_cn;           // host CN
  HAT::t_kkfloat_1d h_dc6;          // host dC6

  DAT::tdual_kkfloat_1d k_r2r4;   // r2r4 table (dual view)
  DAT::tdual_kkfloat_1d k_rcov;   // covalent radii (dual view)
  DAT::tdual_int_1d k_mxci;       // max C6 grid index (dual view)
  DAT::tdual_kkfloat_2d k_r0ab;   // R0 table (dual view)
  Kokkos::DualView<KK_FLOAT*****, LMPDeviceLayout, DeviceType> k_c6ab; // C6 table (dual view)

  typename AT::t_kkfloat_1d d_r2r4;   // device r2r4
  typename AT::t_kkfloat_1d d_rcov;   // device covalent radii
  typename AT::t_int_1d d_mxci;       // device max C6 grid index
  typename AT::t_kkfloat_2d d_r0ab;   // device R0 table
  Kokkos::View<KK_FLOAT*****, LMPDeviceLayout, DeviceType> d_c6ab; // device C6 table

  DAT::ttransform_kkfloat_2d k_cutsq;    // cutoff^2 table (host double, device KK_FLOAT)
  typename AT::t_kkfloat_2d d_cutsq;     // device cutoff^2 table

  typename AT::t_neighbors_2d d_neighbors; // neighbor list
  typename AT::t_int_1d d_ilist;           // neighbor list indices
  typename AT::t_int_1d d_numneigh;        // neighbor counts

  typename AT::t_int_1d d_sendlist, d_recvlist; // comm lists
  typename AT::t_double_1d_um v_buf;            // comm buffer

  int neighflag, newton_pair; // neighbor/newton settings
  int nlocal, nall, eflag, vflag; // local/total counts and flags

  friend void pair_virial_fdotr_compute<PairDispersionD3Kokkos<DeviceType>>(PairDispersionD3Kokkos<DeviceType>*);

  // To make the compute() function cleaner:
  template<int NEIGHFLAG>
  void dispatch_coordination_kernel();

  template<int NEIGHFLAG>
  void dispatch_kernel_A(EV_FLOAT &ev);

  template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
  void launch_kernel_A(EV_FLOAT &ev);

  template<int NEIGHFLAG>
  void dispatch_kernel_B(EV_FLOAT &ev);

  template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
  void launch_kernel_B(EV_FLOAT &ev);

};

}

#endif
#endif
