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

#include "meam_kokkos.h"
#include "math_special_kokkos.h"

using namespace LAMMPS_NS;
using namespace MathSpecialKokkos;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG>
KOKKOS_INLINE_FUNCTION
void MEAMKokkos<DeviceType>::operator()(TagMEAMDensInit<NEIGHFLAG>, const int &i) const {
  int ii, offsetval;
  ii = d_ilist_half[i];
  offsetval = d_offset[i];
  // compute screening function and derivatives
  this->template getscreen<NEIGHFLAG>(ii, offsetval, x, d_numneigh_half,
            d_numneigh_full, ntype, type, d_map);

  // calculate intermediate density terms to be communicated
  this->template calc_rho1<NEIGHFLAG>(ii, ntype, type, d_map, x, d_numneigh_half, offsetval);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void MEAMKokkos<DeviceType>::operator()(TagMEAMZero, const int &i) const {
  d_rho0[i] = 0.0;
  d_arho2b[i] = 0.0;
  d_arho1(i,0) = d_arho1(i,1) = d_arho1(i,2) = 0.0;
  if (msmeamflag) {
    d_arho2mb[i] = 0.0;
    d_arho1m(i,0) = d_arho1m(i,1) = d_arho1m(i,2) = 0.0;
  }
  for (int j = 0; j < 6; j++) {
    d_arho2(i,j) = 0.0;
    if (msmeamflag)
      d_arho2m(i,j) = 0.0;
  }
  for (int j = 0; j < 10; j++) {
    d_arho3(i,j) = 0.0;
    if (msmeamflag)
      d_arho3m(i,j) = 0.0;
  }
  d_arho3b(i,0) = d_arho3b(i,1) = d_arho3b(i,2) = 0.0;
  if (msmeamflag)
    d_arho3mb(i,0) = d_arho3mb(i,1) = d_arho3mb(i,2) = 0.0;
  d_t_ave(i,0) = d_t_ave(i,1) = d_t_ave(i,2) = 0.0;
  d_tsq_ave(i,0) = d_tsq_ave(i,1) = d_tsq_ave(i,2) = 0.0;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void
MEAMKokkos<DeviceType>::meam_dens_setup(int atom_nmax, int nall, int n_neigh)
{
  // grow local arrays if necessary

  if (atom_nmax > nmax) {
    nmax = atom_nmax;

    //memory->create(rho, nmax, "pair:rho");
    k_rho = DAT::tdual_kkfloat_1d("pair:rho",nmax);
    d_rho = k_rho.template view<DeviceType>();
    h_rho = k_rho.view_host();
    //memory->create(rho0, nmax, "pair:rho0");
    k_rho0 = DAT::tdual_kkfloat_1d("pair:rho0",nmax);
    d_rho0 = k_rho0.template view<DeviceType>();
    h_rho0 = k_rho0.view_host();
    //memory->create(rho1, nmax, "pair:rho1");
    k_rho1 = DAT::tdual_kkfloat_1d("pair:rho1",nmax);
    d_rho1 = k_rho1.template view<DeviceType>();
    h_rho1 = k_rho1.view_host();
    //memory->create(rho2, nmax, "pair:rho2");
    k_rho2 = DAT::tdual_kkfloat_1d("pair:rho2",nmax);
    d_rho2 = k_rho2.template view<DeviceType>();
    h_rho2 = k_rho2.view_host();
    //memory->create(rho3, nmax, "pair:rho3");
    k_rho3 = DAT::tdual_kkfloat_1d("pair:rho3",nmax);
    d_rho3 = k_rho3.template view<DeviceType>();
    h_rho3 = k_rho3.view_host();
    //memory->create(frhop, nmax, "pair:frhop");
    k_frhop = DAT::tdual_kkfloat_1d("pair:frhop",nmax);
    d_frhop = k_frhop.template view<DeviceType>();
    h_frhop = k_frhop.view_host();
    //memory->create(gamma, nmax, "pair:gamma");
    k_gamma = DAT::tdual_kkfloat_1d("pair:gamma",nmax);
    d_gamma = k_gamma.template view<DeviceType>();
    h_gamma = k_gamma.view_host();
    //memory->create(dgamma1, nmax, "pair:dgamma1");
    k_dgamma1 = DAT::tdual_kkfloat_1d("pair:dgamma1",nmax);
    d_dgamma1 = k_dgamma1.template view<DeviceType>();
    h_dgamma1 = k_dgamma1.view_host();
    //memory->create(dgamma2, nmax, "pair:dgamma2");
    k_dgamma2 = DAT::tdual_kkfloat_1d("pair:dgamma2",nmax);
    d_dgamma2 = k_dgamma2.template view<DeviceType>();
    h_dgamma2 = k_dgamma2.view_host();
    //memory->create(dgamma3, nmax, "pair:dgamma3");
    k_dgamma3 = DAT::tdual_kkfloat_1d("pair:dgamma3",nmax);
    d_dgamma3 = k_dgamma3.template view<DeviceType>();
    h_dgamma3 = k_dgamma3.view_host();
    //memory->create(arho2b, nmax, "pair:arho2b");
    k_arho2b = DAT::tdual_kkfloat_1d("pair:arho2b",nmax);
    d_arho2b = k_arho2b.template view<DeviceType>();
    h_arho2b = k_arho2b.view_host();
    //memory->create(arho1, nmax, 3, "pair:arho1");
    k_arho1 = DAT::tdual_kkfloat_2d("pair:arho1",nmax, 3);
    d_arho1 = k_arho1.template view<DeviceType>();
    h_arho1 = k_arho1.view_host();
    //memory->create(arho2, nmax, 6, "pair:arho2");
    k_arho2 = DAT::tdual_kkfloat_2d("pair:arho2",nmax, 6);
    d_arho2 = k_arho2.template view<DeviceType>();
    h_arho2 = k_arho2.view_host();
    //memory->create(arho3, nmax, 10, "pair:arho3");
    k_arho3 = DAT::tdual_kkfloat_2d("pair:arho3",nmax, 10);
    d_arho3 = k_arho3.template view<DeviceType>();
    h_arho3 = k_arho3.view_host();
    //memory->create(arho3b, nmax, 3, "pair:arho3b");
    k_arho3b = DAT::tdual_kkfloat_2d("pair:arho3b",nmax, 3);
    d_arho3b = k_arho3b.template view<DeviceType>();
    h_arho3b = k_arho3b.view_host();
    //memory->create(t_ave, nmax, 3, "pair:t_ave");
    k_t_ave = DAT::tdual_kkfloat_2d("pair:t_ave",nmax, 3);
    d_t_ave = k_t_ave.template view<DeviceType>();
    h_t_ave = k_t_ave.view_host();
    //memory->create(tsq_ave, nmax, 3, "pair:tsq_ave");
    k_tsq_ave = DAT::tdual_kkfloat_2d("pair:tsq_ave",nmax, 3);
    d_tsq_ave = k_tsq_ave.template view<DeviceType>();
    h_tsq_ave = k_tsq_ave.view_host();

    // msmeam
    //memory->create(arho2mb, nmax, "pair:arho2mb");
    k_arho2mb = DAT::tdual_kkfloat_1d("pair:arho2mb",nmax);
    d_arho2mb = k_arho2mb.template view<DeviceType>();
    h_arho2mb = k_arho2mb.view_host();
    //memory->create(arho1m, nmax, 3, "pair:arho1m");
    k_arho1m = DAT::tdual_kkfloat_2d("pair:arho1m", nmax, 3);
    d_arho1m = k_arho1m.template view<DeviceType>();
    h_arho1m = k_arho1m.view_host();
    //memory->create(arho2m, nmax, 6, "pair:arho2m");
    k_arho2m = DAT::tdual_kkfloat_2d("pair:arho2m", nmax, 6);
    d_arho2m = k_arho2m.template view<DeviceType>();
    h_arho2m = k_arho2m.view_host();
    //memory->create(arho3m, nmax, 10, "pair:arho3m");
    k_arho3m = DAT::tdual_kkfloat_2d("pair:arho3m", nmax, 10);
    d_arho3m = k_arho3m.template view<DeviceType>();
    h_arho3m = k_arho3m.view_host();
    //memory->create(arho3mb, nmax, 3, "pair:arho3mb");
    k_arho3mb = DAT::tdual_kkfloat_2d("pair:arho3mb", nmax, 3);
    d_arho3mb = k_arho3mb.template view<DeviceType>();
    h_arho3mb = k_arho3mb.view_host();
  }

  if (n_neigh > maxneigh) {
    maxneigh = n_neigh;
   // memory->create(scrfcn, maxneigh, "pair:scrfcn");
    k_scrfcn = DAT::tdual_kkfloat_1d("pair:scrfcn", maxneigh);
    d_scrfcn = k_scrfcn.template view<DeviceType>();
    h_scrfcn = k_scrfcn.view_host();
    //memory->create(dscrfcn, maxneigh, "pair:dscrfcn");
    k_dscrfcn = DAT::tdual_kkfloat_1d("pair:dscrfcn", maxneigh);
    d_dscrfcn = k_dscrfcn.template view<DeviceType>();
    h_dscrfcn = k_dscrfcn.view_host();
    //memory->create(fcpair, maxneigh, "pair:fcpair");
    k_fcpair = DAT::tdual_kkfloat_1d("pair:fcpair", maxneigh);
    d_fcpair = k_fcpair.template view<DeviceType>();
    h_fcpair = k_fcpair.view_host();
  }

  // zero out local arrays

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagMEAMZero>(0, nall),*this);
  copymode = 0;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void
MEAMKokkos<DeviceType>::meam_dens_init(int inum_half, int ntype, typename AT::t_int_1d type, typename AT::t_int_1d d_map, typename AT::t_kkfloat_1d_3_lr x, typename AT::t_int_1d d_numneigh_half, typename AT::t_int_1d d_numneigh_full,
                     typename AT::t_int_1d d_ilist_half, typename AT::t_neighbors_2d d_neighbors_half, typename AT::t_neighbors_2d d_neighbors_full, typename AT::t_int_1d d_offset, int neighflag, int need_dup)
{
  this->ntype = ntype;
  this->type = type;
  this->d_map = d_map;
  this->x = x;
  this->d_numneigh_half = d_numneigh_half;
  this->d_numneigh_full = d_numneigh_full;
  this->d_ilist_half = d_ilist_half;
  this->d_neighbors_half = d_neighbors_half;
  this->d_neighbors_full = d_neighbors_full;
  this->d_offset = d_offset;

  if (need_dup) {
    dup_rho0 = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(d_rho0);
    dup_arho2b = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(d_arho2b);
    dup_arho1 = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(d_arho1);
    dup_arho2 = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(d_arho2);
    dup_arho3 = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(d_arho3);
    dup_arho3b = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(d_arho3b);
    dup_t_ave = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(d_t_ave);
    dup_tsq_ave = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(d_tsq_ave);
    // msmeam
    dup_arho2mb = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(d_arho2mb);
    dup_arho1m = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(d_arho1m);
    dup_arho2m = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(d_arho2m);
    dup_arho3m = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(d_arho3m);
    dup_arho3mb = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(d_arho3mb);
  } else {
    ndup_rho0 = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(d_rho0);
    ndup_arho2b = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(d_arho2b);
    ndup_arho1 = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(d_arho1);
    ndup_arho2 = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(d_arho2);
    ndup_arho3 = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(d_arho3);
    ndup_arho3b = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(d_arho3b);
    ndup_t_ave = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(d_t_ave);
    ndup_tsq_ave = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(d_tsq_ave);
    // msmeam
    ndup_arho2mb = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(d_arho2mb);
    ndup_arho1m = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(d_arho1m);
    ndup_arho2m = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(d_arho2m);
    ndup_arho3m = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(d_arho3m);
    ndup_arho3mb = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(d_arho3mb);
  }

  copymode = 1;
  if (neighflag == HALF)
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagMEAMDensInit<HALF> >(0,inum_half),*this);
  else if (neighflag == HALFTHREAD)
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagMEAMDensInit<HALFTHREAD> >(0,inum_half),*this);
  copymode = 0;

  if (need_dup) {
    Kokkos::Experimental::contribute(d_rho0, dup_rho0);
    Kokkos::Experimental::contribute(d_arho2b, dup_arho2b);
    Kokkos::Experimental::contribute(d_arho1, dup_arho1);
    Kokkos::Experimental::contribute(d_arho2, dup_arho2);
    Kokkos::Experimental::contribute(d_arho3, dup_arho3);
    Kokkos::Experimental::contribute(d_arho3b, dup_arho3b);
    Kokkos::Experimental::contribute(d_t_ave, dup_t_ave);
    Kokkos::Experimental::contribute(d_tsq_ave, dup_tsq_ave);
    // msmeam
    Kokkos::Experimental::contribute(d_arho2mb, dup_arho2mb);
    Kokkos::Experimental::contribute(d_arho1m, dup_arho1m);
    Kokkos::Experimental::contribute(d_arho2m, dup_arho2m);
    Kokkos::Experimental::contribute(d_arho3m, dup_arho3m);
    Kokkos::Experimental::contribute(d_arho3mb, dup_arho3mb);

    // free duplicated memory
    dup_rho0 = {};
    dup_arho2b = {};
    dup_arho1 = {};
    dup_arho2 = {};
    dup_arho3 = {};
    dup_arho3b = {};
    dup_t_ave = {};
    dup_tsq_ave = {};
    // msmeam
    dup_arho2mb = {};
    dup_arho1m = {};
    dup_arho2m = {};
    dup_arho3m = {};
    dup_arho3mb = {};
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG>
KOKKOS_INLINE_FUNCTION
void
MEAMKokkos<DeviceType>::getscreen(int i, int offset, typename AT::t_kkfloat_1d_3_lr x, typename AT::t_int_1d d_numneigh_half,
                typename AT::t_int_1d d_numneigh_full, int /*ntype*/, typename AT::t_int_1d type, typename AT::t_int_1d d_map)
const {
  const KK_FLOAT drinv = 1.0 / delr_meam;
  const int elti = d_map[type[i]];
  if (elti < 0) return;

  const KK_FLOAT xitmp = x(i,0);
  const KK_FLOAT yitmp = x(i,1);
  const KK_FLOAT zitmp = x(i,2);

  for (int jn = 0; jn < d_numneigh_half[i]; jn++) {
    const int j = d_neighbors_half(i,jn);

    const int eltj = d_map[type[j]];
    if (eltj < 0) continue;

    // First compute screening function itself, sij
    const KK_FLOAT xjtmp = x(j,0);
    const KK_FLOAT yjtmp = x(j,1);
    const KK_FLOAT zjtmp = x(j,2);
    const KK_FLOAT delxij = xjtmp - xitmp;
    const KK_FLOAT delyij = yjtmp - yitmp;
    const KK_FLOAT delzij = zjtmp - zitmp;

    const KK_FLOAT rij2 = delxij * delxij + delyij * delyij + delzij * delzij;

    if (rij2 > cutforcesq) {
      d_dscrfcn[offset+jn] = 0.0;
      d_scrfcn[offset+jn] = 0.0;
      d_fcpair[offset+jn] = 0.0;
      continue;
    }

    // Now compute derivatives
    const KK_FLOAT rbound = ebound_meam[elti][eltj] * rij2;
    const KK_FLOAT rij = sqrt(rij2);
    const KK_FLOAT rnorm = (cutforce - rij) * drinv;
    KK_FLOAT sij = 1.0;

    // if rjk2 > ebound*rijsq, atom k is definitely outside the ellipse
    for (int kn = 0; kn < d_numneigh_full[i]; kn++) {
      int k = d_neighbors_full(i,kn);
      if (k == j) continue;
      int eltk = d_map[type[k]];
      if (eltk < 0) continue;

      const KK_FLOAT xktmp = x(k,0);
      const KK_FLOAT yktmp = x(k,1);
      const KK_FLOAT zktmp = x(k,2);

      const KK_FLOAT delxjk = xktmp - xjtmp;
      const KK_FLOAT delyjk = yktmp - yjtmp;
      const KK_FLOAT delzjk = zktmp - zjtmp;
      const KK_FLOAT rjk2 = delxjk * delxjk + delyjk * delyjk + delzjk * delzjk;
      if (rjk2 > rbound) continue;

      const KK_FLOAT delxik = xktmp - xitmp;
      const KK_FLOAT delyik = yktmp - yitmp;
      const KK_FLOAT delzik = zktmp - zitmp;
      const KK_FLOAT rik2 = delxik * delxik + delyik * delyik + delzik * delzik;
      if (rik2 > rbound) continue;

      const KK_FLOAT xik = rik2 / rij2;
      const KK_FLOAT xjk = rjk2 / rij2;
      const KK_FLOAT a = 1 - (xik - xjk) * (xik - xjk);
      // if a < 0, then ellipse equation doesn't describe this case and
      // atom k can't possibly screen i-j
      if (a <= 0.0) continue;

      KK_FLOAT cikj = (2.0 * (xik + xjk) + a - 2.0) / a;
      const KK_FLOAT Cmax = Cmax_meam[elti][eltj][eltk];
      const KK_FLOAT Cmin = Cmin_meam[elti][eltj][eltk];
      KK_FLOAT sikj;
      if (cikj >= Cmax) continue;
      // note that cikj may be slightly negative (within numerical
      // tolerance) if atoms are colinear, so don't reject that case here
      // (other negative cikj cases were handled by the test on "a" above)
      else if (cikj <= Cmin) {
        sij = 0.0;
        break;
      } else {
        const KK_FLOAT delc = Cmax - Cmin;
        cikj = (cikj - Cmin) / delc;
        sikj = fcut(cikj);
      }
      sij *= sikj;
    }

    KK_FLOAT dfc;
    const KK_FLOAT fc = dfcut(rnorm, dfc);
    const KK_FLOAT fcij = fc;
    const KK_FLOAT dfcij = dfc * drinv;

    // Now compute derivatives
    d_dscrfcn[offset+jn] = 0.0;
    const KK_FLOAT sfcij = sij * fcij;
    if (!iszero_kk(sfcij) && !isone_kk(sfcij)) {
      for (int kn = 0; kn < d_numneigh_full[i]; kn++) {
        const int k = d_neighbors_full(i,kn);
        if (k == j) continue;
        const int eltk = d_map[type[k]];
        if (eltk < 0) continue;

        const KK_FLOAT delxjk = x(k,0) - xjtmp;
        const KK_FLOAT delyjk = x(k,1) - yjtmp;
        const KK_FLOAT delzjk = x(k,2) - zjtmp;
        const KK_FLOAT rjk2 = delxjk * delxjk + delyjk * delyjk + delzjk * delzjk;
        if (rjk2 > rbound) continue;

        const KK_FLOAT delxik = x(k,0) - xitmp;
        const KK_FLOAT delyik = x(k,1) - yitmp;
        const KK_FLOAT delzik = x(k,2) - zitmp;
        const KK_FLOAT rik2 = delxik * delxik + delyik * delyik + delzik * delzik;
        if (rik2 > rbound) continue;

        const KK_FLOAT xik = rik2 / rij2;
        const KK_FLOAT xjk = rjk2 / rij2;
        const KK_FLOAT a = 1 - (xik - xjk) * (xik - xjk);
        // if a < 0, then ellipse equation doesn't describe this case and
        // atom k can't possibly screen i-j
        if (a <= 0.0) continue;

        KK_FLOAT cikj = (2.0 * (xik + xjk) + a - 2.0) / a;
        const KK_FLOAT Cmax = Cmax_meam[elti][eltj][eltk];
        const KK_FLOAT Cmin = Cmin_meam[elti][eltj][eltk];
        if (cikj >= Cmax) {
          continue;
          // Note that cikj may be slightly negative (within numerical
          // tolerance) if atoms are colinear, so don't reject that case
          // here
          // (other negative cikj cases were handled by the test on "a"
          // above)
          // Note that we never have 0<cikj<Cmin here, else sij=0
          // (rejected above)
        } else {
          const KK_FLOAT delc = Cmax - Cmin;
          cikj = (cikj - Cmin) / delc;
          KK_FLOAT dfikj;
          const KK_FLOAT sikj = dfcut(cikj, dfikj);
          const KK_FLOAT coef1 = dfikj / (delc * sikj);
          const KK_FLOAT dCikj = dCfunc(rij2, rik2, rjk2);
          d_dscrfcn[offset+jn] += coef1 * dCikj;
        }
      }
      const KK_FLOAT coef1 = sfcij;
      const KK_FLOAT coef2 = sij * dfcij / rij;
      d_dscrfcn[offset+jn] = d_dscrfcn[offset+jn] * coef1 - coef2;
    }

    d_scrfcn[offset+jn] = sij;
    d_fcpair[offset+jn] = fcij;
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG>
KOKKOS_INLINE_FUNCTION
void
MEAMKokkos<DeviceType>::calc_rho1(int i, int /*ntype*/, typename AT::t_int_1d type, typename AT::t_int_1d d_map, typename AT::t_kkfloat_1d_3_lr x, typename AT::t_int_1d d_numneigh,
                int offset) const
{
  // The rho0, etc. arrays are duplicated for OpenMP, atomic for GPU, and neither for Serial
  auto v_rho0 = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_rho0),decltype(ndup_rho0)>::get(dup_rho0,ndup_rho0);
  auto a_rho0 = v_rho0.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();
  auto v_arho2b = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_arho2b),decltype(ndup_arho2b)>::get(dup_arho2b,ndup_arho2b);
  auto a_arho2b = v_arho2b.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();
  auto v_arho1 = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_arho1),decltype(ndup_arho1)>::get(dup_arho1,ndup_arho1);
  auto a_arho1 = v_arho1.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();
  auto v_arho2 = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_arho2),decltype(ndup_arho2)>::get(dup_arho2,ndup_arho2);
  auto a_arho2 = v_arho2.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();
  auto v_arho3 = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_arho3),decltype(ndup_arho3)>::get(dup_arho3,ndup_arho3);
  auto a_arho3 = v_arho3.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();
  auto v_arho3b = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_arho3b),decltype(ndup_arho3b)>::get(dup_arho3b,ndup_arho3b);
  auto a_arho3b = v_arho3b.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();
  auto v_t_ave = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_t_ave),decltype(ndup_t_ave)>::get(dup_t_ave,ndup_t_ave);
  auto a_t_ave = v_t_ave.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();
  auto v_tsq_ave = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_tsq_ave),decltype(ndup_tsq_ave)>::get(dup_tsq_ave,ndup_tsq_ave);
  auto a_tsq_ave = v_tsq_ave.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();
  // msmeam
  auto v_arho2mb = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_arho2mb),decltype(ndup_arho2mb)>::get(dup_arho2mb,ndup_arho2mb);
  auto a_arho2mb = v_arho2mb.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();
  auto v_arho1m = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_arho1m),decltype(ndup_arho1m)>::get(dup_arho1m,ndup_arho1m);
  auto a_arho1m = v_arho1m.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();
  auto v_arho2m = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_arho2m),decltype(ndup_arho2m)>::get(dup_arho2m,ndup_arho2m);
  auto a_arho2m = v_arho2m.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();
  auto v_arho3m = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_arho3m),decltype(ndup_arho3m)>::get(dup_arho3m,ndup_arho3m);
  auto a_arho3m = v_arho3m.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();
  auto v_arho3mb = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_arho3mb),decltype(ndup_arho3mb)>::get(dup_arho3mb,ndup_arho3mb);
  auto a_arho3mb = v_arho3mb.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  const int elti = d_map[type[i]];
  const KK_FLOAT xtmp = x(i,0);
  const KK_FLOAT ytmp = x(i,1);
  const KK_FLOAT ztmp = x(i,2);
  for (int jn = 0; jn < d_numneigh[i]; jn++) {
    if (!iszero_kk(d_scrfcn[offset+jn])) {
      const int j = d_neighbors_half(i,jn);
      const KK_FLOAT sij = d_scrfcn[offset+jn] * d_fcpair[offset+jn];
      KK_FLOAT delij[3];
      delij[0] = x(j,0) - xtmp;
      delij[1] = x(j,1) - ytmp;
      delij[2] = x(j,2) - ztmp;
      const KK_FLOAT rij2 = delij[0] * delij[0] + delij[1] * delij[1] + delij[2] * delij[2];
      if (rij2 < cutforcesq) {
        const int eltj = d_map[type[j]];
        const KK_FLOAT rij = sqrt(rij2);
        const KK_FLOAT ai = rij / re_meam[elti][elti] - 1.0;
        const KK_FLOAT aj = rij / re_meam[eltj][eltj] - 1.0;
        const KK_FLOAT ro0i = rho0_meam[elti];
        const KK_FLOAT ro0j = rho0_meam[eltj];
        const KK_FLOAT rhoa0j = ro0j * MathSpecialKokkos::fm_exp(-beta0_meam[eltj] * aj) * sij;
        KK_FLOAT rhoa1j = ro0j * MathSpecialKokkos::fm_exp(-beta1_meam[eltj] * aj) * sij;
        KK_FLOAT rhoa2j = ro0j * MathSpecialKokkos::fm_exp(-beta2_meam[eltj] * aj) * sij;
        KK_FLOAT rhoa3j = ro0j * MathSpecialKokkos::fm_exp(-beta3_meam[eltj] * aj) * sij;
        const KK_FLOAT rhoa0i = ro0i * MathSpecialKokkos::fm_exp(-beta0_meam[elti] * ai) * sij;
        KK_FLOAT rhoa1i = ro0i * MathSpecialKokkos::fm_exp(-beta1_meam[elti] * ai) * sij;
        KK_FLOAT rhoa2i = ro0i * MathSpecialKokkos::fm_exp(-beta2_meam[elti] * ai) * sij;
        KK_FLOAT rhoa3i = ro0i * MathSpecialKokkos::fm_exp(-beta3_meam[elti] * ai) * sij;
        // msmeam
        KK_FLOAT rhoa1mj, rhoa2mj, rhoa3mj, rhoa1mi, rhoa2mi, rhoa3mi;
        if (msmeamflag) {
          rhoa1mj = ro0j * t1m_meam[eltj] * MathSpecialKokkos::fm_exp(-beta1m_meam[eltj] * aj) * sij;
          rhoa2mj = ro0j * t2m_meam[eltj] * MathSpecialKokkos::fm_exp(-beta2m_meam[eltj] * aj) * sij;
          rhoa3mj = ro0j * t3m_meam[eltj] * MathSpecialKokkos::fm_exp(-beta3m_meam[eltj] * aj) * sij;
          rhoa1mi = ro0i * t1m_meam[elti] * MathSpecialKokkos::fm_exp(-beta1m_meam[elti] * ai) * sij;
          rhoa2mi = ro0i * t2m_meam[elti] * MathSpecialKokkos::fm_exp(-beta2m_meam[elti] * ai) * sij;
          rhoa3mi = ro0i * t3m_meam[elti] * MathSpecialKokkos::fm_exp(-beta3m_meam[elti] * ai) * sij;
        }
        if (ialloy == 1) {
          rhoa1j *= t1_meam[eltj];
          rhoa2j *= t2_meam[eltj];
          rhoa3j *= t3_meam[eltj];
          rhoa1i *= t1_meam[elti];
          rhoa2i *= t2_meam[elti];
          rhoa3i *= t3_meam[elti];
        }
        a_rho0[i] += rhoa0j;
        a_rho0[j] += rhoa0i;
        // For ialloy = 2, use single-element value (not average)
        if (ialloy != 2) {
          a_t_ave(i,0) += t1_meam[eltj] * rhoa0j;
          a_t_ave(i,1) += t2_meam[eltj] * rhoa0j;
          a_t_ave(i,2) += t3_meam[eltj] * rhoa0j;
          a_t_ave(j,0) += t1_meam[elti] * rhoa0i;
          a_t_ave(j,1) += t2_meam[elti] * rhoa0i;
          a_t_ave(j,2) += t3_meam[elti] * rhoa0i;
        }
        if (ialloy == 1) {
          a_tsq_ave(i,0) += t1_meam[eltj] * t1_meam[eltj] * rhoa0j;
          a_tsq_ave(i,1) += t2_meam[eltj] * t2_meam[eltj] * rhoa0j;
          a_tsq_ave(i,2) += t3_meam[eltj] * t3_meam[eltj] * rhoa0j;
          a_tsq_ave(j,0) += t1_meam[elti] * t1_meam[elti] * rhoa0i;
          a_tsq_ave(j,1) += t2_meam[elti] * t2_meam[elti] * rhoa0i;
          a_tsq_ave(j,2) += t3_meam[elti] * t3_meam[elti] * rhoa0i;
        }
        a_arho2b[i] += rhoa2j;
        a_arho2b[j] += rhoa2i;

        const KK_FLOAT A1j = rhoa1j / rij;
        const KK_FLOAT A2j = rhoa2j / rij2;
        const KK_FLOAT A3j = rhoa3j / (rij2 * rij);
        const KK_FLOAT A1i = rhoa1i / rij;
        const KK_FLOAT A2i = rhoa2i / rij2;
        const KK_FLOAT A3i = rhoa3i / (rij2 * rij);
        KK_FLOAT A1mj, A2mj, A3mj, A1mi, A2mi, A3mi;
        if (msmeamflag) {
          a_arho2mb[i] += rhoa2mj;
          a_arho2mb[j] += rhoa2mi;
          A1mj = rhoa1mj / rij;
          A2mj = rhoa2mj / rij2;
          A3mj = rhoa3mj / (rij2 * rij);
          A1mi = rhoa1mi / rij;
          A2mi = rhoa2mi / rij2;
          A3mi = rhoa3mi / (rij2 * rij);
        }
        int nv2 = 0;
        int nv3 = 0;
        for (int m = 0; m < 3; m++) {
          a_arho1(i,m) +=  A1j * delij[m];
          a_arho1(j,m) += -A1i * delij[m];
          a_arho3b(i,m) +=  rhoa3j * delij[m] / rij;
          a_arho3b(j,m) += -rhoa3i * delij[m] / rij;
          if (msmeamflag) {
            a_arho1m(i,m) +=  A1mj * delij[m];
            a_arho1m(j,m) += -A1mi * delij[m];
            a_arho3mb(i,m) +=  rhoa3mj * delij[m] / rij;
            a_arho3mb(j,m) += -rhoa3mi * delij[m] / rij;
          }
          for (int n = m; n < 3; n++) {
            a_arho2(i,nv2) += A2j * delij[m] * delij[n];
            a_arho2(j,nv2) += A2i * delij[m] * delij[n];
            if (msmeamflag) {
              a_arho2m(i,nv2) += A2mj * delij[m] * delij[n];
              a_arho2m(j,nv2) += A2mi * delij[m] * delij[n];
            }
            nv2++;
            for (int p = n; p < 3; p++) {
              a_arho3(i,nv3) +=  A3j * delij[m] * delij[n] * delij[p];
              a_arho3(j,nv3) += -A3i * delij[m] * delij[n] * delij[p];
              if (msmeamflag) {
                a_arho3m(i,nv3) +=  A3mj * delij[m] * delij[n] * delij[p];
                a_arho3m(j,nv3) += -A3mi * delij[m] * delij[n] * delij[p];
              }
              nv3++;
            }
          }
        }
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

//Cutoff function and derivative

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
KK_FLOAT MEAMKokkos<DeviceType>::dfcut(const KK_FLOAT xi, KK_FLOAT& dfc) const
{
  if (xi >= 1.0) {
    dfc = 0.0;
    return 1.0;
  } else if (xi <= 0.0) {
    dfc = 0.0;
    return 0.0;
  } else {
    const KK_FLOAT a = 1.0 - xi;
    const KK_FLOAT a3 = a * a * a;
    const KK_FLOAT a4 = a * a3;
    const KK_FLOAT a1m4 = 1.0 - a4;

    dfc = 8 * a1m4 * a3;
    return a1m4*a1m4;
  }
}

  //-----------------------------------------------------------------------------
  // Derivative of Cikj w.r.t. rij
  // Inputs: rij,rij2,rik2,rjk2

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
KK_FLOAT MEAMKokkos<DeviceType>::dCfunc(const KK_FLOAT rij2, const KK_FLOAT rik2, const KK_FLOAT rjk2) const
{
  const KK_FLOAT rij4 = rij2 * rij2;
  const KK_FLOAT a = rik2 - rjk2;
  const KK_FLOAT b = rik2 + rjk2;
  const KK_FLOAT asq = a*a;
  KK_FLOAT denom = rij4 - asq;
  denom = denom * denom;
  return -4 * (-2 * rij2 * asq + rij4 * b + asq * b) / denom;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void MEAMKokkos<DeviceType>::dCfunc2(const KK_FLOAT rij2, const KK_FLOAT rik2, const KK_FLOAT rjk2, KK_FLOAT& dCikj1, KK_FLOAT& dCikj2) const
{
  const KK_FLOAT rij4 = rij2 * rij2;
  const KK_FLOAT rik4 = rik2 * rik2;
  const KK_FLOAT rjk4 = rjk2 * rjk2;
  const KK_FLOAT a = rik2 - rjk2;
  KK_FLOAT denom = rij4 - a * a;
  denom = denom * denom;
  dCikj1 = 4 * rij2 * (rij4 + rik4 + 2 * rik2 * rjk2 - 3 * rjk4 - 2 * rij2 * a) / denom;
  dCikj2 = 4 * rij2 * (rij4 - 3 * rik4 + 2 * rik2 * rjk2 + rjk4 + 2 * rij2 * a) / denom;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
KK_FLOAT MEAMKokkos<DeviceType>::fcut(const KK_FLOAT xi) const
{
  KK_FLOAT a;
  if (xi >= 1.0)
    return 1.0;
  else if (xi <= 0.0)
    return 0.0;
  else {
    // ( 1.d0 - (1.d0 - xi)**4 )**2, but with better codegen
    a = 1.0 - xi;
    a *= a; a *= a;
    a = 1.0 - a;
    return a * a;
  }
}

