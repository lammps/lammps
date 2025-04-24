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

#ifndef LMP_PAIR_REAXFF_NEIGHBOR_KOKKOS_HPP
#define LMP_PAIR_REAXFF_NEIGHBOR_KOKKOS_HPP

#include "pair_reaxff_kokkos.h"

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::operator()(TagPairReaxBuildListsHalfBlocking<NEIGHFLAG>, const int &ii) const {
  constexpr int blocksize = PairReaxFFKokkos<DeviceType>::build_lists_half_blocksize;

  const auto v_dDeltap_self = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_dDeltap_self),decltype(ndup_dDeltap_self)>::get(dup_dDeltap_self,ndup_dDeltap_self);
  const auto a_dDeltap_self = v_dDeltap_self.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  const auto v_total_bo = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_total_bo),decltype(ndup_total_bo)>::get(dup_total_bo,ndup_total_bo);
  const auto a_total_bo = v_total_bo.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  const int i = d_ilist[ii];
  const X_FLOAT xtmp = x(i,0);
  const X_FLOAT ytmp = x(i,1);
  const X_FLOAT ztmp = x(i,2);
  const int itype = type(i);
  const int jnum = d_numneigh[i];

  F_FLOAT C12, C34, C56, BO_s, BO_pi, BO_pi2, BO, delij[3], dBOp_i[3];
  F_FLOAT dDeltap_self_i[3] = {0.0,0.0,0.0};
  F_FLOAT total_bo_i = 0.0;

  int ihb = -1;

  if (cut_hbsq > 0.0)
    ihb = paramssing(itype).p_hbond;

  int nnz;
  blocking_t selected_jj[blocksize];
  int jj_current = 0;

  while (jj_current < jnum) {
    nnz = 0;

    while (nnz < blocksize) {
      int jj = jj_current;
      int j = d_neighbors(i,jj);
      j &= NEIGHMASK;

      delij[0] = x(j,0) - xtmp;
      delij[1] = x(j,1) - ytmp;
      delij[2] = x(j,2) - ztmp;
      const F_FLOAT rsq = delij[0]*delij[0] + delij[1]*delij[1] + delij[2]*delij[2];

      double cutoffsq;
      if (i < nlocal) cutoffsq = MAX(cut_bosq,cut_hbsq);
      else cutoffsq = cut_bosq;
      if (rsq <= cutoffsq) {
        selected_jj[nnz] = jj_current;
        nnz++;
      }
      jj_current++;

      if (jj_current == jnum) break;
    }

    for (int jj_inner = 0; jj_inner < nnz; jj_inner++) {
      const int jj = selected_jj[jj_inner];
      int j = d_neighbors(i,jj);
      j &= NEIGHMASK;
      const int jtype = type(j);
      delij[0] = x(j,0) - xtmp;
      delij[1] = x(j,1) - ytmp;
      delij[2] = x(j,2) - ztmp;
      const F_FLOAT rsq = delij[0]*delij[0] + delij[1]*delij[1] + delij[2]*delij[2];

      // hbond list
      build_hb_list<NEIGHFLAG>(rsq, i, ihb, j, jtype);

      if (rsq > cut_bosq) continue;

      // bond_list
      const F_FLOAT rij = sqrt(rsq);
      const F_FLOAT p_bo2 = paramstwbp(itype,jtype).p_bo2;
      const F_FLOAT p_bo4 = paramstwbp(itype,jtype).p_bo4;
      const F_FLOAT p_bo6 = paramstwbp(itype,jtype).p_bo6;

      // returns BO_*, C** by reference
      compute_bo(rij, itype, jtype, p_bo2, p_bo4, p_bo6,
        BO_s, BO_pi, BO_pi2, C12, C34, C56);

      BO = BO_s + BO_pi + BO_pi2;
      if (BO < bo_cut) continue;

      int i_index = -1;
      int j_index = -1;
      if (build_bo_list<NEIGHFLAG>(i, j, i_index, j_index)) {

        // from BondOrder1

        d_BO(i,j_index) = BO;
        d_BO_s(i,j_index) = BO_s;

        d_BO(j,i_index) = BO;
        d_BO_s(j,i_index) = BO_s;

        d_BO_pi(j,i_index) = BO_pi;
        d_BO_pi2(j,i_index) = BO_pi2;

        d_BO_pi(i,j_index) = BO_pi;
        d_BO_pi2(i,j_index) = BO_pi2;

        F_FLOAT Cln_BOp_s = p_bo2 * C12 / rij / rij;
        F_FLOAT Cln_BOp_pi = p_bo4 * C34 / rij / rij;
        F_FLOAT Cln_BOp_pi2 = p_bo6 * C56 / rij / rij;

        if (nlocal == 0)
          Cln_BOp_s = Cln_BOp_pi = Cln_BOp_pi2 = 0.0;

        for (int d = 0; d < 3; d++) dBOp_i[d] = -(BO_s*Cln_BOp_s+BO_pi*Cln_BOp_pi+BO_pi2*Cln_BOp_pi2)*delij[d];
        for (int d = 0; d < 3; d++) dDeltap_self_i[d] += dBOp_i[d];
        for (int d = 0; d < 3; d++) a_dDeltap_self(j,d) += -dBOp_i[d];

        d_dln_BOp_pi(i,j_index) = -(BO_pi*Cln_BOp_pi);
        d_dln_BOp_pi(j,i_index) = -(BO_pi*Cln_BOp_pi);

        d_dln_BOp_pi2(i,j_index) = -(BO_pi2*Cln_BOp_pi2);
        d_dln_BOp_pi2(j,i_index) = -(BO_pi2*Cln_BOp_pi2);

        d_dBOp(i,j_index) = -(BO_s*Cln_BOp_s+BO_pi*Cln_BOp_pi+BO_pi2*Cln_BOp_pi2);
        d_dBOp(j,i_index) = -(BO_s*Cln_BOp_s+BO_pi*Cln_BOp_pi+BO_pi2*Cln_BOp_pi2);
        d_BO(i,j_index) = BO - bo_cut;
        d_BO(j,i_index) = BO - bo_cut;
        d_BO_s(i,j_index) = BO_s - bo_cut;
        d_BO_s(j,i_index) = BO_s - bo_cut;
        total_bo_i += (BO - bo_cut);
        a_total_bo[j] += (BO - bo_cut);
      }
    }
  }

  for (int d = 0; d < 3; d++)
    a_dDeltap_self(i,d) += dDeltap_self_i[d];

  a_total_bo[i] += total_bo_i;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::operator()(TagPairReaxBuildListsHalfBlockingPreview<NEIGHFLAG>, const int &ii) const {
  constexpr int blocksize = PairReaxFFKokkos<DeviceType>::build_lists_half_blocksize;

  const int i = d_ilist[ii];
  const X_FLOAT xtmp = x(i,0);
  const X_FLOAT ytmp = x(i,1);
  const X_FLOAT ztmp = x(i,2);
  const int itype = type(i);
  const int jnum = d_numneigh[i];

  F_FLOAT C12, C34, C56, BO_s, BO_pi, BO_pi2, BO, delij[3];

  int ihb = -1;

  if (cut_hbsq > 0.0)
    ihb = paramssing(itype).p_hbond;

  int nnz;
  blocking_t selected_jj[blocksize];
  int jj_current = 0;

  double cutoffsq;
  if (i < nlocal) cutoffsq = MAX(cut_bosq,cut_hbsq);
  else cutoffsq = cut_bosq;

  while (jj_current < jnum) {
    nnz = 0;

    while (nnz < blocksize) {
      int jj = jj_current;
      int j = d_neighbors(i,jj);
      j &= NEIGHMASK;

      delij[0] = x(j,0) - xtmp;
      delij[1] = x(j,1) - ytmp;
      delij[2] = x(j,2) - ztmp;
      const F_FLOAT rsq = delij[0]*delij[0] + delij[1]*delij[1] + delij[2]*delij[2];

      if (rsq <= cutoffsq) {
        selected_jj[nnz] = jj_current;
        nnz++;
      }
      jj_current++;

      if (jj_current == jnum) break;
    }

    for (int jj_inner = 0; jj_inner < nnz; jj_inner++) {
      const int jj = selected_jj[jj_inner];
      int j = d_neighbors(i,jj);
      j &= NEIGHMASK;
      const int jtype = type(j);
      delij[0] = x(j,0) - xtmp;
      delij[1] = x(j,1) - ytmp;
      delij[2] = x(j,2) - ztmp;
      const F_FLOAT rsq = delij[0]*delij[0] + delij[1]*delij[1] + delij[2]*delij[2];

      // hbond list
      build_hb_list<NEIGHFLAG>(rsq, i, ihb, j, jtype);

      if (rsq > cut_bosq) continue;

      // bond_list
      const F_FLOAT rij = sqrt(rsq);
      const F_FLOAT p_bo2 = paramstwbp(itype,jtype).p_bo2;
      const F_FLOAT p_bo4 = paramstwbp(itype,jtype).p_bo4;
      const F_FLOAT p_bo6 = paramstwbp(itype,jtype).p_bo6;

      // returns BO_*, C** by reference
      compute_bo(rij, itype, jtype, p_bo2, p_bo4, p_bo6,
        BO_s, BO_pi, BO_pi2, C12, C34, C56);

      BO = BO_s + BO_pi + BO_pi2;
      if (BO < bo_cut) continue;

      int i_index = -1;
      int j_index = -1;
      build_bo_list<NEIGHFLAG>(i, j, i_index, j_index);
    }
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::operator()(TagPairReaxBuildListsHalfPreview<NEIGHFLAG>, const int &ii) const {

  const int i = d_ilist[ii];
  const X_FLOAT xtmp = x(i,0);
  const X_FLOAT ytmp = x(i,1);
  const X_FLOAT ztmp = x(i,2);
  const int itype = type(i);
  const int jnum = d_numneigh[i];

  F_FLOAT C12, C34, C56, BO_s, BO_pi, BO_pi2, BO, delij[3];

  int ihb = -1;

  if (cut_hbsq > 0.0)
    ihb = paramssing(itype).p_hbond;

  for (int jj = 0; jj < jnum; jj++) {
    int j = d_neighbors(i,jj);
    j &= NEIGHMASK;

    const int jtype = type(j);

    delij[0] = x(j,0) - xtmp;
    delij[1] = x(j,1) - ytmp;
    delij[2] = x(j,2) - ztmp;
    const F_FLOAT rsq = delij[0]*delij[0] + delij[1]*delij[1] + delij[2]*delij[2];

    // hbond list
    build_hb_list<NEIGHFLAG>(rsq, i, ihb, j, jtype);

    if (rsq > cut_bosq) continue;

    // bond_list
    const F_FLOAT rij = sqrt(rsq);
    const F_FLOAT p_bo2 = paramstwbp(itype,jtype).p_bo2;
    const F_FLOAT p_bo4 = paramstwbp(itype,jtype).p_bo4;
    const F_FLOAT p_bo6 = paramstwbp(itype,jtype).p_bo6;

    // returns BO_*, C** by reference
    compute_bo(rij, itype, jtype, p_bo2, p_bo4, p_bo6,
      BO_s, BO_pi, BO_pi2, C12, C34, C56);

    BO = BO_s + BO_pi + BO_pi2;
    if (BO < bo_cut) continue;

    int i_index = -1;
    int j_index = -1;

    build_bo_list<NEIGHFLAG>(i, j, i_index, j_index);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::build_hb_list(F_FLOAT rsq, int i, int ihb, int j, int jtype) const {

  int i_index, j_index;
  int jhb = -1;
  if (i < nlocal && cut_hbsq > 0.0 && (ihb == 1 || ihb == 2) && rsq <= cut_hbsq) {
    jhb = paramssing(jtype).p_hbond;
    if (ihb == 1 && jhb == 2) {
      if (NEIGHFLAG == HALF) {
        j_index = d_hb_num[i];
        d_hb_num[i]++;
      } else
        j_index = Kokkos::atomic_fetch_add(&d_hb_num[i],1);

      if (j_index >= maxhb)
        d_resize_hb() = MAX(d_resize_hb(), j_index+1);
      else
        d_hb_list(i, j_index) = j;
    } else if (j < nlocal && ihb == 2 && jhb == 1) {
      if (NEIGHFLAG == HALF) {
        i_index = d_hb_num[j];
        d_hb_num[j]++;
      } else
        i_index = Kokkos::atomic_fetch_add(&d_hb_num[j],1);

      if (i_index >= maxhb)
        d_resize_hb() = MAX(d_resize_hb(), i_index+1);
      else
        d_hb_list(j, i_index) = i;
    }
  }

}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG>
KOKKOS_INLINE_FUNCTION
bool PairReaxFFKokkos<DeviceType>::build_bo_list(int i, int j, int& i_index, int& j_index) const {

  if (NEIGHFLAG == HALF) {
    j_index = d_bo_num[i];
    i_index = d_bo_num[j];
    d_bo_num[i]++;
    d_bo_num[j]++;
  } else {
    j_index = Kokkos::atomic_fetch_add(&d_bo_num[i],1);
    i_index = Kokkos::atomic_fetch_add(&d_bo_num[j],1);
  }

  bool set_dB_flag = true;

  if (j_index >= maxbo || i_index >= maxbo) {
    const int max_val = MAX(i_index + 1, j_index + 1);
    d_resize_bo() = MAX(d_resize_bo(),max_val);
    set_dB_flag = false;
  } else {
    d_bo_list(i, j_index) = j;
    d_bo_list(j, i_index) = i;
    set_dB_flag = true;
  }

  return set_dB_flag;

}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::operator()(TagPairReaxBuildListsFull, const int &ii) const {

  const int i = d_ilist[ii];
  const X_FLOAT xtmp = x(i,0);
  const X_FLOAT ytmp = x(i,1);
  const X_FLOAT ztmp = x(i,2);
  const int itype = type(i);

  F_FLOAT C12, C34, C56, BO_s, BO_pi, BO_pi2, BO, delij[3], dBOp_i[3];
  F_FLOAT dDeltap_self_i[3] = {0.0,0.0,0.0};
  F_FLOAT total_bo_i = 0.0;

  const int jnum = d_bo_num[i];
  for (int j_index = 0; j_index < jnum; j_index++) {
    int j = d_bo_list(i, j_index);
    j &= NEIGHMASK;
    const int jtype = type(j);
    delij[0] = x(j,0) - xtmp;
    delij[1] = x(j,1) - ytmp;
    delij[2] = x(j,2) - ztmp;
    const F_FLOAT rsq = delij[0]*delij[0] + delij[1]*delij[1] + delij[2]*delij[2];
    const F_FLOAT rsq_inv = 1.0 / rsq;

    // bond_list
    const F_FLOAT rij = sqrt(rsq);
    const F_FLOAT p_bo2 = paramstwbp(itype,jtype).p_bo2;
    const F_FLOAT p_bo4 = paramstwbp(itype,jtype).p_bo4;
    const F_FLOAT p_bo6 = paramstwbp(itype,jtype).p_bo6;

    // returns BO_*, C** by reference
    compute_bo(rij, itype, jtype, p_bo2, p_bo4, p_bo6,
      BO_s, BO_pi, BO_pi2, C12, C34, C56);

    BO = BO_s + BO_pi + BO_pi2;

    // from BondOrder1

    d_BO(i,j_index) = BO;
    d_BO_s(i,j_index) = BO_s;
    d_BO_pi(i,j_index) = BO_pi;
    d_BO_pi2(i,j_index) = BO_pi2;

    F_FLOAT Cln_BOp_s = p_bo2 * C12 * rsq_inv;
    F_FLOAT Cln_BOp_pi = p_bo4 * C34 * rsq_inv;
    F_FLOAT Cln_BOp_pi2 = p_bo6 * C56 * rsq_inv;

    if (nlocal == 0)
      Cln_BOp_s = Cln_BOp_pi = Cln_BOp_pi2 = 0.0;

    for (int d = 0; d < 3; d++) dBOp_i[d] = -(BO_s*Cln_BOp_s+BO_pi*Cln_BOp_pi+BO_pi2*Cln_BOp_pi2)*delij[d];
    for (int d = 0; d < 3; d++) dDeltap_self_i[d] += dBOp_i[d];


    d_dln_BOp_pi(i,j_index) = -(BO_pi*Cln_BOp_pi);
    d_dln_BOp_pi2(i,j_index) = -(BO_pi2*Cln_BOp_pi2);
    d_dBOp(i,j_index) = -(BO_s*Cln_BOp_s+BO_pi*Cln_BOp_pi+BO_pi2*Cln_BOp_pi2);

    d_BO(i,j_index) = BO - bo_cut;
    d_BO_s(i,j_index) = BO_s - bo_cut;
    total_bo_i += (BO - bo_cut);
  }

  for (int d = 0; d < 3; d++)
    d_dDeltap_self(i,d) = dDeltap_self_i[d];

  d_total_bo[i] = total_bo_i;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::compute_bo(F_FLOAT rij, int itype, int jtype, F_FLOAT p_bo2, F_FLOAT p_bo4, F_FLOAT p_bo6,
  F_FLOAT& BO_s, F_FLOAT& BO_pi, F_FLOAT& BO_pi2, F_FLOAT& C12, F_FLOAT& C34, F_FLOAT& C56) const {

  const F_FLOAT p_bo1 = paramstwbp(itype,jtype).p_bo1;
  const F_FLOAT p_bo3 = paramstwbp(itype,jtype).p_bo3;
  const F_FLOAT p_bo5 = paramstwbp(itype,jtype).p_bo5;
  const F_FLOAT r_s = paramstwbp(itype,jtype).r_s;
  const F_FLOAT r_p = paramstwbp(itype,jtype).r_p;
  const F_FLOAT r_pp = paramstwbp(itype,jtype).r_pp;

  if (paramssing(itype).r_s > 0.0  && paramssing(jtype).r_s > 0.0) {
    C12 = p_bo1 * ((p_bo2 != 0) ? (pow(rij/r_s,p_bo2)) : 1.0);
    BO_s = (1.0+bo_cut)*exp(C12);
  } else BO_s = C12 = 0.0;

  if (paramssing(itype).r_p > 0.0  && paramssing(jtype).r_p > 0.0) {
    C34 = p_bo3 * ((p_bo4 != 0) ? (pow(rij/r_p,p_bo4)) : 1.0);
    BO_pi = exp(C34);
  } else BO_pi = C34 = 0.0;

  if (paramssing(itype).r_pp > 0.0  && paramssing(jtype).r_pp > 0.0) {
    C56 = p_bo5 * ((p_bo6 != 0) ? (pow(rij/r_pp,p_bo6)) : 1.0);
    BO_pi2 = exp(C56);
  } else BO_pi2 = C56 = 0.0;

}





#endif
