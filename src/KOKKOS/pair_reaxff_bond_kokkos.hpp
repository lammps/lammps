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

#ifndef LMP_PAIR_REAXFF_BOND_KOKKOS_HPP
#define LMP_PAIR_REAXFF_BOND_KOKKOS_HPP

#include "pair_reaxff_kokkos.h"

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::operator()(TagPairReaxUpdateBond<NEIGHFLAG>, const int &ii) const {

  Kokkos::View<F_FLOAT**, typename decltype(d_Cdbo)::array_layout,KKDeviceType,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value>> a_Cdbo = d_Cdbo;
  Kokkos::View<F_FLOAT**, typename decltype(d_Cdbopi)::array_layout,KKDeviceType,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value>> a_Cdbopi = d_Cdbopi;
  Kokkos::View<F_FLOAT**, typename decltype(d_Cdbopi2)::array_layout,KKDeviceType,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value>> a_Cdbopi2 = d_Cdbopi2;

  const int i = d_ilist[ii];
  const X_FLOAT xtmp = x(i,0);
  const X_FLOAT ytmp = x(i,1);
  const X_FLOAT ztmp = x(i,2);
  const tagint itag = tag(i);
  const int jnum = d_bo_num[i];

  for (int j_index = 0; j_index < jnum; j_index++) {
    int j = d_bo_list(i, j_index);
    j &= NEIGHMASK;
    const tagint jtag = tag(j);

    int flag = 0;

    if (itag > jtag) {
      if ((itag+jtag) % 2 == 0) flag = 1;
    } else if (itag < jtag) {
      if ((itag+jtag) % 2 == 1) flag = 1;
    } else {
      if (x(j,2)  < ztmp) flag = 1;
      if (x(j,2) == ztmp && x(j,1)  < ytmp) flag = 1;
      if (x(j,2) == ztmp && x(j,1) == ytmp && x(j,0) < xtmp) flag = 1;
    }

    if (!flag) continue;

    const F_FLOAT Cdbo_i = d_Cdbo(i,j_index);
    const F_FLOAT Cdbopi_i = d_Cdbopi(i,j_index);
    const F_FLOAT Cdbopi2_i = d_Cdbopi2(i,j_index);

    const int knum = d_bo_num[j];

    for (int k_index = 0; k_index < knum; k_index++) {
      int k = d_bo_list(j, k_index);
      k &= NEIGHMASK;
      if (k != i) continue;

      a_Cdbo(j,k_index) += Cdbo_i;
      a_Cdbopi(j,k_index) += Cdbopi_i;
      a_Cdbopi2(j,k_index) += Cdbopi2_i;
    }
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG, int EFLAG>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::operator()(TagPairReaxComputeBond1<NEIGHFLAG,EFLAG>, const int &ii, EV_FLOAT_REAX& ev) const {

  auto v_CdDelta = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_CdDelta),decltype(ndup_CdDelta)>::get(dup_CdDelta,ndup_CdDelta);
  auto a_CdDelta = v_CdDelta.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  F_FLOAT p_be1, p_be2, De_s, De_p, De_pp, pow_BOs_be2, exp_be12, CEbo, ebond;

  const int i = d_ilist[ii];
  const X_FLOAT xtmp = x(i,0);
  const X_FLOAT ytmp = x(i,1);
  const X_FLOAT ztmp = x(i,2);
  const int itype = type(i);
  const tagint itag = tag(i);
  const F_FLOAT imass = paramssing(itype).mass;
  const int jnum = d_bo_num[i];

  F_FLOAT CdDelta_i = 0.0;

  for (int j_index = 0; j_index < jnum; j_index++) {
    int j = d_bo_list(i, j_index);
    j &= NEIGHMASK;
    const tagint jtag = tag(j);

    if (itag > jtag) {
      if ((itag+jtag) % 2 == 0) continue;
    } else if (itag < jtag) {
      if ((itag+jtag) % 2 == 1) continue;
    } else {
      if (x(j,2)  < ztmp) continue;
      if (x(j,2) == ztmp && x(j,1)  < ytmp) continue;
      if (x(j,2) == ztmp && x(j,1) == ytmp && x(j,0) < xtmp) continue;
    }

    const int jtype = type(j);
    const F_FLOAT jmass = paramssing(jtype).mass;

    // bond energy (nlocal only)
    p_be1 = paramstwbp(itype,jtype).p_be1;
    p_be2 = paramstwbp(itype,jtype).p_be2;
    De_s = paramstwbp(itype,jtype).De_s;
    De_p = paramstwbp(itype,jtype).De_p;
    De_pp = paramstwbp(itype,jtype).De_pp;

    const F_FLOAT BO_i = d_BO(i,j_index);
    const F_FLOAT BO_s_i = d_BO_s(i,j_index);
    const F_FLOAT BO_pi_i = d_BO_pi(i,j_index);
    const F_FLOAT BO_pi2_i = d_BO_pi2(i,j_index);

    if (BO_s_i == 0.0) pow_BOs_be2 = 0.0;
    else pow_BOs_be2 = pow(BO_s_i,p_be2);
    exp_be12 = exp(p_be1*(1.0-pow_BOs_be2));
    CEbo = -De_s*exp_be12*(1.0-p_be1*p_be2*pow_BOs_be2);
    ebond = -De_s*BO_s_i*exp_be12
                              -De_p*BO_pi_i
                          -De_pp*BO_pi2_i;

    if (EFLAG && eflag_global) ev.evdwl += ebond;
    //if (eflag_atom) this->template ev_tally<NEIGHFLAG>(ev,i,j,ebond,0.0,0.0,0.0,0.0);
    //if (eflag_atom) this->template e_tally<NEIGHFLAG>(ev,i,j,ebond);

    // calculate derivatives of Bond Orders
    d_Cdbo(i,j_index) += CEbo;
    d_Cdbopi(i,j_index) -= (CEbo + De_p);
    d_Cdbopi2(i,j_index) -= (CEbo + De_pp);

    // Stabilisation terminal triple bond
    F_FLOAT estriph = 0.0;

    if (BO_i >= 1.00) {
      if (gp[37] == 2 || (imass == 12.0000 && jmass == 15.9990) ||
                         (jmass == 12.0000 && imass == 15.9990)) {
        const F_FLOAT exphu = exp(-gp[7] * SQR(BO_i - 2.50));
        const F_FLOAT exphua1 = exp(-gp[3] * (d_total_bo[i]-BO_i));
        const F_FLOAT exphub1 = exp(-gp[3] * (d_total_bo[j]-BO_i));
        const F_FLOAT exphuov = exp(gp[4] * (d_Delta[i] + d_Delta[j]));
        const F_FLOAT hulpov = 1.0 / (1.0 + 25.0 * exphuov);
        estriph = gp[10] * exphu * hulpov * (exphua1 + exphub1);

        if (EFLAG && eflag_global) ev.evdwl += estriph;
        //if (eflag_atom) this->template ev_tally<NEIGHFLAG>(ev,i,j,estriph,0.0,0.0,0.0,0.0);
        //if (eflag_atom) this->template e_tally<NEIGHFLAG>(ev,i,j,estriph);

        const F_FLOAT decobdbo = gp[10] * exphu * hulpov * (exphua1 + exphub1) *
            (gp[3] - 2.0 * gp[7] * (BO_i-2.50));
        const F_FLOAT decobdboua = -gp[10] * exphu * hulpov *
            (gp[3]*exphua1 + 25.0*gp[4]*exphuov*hulpov*(exphua1+exphub1));
        const F_FLOAT decobdboub = -gp[10] * exphu * hulpov *
            (gp[3]*exphub1 + 25.0*gp[4]*exphuov*hulpov*(exphua1+exphub1));

        d_Cdbo(i,j_index) += decobdbo;
        CdDelta_i += decobdboua;
        a_CdDelta[j] += decobdboub;
      }
    }
    const F_FLOAT eng_tmp = ebond + estriph;
    if (eflag_atom) this->template e_tally<NEIGHFLAG>(ev,i,j,eng_tmp);
  }
  a_CdDelta[i] += CdDelta_i;
}

template<class DeviceType>
template<int NEIGHFLAG, int EFLAG>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::operator()(TagPairReaxComputeBond1<NEIGHFLAG,EFLAG>, const int &ii) const {
  EV_FLOAT_REAX ev;
  this->template operator()<NEIGHFLAG,EFLAG>(TagPairReaxComputeBond1<NEIGHFLAG,EFLAG>(), ii, ev);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG, int VFLAG>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::operator()(TagPairReaxComputeBond2<NEIGHFLAG,VFLAG>, const int &ii, EV_FLOAT_REAX& ev) const {

  auto v_f = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_f),decltype(ndup_f)>::get(dup_f,ndup_f);
  auto a_f = v_f.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  F_FLOAT delij[3], delik[3], deljk[3], tmpvec[3];
  F_FLOAT dBOp_i[3], dBOp_k[3], dln_BOp_pi[3], dln_BOp_pi2[3];

  const int i = d_ilist[ii];
  const X_FLOAT xtmp = x(i,0);
  const X_FLOAT ytmp = x(i,1);
  const X_FLOAT ztmp = x(i,2);
  const tagint itag = tag(i);
  const int jknum = d_bo_num[i];

  F_FLOAT CdDelta_i = d_CdDelta[i];
  F_FLOAT fitmp[3];
  for (int j = 0; j < 3; j++) fitmp[j] = 0.0;

  for (int j_index = 0; j_index < jknum; j_index++) {
    int j = d_bo_list(i, j_index);
    j &= NEIGHMASK;
    const tagint jtag = tag(j);

    if (itag > jtag) {
      if ((itag+jtag) % 2 == 0) continue;
    } else if (itag < jtag) {
      if ((itag+jtag) % 2 == 1) continue;
    } else {
      if (x(j,2)  < ztmp) continue;
      if (x(j,2) == ztmp && x(j,1)  < ytmp) continue;
      if (x(j,2) == ztmp && x(j,1) == ytmp && x(j,0) < xtmp) continue;
    }

    F_FLOAT CdDelta_j = d_CdDelta[j];

    delij[0] = x(j,0) - xtmp;
    delij[1] = x(j,1) - ytmp;
    delij[2] = x(j,2) - ztmp;

    const int knum = d_bo_num[j];

    F_FLOAT coef_C1dbo, coef_C2dbo, coef_C3dbo, coef_C1dbopi, coef_C2dbopi, coef_C3dbopi, coef_C4dbopi;
    F_FLOAT coef_C1dbopi2, coef_C2dbopi2, coef_C3dbopi2, coef_C4dbopi2, coef_C1dDelta, coef_C2dDelta, coef_C3dDelta;

    coef_C1dbo = coef_C2dbo = coef_C3dbo = 0.0;
    coef_C1dbopi = coef_C2dbopi = coef_C3dbopi = coef_C4dbopi = 0.0;
    coef_C1dbopi2 = coef_C2dbopi2 = coef_C3dbopi2 = coef_C4dbopi2 = 0.0;
    coef_C1dDelta = coef_C2dDelta = coef_C3dDelta = 0.0;

    // total forces on i, j, k (nlocal + nghost, from Add_dBond_to_Forces))
    const F_FLOAT Cdbo_ij = d_Cdbo(i,j_index);
    coef_C1dbo = d_C1dbo(i,j_index) * (Cdbo_ij);
    coef_C2dbo = d_C2dbo(i,j_index) * (Cdbo_ij);
    coef_C3dbo = d_C3dbo(i,j_index) * (Cdbo_ij);

    const F_FLOAT Cdbopi_ij = d_Cdbopi(i,j_index);
    coef_C1dbopi = d_C1dbopi(i,j_index) * (Cdbopi_ij);
    coef_C2dbopi = d_C2dbopi(i,j_index) * (Cdbopi_ij);
    coef_C3dbopi = d_C3dbopi(i,j_index) * (Cdbopi_ij);
    coef_C4dbopi = d_C4dbopi(i,j_index) * (Cdbopi_ij);

    const F_FLOAT Cdbopi2_ij = d_Cdbopi2(i,j_index);
    coef_C1dbopi2 = d_C1dbopi2(i,j_index) * (Cdbopi2_ij);
    coef_C2dbopi2 = d_C2dbopi2(i,j_index) * (Cdbopi2_ij);
    coef_C3dbopi2 = d_C3dbopi2(i,j_index) * (Cdbopi2_ij);
    coef_C4dbopi2 = d_C4dbopi2(i,j_index) * (Cdbopi2_ij);

    const F_FLOAT coeff_CdDelta_ij = CdDelta_i + CdDelta_j;
    coef_C1dDelta = d_C1dbo(i,j_index) * (coeff_CdDelta_ij);
    coef_C2dDelta = d_C2dbo(i,j_index) * (coeff_CdDelta_ij);
    coef_C3dDelta = d_C3dbo(i,j_index) * (coeff_CdDelta_ij);

    F_FLOAT temp[3];

    F_FLOAT d_dln_BOp_pi_local = d_dln_BOp_pi(i,j_index);
    for (int d = 0; d < 3; d++) dln_BOp_pi[d] = d_dln_BOp_pi_local * delij[d];

    F_FLOAT d_dln_BOp_pi2_local = d_dln_BOp_pi2(i,j_index);
    for (int d = 0; d < 3; d++) dln_BOp_pi2[d] = d_dln_BOp_pi2_local * delij[d];

    F_FLOAT d_dBOp_local = d_dBOp(i,j_index);
    for (int d = 0; d < 3; d++) dBOp_i[d] = d_dBOp_local * delij[d];

    // forces on i
    for (int d = 0; d < 3; d++) temp[d] =  coef_C1dbo * dBOp_i[d];
    for (int d = 0; d < 3; d++) temp[d] += coef_C2dbo * d_dDeltap_self(i,d);
    for (int d = 0; d < 3; d++) temp[d] += coef_C1dDelta * dBOp_i[d];
    for (int d = 0; d < 3; d++) temp[d] += coef_C2dDelta * d_dDeltap_self(i,d);
    for (int d = 0; d < 3; d++) temp[d] += coef_C1dbopi * dln_BOp_pi[d];
    for (int d = 0; d < 3; d++) temp[d] += coef_C2dbopi * dBOp_i[d];
    for (int d = 0; d < 3; d++) temp[d] += coef_C3dbopi * d_dDeltap_self(i,d);
    for (int d = 0; d < 3; d++) temp[d] += coef_C1dbopi2 * dln_BOp_pi2[d];
    for (int d = 0; d < 3; d++) temp[d] += coef_C2dbopi2 * dBOp_i[d];
    for (int d = 0; d < 3; d++) temp[d] += coef_C3dbopi2 * d_dDeltap_self(i,d);

    if (VFLAG && vflag_either) this->template v_tally<NEIGHFLAG>(ev,i,temp,delij);

    fitmp[0] -= temp[0];
    fitmp[1] -= temp[1];
    fitmp[2] -= temp[2];

    // forces on j
    for (int d = 0; d < 3; d++) temp[d] = -coef_C1dbo * dBOp_i[d];
    for (int d = 0; d < 3; d++) temp[d] += coef_C3dbo * d_dDeltap_self(j,d);
    for (int d = 0; d < 3; d++) temp[d] -= coef_C1dDelta * dBOp_i[d];
    for (int d = 0; d < 3; d++) temp[d] += coef_C3dDelta * d_dDeltap_self(j,d);
    for (int d = 0; d < 3; d++) temp[d] -= coef_C1dbopi * dln_BOp_pi[d];
    for (int d = 0; d < 3; d++) temp[d] -= coef_C2dbopi * dBOp_i[d];
    for (int d = 0; d < 3; d++) temp[d] += coef_C4dbopi * d_dDeltap_self(j,d);
    for (int d = 0; d < 3; d++) temp[d] -= coef_C1dbopi2 * dln_BOp_pi2[d];
    for (int d = 0; d < 3; d++) temp[d] -= coef_C2dbopi2 * dBOp_i[d];
    for (int d = 0; d < 3; d++) temp[d] += coef_C4dbopi2 * d_dDeltap_self(j,d);

    a_f(j,0) -= temp[0];
    a_f(j,1) -= temp[1];
    a_f(j,2) -= temp[2];

    if (VFLAG && vflag_either) {
      for (int d = 0; d < 3; d++) tmpvec[d] = -delij[d];
      this->template v_tally<NEIGHFLAG>(ev,j,temp,tmpvec);
    }

    // forces on k: i neighbor
    for (int k_index = 0; k_index < jknum; k_index++) {
      int k = d_bo_list(i, k_index);
      k &= NEIGHMASK;

      delik[0] = x(k,0) - xtmp;
      delik[1] = x(k,1) - ytmp;
      delik[2] = x(k,2) - ztmp;

      d_dBOp_local = d_dBOp(i,k_index);
      for (int d = 0; d < 3; d++) dBOp_k[d] = d_dBOp_local * delik[d];

      const F_FLOAT coef_all = -coef_C2dbo - coef_C2dDelta - coef_C3dbopi - coef_C3dbopi2;
      for (int d = 0; d < 3; d++) temp[d] = coef_all * dBOp_k[d];

      a_f(k,0) -= temp[0];
      a_f(k,1) -= temp[1];
      a_f(k,2) -= temp[2];

      if (VFLAG && vflag_either) {
        delik[0] = x(k,0) - xtmp;
        delik[1] = x(k,1) - ytmp;
        delik[2] = x(k,2) - ztmp;
        for (int d = 0; d < 3; d++) tmpvec[d] = x(j,d) - x(k,d) - delik[d];
        this->template v_tally<NEIGHFLAG>(ev,k,temp,tmpvec);
      }

    }

    // forces on k: j neighbor
    for (int k_index = 0; k_index < knum; k_index++) {
      int k = d_bo_list(j, k_index);
      k &= NEIGHMASK;

      for (int d = 0; d < 3; d++) deljk[d] = x(k,d) - x(j,d);

      d_dBOp_local = d_dBOp(j,k_index);
      for (int d = 0; d < 3; d++) dBOp_k[d] = d_dBOp_local * deljk[d];
      const F_FLOAT coef_all = -coef_C3dbo - coef_C3dDelta - coef_C4dbopi - coef_C4dbopi2;
      for (int d = 0; d < 3; d++) temp[d] = coef_all * dBOp_k[d];

      a_f(k,0) -= temp[0];
      a_f(k,1) -= temp[1];
      a_f(k,2) -= temp[2];

      if (VFLAG && vflag_either) {
        for (int d = 0; d < 3; d++) deljk[d] = x(k,d) - x(j,d);
        for (int d = 0; d < 3; d++) tmpvec[d] = x(i,d) - x(k,d) - deljk[d];
        this->template v_tally<NEIGHFLAG>(ev,k,temp,tmpvec);
      }
    }
  }
  for (int d = 0; d < 3; d++) a_f(i,d) += fitmp[d];
}

template<class DeviceType>
template<int NEIGHFLAG, int VFLAG>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::operator()(TagPairReaxComputeBond2<NEIGHFLAG,VFLAG>, const int &ii) const {
  EV_FLOAT_REAX ev;
  this->template operator()<NEIGHFLAG,VFLAG>(TagPairReaxComputeBond2<NEIGHFLAG,VFLAG>(), ii, ev);
}


#endif
