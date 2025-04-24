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

#ifndef LMP_PAIR_REAXFF_ANGULAR_KOKKOS_HPP
#define LMP_PAIR_REAXFF_ANGULAR_KOKKOS_HPP

#include "pair_reaxff_kokkos.h"

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<bool POPULATE>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::operator()(TagPairReaxCountAngularTorsion<POPULATE>, const int &ii) const {

  const int i = d_ilist[ii];
  const int itype = type(i);

  const int jnum = d_bo_num[i];

  if (POPULATE) {
    // Computes and stores SBO2, CSBO2, dSBO1, dSBO2
    compute_angular_sbo(i, itype, jnum);
  }

  // Angular

  // Count buffer size for `i`
  int location_angular = 0; // dummy declaration
  int count_angular = preprocess_angular<false>(i, itype, jnum, location_angular);
  location_angular = Kokkos::atomic_fetch_add(&d_count_angular_torsion(0), count_angular);

  if (POPULATE) {
    // Fill buffer for `i`
    preprocess_angular<true>(i, itype, jnum, location_angular);
  }

  // Torsion

  const tagint itag = tag(i);
  const X_FLOAT xtmp = x(i,0);
  const X_FLOAT ytmp = x(i,1);
  const X_FLOAT ztmp = x(i,2);

  // Count buffer size for `i`
  int location_torsion = 0; // dummy declaration
  int count_torsion = preprocess_torsion<false>(i, itype, itag, xtmp, ytmp, ztmp, jnum, location_torsion);
  location_torsion = Kokkos::atomic_fetch_add(&d_count_angular_torsion(1), count_torsion);

  if (POPULATE) {
    // Fill buffer for `i`
    preprocess_torsion<true>(i, itype, itag, xtmp, ytmp, ztmp, jnum, location_torsion);
  }

}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::compute_angular_sbo(int i, int itype, int jnum) const {

  F_FLOAT SBO2, CSBO2, dSBO1, dSBO2;

  const F_FLOAT p_val8 = gp[33];
  const F_FLOAT p_val9 = gp[16];

  F_FLOAT SBOp = 0.0;
  F_FLOAT prod_SBO = 1.0;

  for (int j_index = 0; j_index < jnum; j_index++) {
    const F_FLOAT bo_ij = d_BO(i,j_index);

    SBOp += (d_BO_pi(i,j_index) + d_BO_pi2(i,j_index));
    F_FLOAT temp = SQR(bo_ij);
    temp *= temp;
    temp *= temp;
    prod_SBO *= exp(-temp);
  }

  F_FLOAT vlpadj;

  const F_FLOAT Delta_e = d_total_bo[i] - paramssing(itype).valency_e;
  const F_FLOAT vlpex = Delta_e - 2.0 * (int)(Delta_e/2.0);
  const F_FLOAT explp1 = exp(-gp[15] * SQR(2.0 + vlpex));
  const F_FLOAT nlp = explp1 - (int)(Delta_e / 2.0);
  if (vlpex >= 0.0) {
    vlpadj = 0.0;
    dSBO2 = prod_SBO - 1.0;
  } else {
    vlpadj = nlp;
    dSBO2 = (prod_SBO - 1.0) * (1.0 - p_val8 * d_dDelta_lp[i]);
  }

  const F_FLOAT SBO = SBOp + (1.0 - prod_SBO) * (-d_Delta_boc[i] - p_val8 * vlpadj);
  dSBO1 = -8.0 * prod_SBO * (d_Delta_boc[i] + p_val8 * vlpadj);

  if (SBO <= 0.0) {
    SBO2 = 0.0;
    CSBO2 = 0.0;
  } else if (SBO > 0.0 && SBO <= 1.0) {
    SBO2 = pow(SBO, p_val9);
    CSBO2 = p_val9 * pow(SBO, p_val9 - 1.0);
  } else if (SBO > 1.0 && SBO < 2.0) {
    SBO2 = 2.0 - pow(2.0-SBO, p_val9);
    CSBO2 = p_val9 * pow(2.0 - SBO, p_val9 - 1.0);
  } else {
    SBO2 = 2.0;
    CSBO2 = 0.0;
  }

  d_angular_intermediates(i,0) = SBO2;
  d_angular_intermediates(i,1) = CSBO2;
  d_angular_intermediates(i,2) = dSBO1;
  d_angular_intermediates(i,3) = dSBO2;

}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<bool POPULATE>
KOKKOS_INLINE_FUNCTION
int PairReaxFFKokkos<DeviceType>::preprocess_angular(int i, int itype, int jnum, int location_angular) const {

  int count_angular = 0;

  for (int j_index = 0; j_index < jnum; j_index++) {
    int j = d_bo_list(i, j_index);
    j &= NEIGHMASK;
    const F_FLOAT bo_ij = d_BO(i,j_index);

    if (bo_ij <= thb_cut) continue;
    if (i >= nlocal && j >= nlocal) continue;

    const int jtype = type(j);

    for (int k_index = j_index + 1; k_index < jnum; k_index++) {
    //for (int kk = j_start; kk < j_end; kk++) {
      int k = d_bo_list(i, k_index);
      k &= NEIGHMASK;
      if (k == j) continue;

      const F_FLOAT bo_ik = d_BO(i,k_index);

      if (bo_ij <= thb_cut || bo_ik <= thb_cut || bo_ij * bo_ik <= thb_cutsq) continue;

      const int ktype = type(k);

      F_FLOAT p_val1 = paramsthbp(jtype,itype,ktype).p_val1;

      if (fabs(p_val1) <= 0.001) continue;

      if (POPULATE) {
        reax_int4 pack;

        // First pack stores i, j, k, and j_start
        pack.i0 = i;
        pack.i1 = j;
        pack.i2 = k;
        pack.i3 = jnum;
        d_angular_pack(location_angular, 0) = pack;

        // Second pack stores j_index and k_index
        // i0 is unused because there's no i_index
        pack.i1 = j_index;
        pack.i2 = k_index;
        // i3 is unused
        d_angular_pack(location_angular, 1) = pack;

        location_angular++;
      } else {
        count_angular++;
      }
    }
  }

  return count_angular;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<bool POPULATE>
KOKKOS_INLINE_FUNCTION
int PairReaxFFKokkos<DeviceType>::preprocess_torsion(int i, int /*itype*/, tagint itag,
  F_FLOAT xtmp, F_FLOAT ytmp, F_FLOAT ztmp, int jknum, int location_torsion) const {

  // in reaxff_torsion_angles: j = i, k = j, i = k;

  int count_torsion = 0;

  for (int j_index = 0; j_index < jknum; j_index++) {
    int j = d_bo_list(i, j_index);
    j &= NEIGHMASK;
    const tagint jtag = tag(j);

    // skip half of the interactions
    if (itag > jtag) {
      if ((itag+jtag) % 2 == 0) continue;
    } else if (itag < jtag) {
      if ((itag+jtag) % 2 == 1) continue;
    } else {
      if (x(j,2)  < ztmp) continue;
      if (x(j,2) == ztmp && x(j,1)  < ytmp) continue;
      if (x(j,2) == ztmp && x(j,1) == ytmp && x(j,0) < xtmp) continue;
    }

    const F_FLOAT bo_ij = d_BO(i,j_index);
    if (bo_ij < thb_cut) continue;

    const int lnum = d_bo_num[j];

    for (int k_index = 0; k_index < jknum; k_index++) {
      int k = d_bo_list(i, k_index);
      k &= NEIGHMASK;
      if (k == j) continue;

      const F_FLOAT bo_ik = d_BO(i,k_index);
      if (bo_ik < thb_cut) continue;

      for (int l_index = 0; l_index < lnum; l_index++) {
        int l = d_bo_list(j, l_index);
        l &= NEIGHMASK;
        if (l == i) continue;

        const F_FLOAT bo_jl = d_BO(j,l_index);
        if (l == k || bo_jl < thb_cut || bo_ij*bo_ik*bo_jl < thb_cut) continue;

        if (POPULATE) {
          reax_int4 pack;

          pack.i0 = i;
          pack.i1 = j;
          pack.i2 = k;
          pack.i3 = l;
          d_torsion_pack(location_torsion, 0) = pack;

          pack.i0 = 0; // no i_index
          pack.i1 = j_index;
          pack.i2 = k_index;
          pack.i3 = l_index;
          d_torsion_pack(location_torsion, 1) = pack;

          location_torsion++;
        } else {
          count_torsion++;
        }
      }
    }
  }

  return count_torsion;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::operator()(TagPairReaxComputeAngularPreprocessed<NEIGHFLAG,EVFLAG>, const int &apack, EV_FLOAT_REAX& ev) const {

  auto v_f = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_f),decltype(ndup_f)>::get(dup_f,ndup_f);
  auto a_f = v_f.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();
  Kokkos::View<F_FLOAT**, typename decltype(d_Cdbo)::array_layout,KKDeviceType,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value>> a_Cdbo = d_Cdbo;
  Kokkos::View<F_FLOAT**, typename decltype(d_Cdbopi)::array_layout,KKDeviceType,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value>> a_Cdbopi = d_Cdbopi;
  Kokkos::View<F_FLOAT**, typename decltype(d_Cdbopi2)::array_layout,KKDeviceType,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value>> a_Cdbopi2 = d_Cdbopi2;

  auto v_CdDelta = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_CdDelta),decltype(ndup_CdDelta)>::get(dup_CdDelta,ndup_CdDelta);
  auto a_CdDelta = v_CdDelta.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  F_FLOAT temp, temp_bo_jt, pBOjt7;
  F_FLOAT p_val1, p_val2, p_val3, p_val4, p_val5;
  F_FLOAT p_val6, p_val7, p_val10;
  F_FLOAT p_pen1, p_pen2, p_pen3, p_pen4;
  F_FLOAT p_coa1, p_coa2, p_coa3, p_coa4;
  F_FLOAT trm8, expval6, expval7, expval2theta, expval12theta, exp3ij, exp3jk;
  F_FLOAT exp_pen2ij, exp_pen2jk, exp_pen3, exp_pen4, trm_pen34, exp_coa2;
  F_FLOAT dSBO1, dSBO2, SBO2, CSBO2;
  F_FLOAT CEval1, CEval2, CEval3, CEval4, CEval5, CEval6, CEval7, CEval8;
  F_FLOAT CEpen1, CEpen2, CEpen3;
  F_FLOAT e_ang, e_coa, e_pen;
  F_FLOAT CEcoa1, CEcoa2, CEcoa3, CEcoa4, CEcoa5;
  F_FLOAT Cf7ij, Cf7jk, Cf8j, Cf9j;
  F_FLOAT f7_ij, f7_jk, f8_Dj, f9_Dj;
  F_FLOAT Ctheta_0, theta_0, theta_00, theta, cos_theta, sin_theta;
  F_FLOAT BOA_ij, BOA_ik, rij, bo_ij, bo_ik;
  F_FLOAT dcos_theta_di[3], dcos_theta_dj[3], dcos_theta_dk[3];
  F_FLOAT eng_tmp, fi_tmp[3], fj_tmp[3], fk_tmp[3];
  F_FLOAT delij[3], delik[3], delji[3], delki[3];

  p_val6 = gp[14];
  p_val10 = gp[17];

  p_pen2 = gp[19];
  p_pen3 = gp[20];
  p_pen4 = gp[21];

  p_coa2 = gp[2];
  p_coa3 = gp[38];
  p_coa4 = gp[30];

  reax_int4 pack = d_angular_pack(apack,0);
  const int i = pack.i0;
  const int j = pack.i1;
  const int k = pack.i2;
  const int jnum = pack.i3;

  pack = d_angular_pack(apack, 1);
  // i0 is unused
  const int j_index = pack.i1;
  const int k_index = pack.i2;
  // i3 is unused

  const int itype = type(i);
  const X_FLOAT xtmp = x(i,0);
  const X_FLOAT ytmp = x(i,1);
  const X_FLOAT ztmp = x(i,2);

  p_val3 = paramssing(itype).p_val3;
  p_val5 = paramssing(itype).p_val5;

  const F_FLOAT Delta_val = d_total_bo[i] - paramssing(itype).valency_val;

  SBO2 = d_angular_intermediates(i, 0);
  CSBO2 = d_angular_intermediates(i, 1);
  dSBO1 = d_angular_intermediates(i, 2);
  dSBO2 = d_angular_intermediates(i, 3);

  expval6 = exp(p_val6 * d_Delta_boc[i]);

  F_FLOAT CdDelta_i = 0.0;
  F_FLOAT fitmp[3],fjtmp[3];
  for (int j = 0; j < 3; j++) fitmp[j] = 0.0;

  delij[0] = x(j,0) - xtmp;
  delij[1] = x(j,1) - ytmp;
  delij[2] = x(j,2) - ztmp;
  const F_FLOAT rsqij = delij[0]*delij[0] + delij[1]*delij[1] + delij[2]*delij[2];
  rij = sqrt(rsqij);
  bo_ij = d_BO(i,j_index);

  BOA_ij = bo_ij - thb_cut;

  const int jtype = type(j);

  F_FLOAT CdDelta_j = 0.0;
  for (int k = 0; k < 3; k++) fjtmp[k] = 0.0;

  delik[0] = x(k,0) - xtmp;
  delik[1] = x(k,1) - ytmp;
  delik[2] = x(k,2) - ztmp;
  const F_FLOAT rsqik = delik[0]*delik[0] + delik[1]*delik[1] + delik[2]*delik[2];
  const F_FLOAT rik = sqrt(rsqik);
  bo_ik = d_BO(i,k_index);
  BOA_ik   = bo_ik - thb_cut;

  const int ktype = type(k);

  // theta and derivatives

  cos_theta = (delij[0]*delik[0]+delij[1]*delik[1]+delij[2]*delik[2])/(rij*rik);
  if (cos_theta > 1.0) cos_theta  = 1.0;
  if (cos_theta < -1.0) cos_theta  = -1.0;
  theta = acos(cos_theta);

  const F_FLOAT inv_dists = 1.0 / (rij * rik);
  const F_FLOAT Cdot_inv3 = cos_theta * inv_dists * inv_dists;

  for (int t = 0; t < 3; t++) {
    dcos_theta_di[t] = -(delik[t] + delij[t]) * inv_dists + Cdot_inv3 * (rsqik * delij[t] + rsqij * delik[t]);
    dcos_theta_dj[t] = delik[t] * inv_dists - Cdot_inv3 * rsqik * delij[t];
    dcos_theta_dk[t] = delij[t] * inv_dists - Cdot_inv3 * rsqij * delik[t];
  }

  sin_theta = sin(theta);
  if (sin_theta < 1.0e-5) sin_theta = 1.0e-5;
  p_val1 = paramsthbp(jtype,itype,ktype).p_val1;

  // ANGLE ENERGY

  p_val1 = paramsthbp(jtype,itype,ktype).p_val1;
  p_val2 = paramsthbp(jtype,itype,ktype).p_val2;
  p_val4 = paramsthbp(jtype,itype,ktype).p_val4;
  p_val7 = paramsthbp(jtype,itype,ktype).p_val7;
  theta_00 = paramsthbp(jtype,itype,ktype).theta_00;

  exp3ij = exp(-p_val3 * pow(BOA_ij, p_val4));
  f7_ij = 1.0 - exp3ij;
  Cf7ij = p_val3 * p_val4 * pow(BOA_ij, p_val4 - 1.0) * exp3ij;
  exp3jk = exp(-p_val3 * pow(BOA_ik, p_val4));
  f7_jk = 1.0 - exp3jk;
  Cf7jk = p_val3 * p_val4 * pow(BOA_ik, p_val4 - 1.0) * exp3jk;
  expval7 = exp(-p_val7 * d_Delta_boc[i]);
  trm8 = 1.0 + expval6 + expval7;
  f8_Dj = p_val5 - ((p_val5 - 1.0) * (2.0 + expval6) / trm8);
  Cf8j = ((1.0 - p_val5) / (trm8*trm8)) *
   (p_val6 * expval6 * trm8 - (2.0 + expval6) * (p_val6*expval6 - p_val7*expval7));
  theta_0 = 180.0 - theta_00 * (1.0 - exp(-p_val10 * (2.0 - SBO2)));
  theta_0 = theta_0*constPI/180.0;

  expval2theta  = exp(-p_val2 * (theta_0-theta)*(theta_0-theta));
  if (p_val1 >= 0)
    expval12theta = p_val1 * (1.0 - expval2theta);
  else // To avoid linear Me-H-Me angles (6/6/06)
    expval12theta = p_val1 * -expval2theta;

  CEval1 = Cf7ij * f7_jk * f8_Dj * expval12theta;
  CEval2 = Cf7jk * f7_ij * f8_Dj * expval12theta;
  CEval3 = Cf8j  * f7_ij * f7_jk * expval12theta;
  CEval4 = -2.0 * p_val1 * p_val2 * f7_ij * f7_jk * f8_Dj * expval2theta * (theta_0 - theta);
  Ctheta_0 = p_val10 * theta_00*constPI/180.0 * exp(-p_val10 * (2.0 - SBO2));
  CEval5 = -CEval4 * Ctheta_0 * CSBO2;
  CEval6 = CEval5 * dSBO1;
  CEval7 = CEval5 * dSBO2;
  CEval8 = -CEval4 / sin_theta;

  e_ang = f7_ij * f7_jk * f8_Dj * expval12theta;
  if (eflag) ev.ereax[3] += e_ang;

  // Penalty energy

  p_pen1 = paramsthbp(jtype,itype,ktype).p_pen1;

  exp_pen2ij = exp(-p_pen2 * (BOA_ij - 2.0)*(BOA_ij - 2.0));
  exp_pen2jk = exp(-p_pen2 * (BOA_ik - 2.0)*(BOA_ik - 2.0));
  exp_pen3 = exp(-p_pen3 * d_Delta[i]);
  exp_pen4 = exp(p_pen4 * d_Delta[i]);
  trm_pen34 = 1.0 + exp_pen3 + exp_pen4;
  f9_Dj = (2.0 + exp_pen3) / trm_pen34;
  Cf9j = (-p_pen3 * exp_pen3 * trm_pen34 - (2.0 + exp_pen3) *
   (-p_pen3 * exp_pen3 + p_pen4 * exp_pen4))/(trm_pen34*trm_pen34);

  e_pen = p_pen1 * f9_Dj * exp_pen2ij * exp_pen2jk;
  if (eflag) ev.ereax[4] += e_pen;

  CEpen1 = e_pen * Cf9j / f9_Dj;
  temp   = -2.0 * p_pen2 * e_pen;
  CEpen2 = temp * (BOA_ij - 2.0);
  CEpen3 = temp * (BOA_ik - 2.0);

  // ConjAngle energy

  p_coa1 = paramsthbp(jtype,itype,ktype).p_coa1;
  exp_coa2 = exp(p_coa2 * Delta_val);
  e_coa = p_coa1 / (1. + exp_coa2) *
          exp(-p_coa3 * SQR(d_total_bo[j]-BOA_ij)) *
          exp(-p_coa3 * SQR(d_total_bo[k]-BOA_ik)) *
          exp(-p_coa4 * SQR(BOA_ij - 1.5)) *
          exp(-p_coa4 * SQR(BOA_ik - 1.5));

  CEcoa1 = -2 * p_coa4 * (BOA_ij - 1.5) * e_coa;
  CEcoa2 = -2 * p_coa4 * (BOA_ik - 1.5) * e_coa;
  CEcoa3 = -p_coa2 * exp_coa2 * e_coa / (1 + exp_coa2);
  CEcoa4 = -2 * p_coa3 * (d_total_bo[j]-BOA_ij) * e_coa;
  CEcoa5 = -2 * p_coa3 * (d_total_bo[k]-BOA_ik) * e_coa;

  if (eflag) ev.ereax[5] += e_coa;

  // Forces

  a_Cdbo(i,j_index) += (CEval1 + CEpen2 + (CEcoa1 - CEcoa4));
  a_Cdbo(i,k_index) += (CEval2 + CEpen3 + (CEcoa2 - CEcoa5));

  CdDelta_i += ((CEval3 + CEval7) + CEpen1 + CEcoa3);
  CdDelta_j += CEcoa4;
  a_CdDelta[k] += CEcoa5;

  for (int l_index = 0; l_index < jnum; l_index++) {
    temp_bo_jt = d_BO(i,l_index);
    temp = temp_bo_jt * temp_bo_jt * temp_bo_jt;
    pBOjt7 = temp * temp * temp_bo_jt;

    a_Cdbo(i,l_index) += (CEval6 * pBOjt7);
    a_Cdbopi(i,l_index) += CEval5;
    a_Cdbopi2(i,l_index) += CEval5;
  }

  for (int d = 0; d < 3; d++) fi_tmp[d] = CEval8 * dcos_theta_di[d];
  for (int d = 0; d < 3; d++) fj_tmp[d] = CEval8 * dcos_theta_dj[d];
  for (int d = 0; d < 3; d++) fk_tmp[d] = CEval8 * dcos_theta_dk[d];
  for (int d = 0; d < 3; d++) fitmp[d] -= fi_tmp[d];
  for (int d = 0; d < 3; d++) fjtmp[d] -= fj_tmp[d];
  for (int d = 0; d < 3; d++) a_f(k,d) -= fk_tmp[d];

  // energy/virial tally
  if (EVFLAG) {
    eng_tmp = e_ang + e_pen + e_coa;
    //if (eflag_atom) this->template ev_tally<NEIGHFLAG>(ev,i,j,eng_tmp,0.0,0.0,0.0,0.0);
    for (int d = 0; d < 3; d++) delki[d] = -1.0 * delik[d];
    for (int d = 0; d < 3; d++) delji[d] = -1.0 * delij[d];
    if (eflag_atom) this->template e_tally<NEIGHFLAG>(ev,i,j,eng_tmp);
    if (vflag_either) this->template v_tally3<NEIGHFLAG>(ev,i,j,k,fj_tmp,fk_tmp,delji,delki);
  }

  a_CdDelta[j] += CdDelta_j;
  for (int d = 0; d < 3; d++) a_f(j,d) += fjtmp[d];
  a_CdDelta[i] += CdDelta_i;
  for (int d = 0; d < 3; d++) a_f(i,d) += fitmp[d];
}

template<class DeviceType>
template<int NEIGHFLAG, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::operator()(TagPairReaxComputeAngularPreprocessed<NEIGHFLAG,EVFLAG>, const int &apack) const {
  EV_FLOAT_REAX ev;
  this->template operator()<NEIGHFLAG,EVFLAG>(TagPairReaxComputeAngularPreprocessed<NEIGHFLAG,EVFLAG>(), apack, ev);

}





#endif
