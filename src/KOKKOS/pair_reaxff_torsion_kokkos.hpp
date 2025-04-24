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

#ifndef LMP_PAIR_REAXFF_TORSION_KOKKOS_HPP
#define LMP_PAIR_REAXFF_TORSION_KOKKOS_HPP

#include "pair_reaxff_kokkos.h"

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::operator()(TagPairReaxComputeTorsionPreprocessed<NEIGHFLAG,EVFLAG>, const int &tpack, EV_FLOAT_REAX& ev) const {

  auto v_f = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_f),decltype(ndup_f)>::get(dup_f,ndup_f);
  auto a_f = v_f.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  auto v_CdDelta = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_CdDelta),decltype(ndup_CdDelta)>::get(dup_CdDelta,ndup_CdDelta);
  auto a_CdDelta = v_CdDelta.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();
  Kokkos::View<F_FLOAT**, typename decltype(d_Cdbo)::array_layout,KKDeviceType,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value>> a_Cdbo = d_Cdbo;
  Kokkos::View<F_FLOAT**, typename decltype(d_Cdbopi)::array_layout,KKDeviceType,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value>> a_Cdbopi = d_Cdbopi;
  //auto a_Cdbo = dup_Cdbo.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  // in reaxff_torsion_angles: j = i, k = j, i = k;

  F_FLOAT Delta_i, Delta_j, bo_ij, bo_ik, bo_jl, BOA_ij, BOA_ik, BOA_jl;
  F_FLOAT p_tor1, p_cot1, V1, V2, V3;
  F_FLOAT exp_tor2_ij, exp_tor2_ik, exp_tor2_jl, exp_tor1, exp_tor3_DiDj, exp_tor4_DiDj, exp_tor34_inv;
  F_FLOAT exp_cot2_ij, exp_cot2_ik, exp_cot2_jl, fn10, f11_DiDj, dfn11, fn12;
  F_FLOAT theta_ijk, theta_jil, sin_ijk, sin_jil, cos_ijk, cos_jil, tan_ijk_i, tan_jil_i;
  F_FLOAT cos_omega, cos2omega, cos3omega;
  F_FLOAT CV, cmn, CEtors1, CEtors2, CEtors3, CEtors4;
  F_FLOAT CEtors5, CEtors6, CEtors7, CEtors8, CEtors9;
  F_FLOAT Cconj, CEconj1, CEconj2, CEconj3, CEconj4, CEconj5, CEconj6;
  F_FLOAT e_tor, e_con, eng_tmp;

  F_FLOAT delij[3], delik[3], deljl[3], dellk[3], delil[3], delkl[3];
  F_FLOAT fi_tmp[3], fj_tmp[3], fk_tmp[3];
  F_FLOAT dcos_ijk_di[3], dcos_ijk_dj[3], dcos_ijk_dk[3], dcos_jil_di[3], dcos_jil_dj[3], dcos_jil_dk[3];

  F_FLOAT p_tor2 = gp[23];
  F_FLOAT p_tor3 = gp[24];
  F_FLOAT p_tor4 = gp[25];
  F_FLOAT p_cot2 = gp[27];

  reax_int4 pack = d_torsion_pack(tpack,0);
  const int i = pack.i0;
  const int j = pack.i1;
  const int k = pack.i2;
  const int l = pack.i3;

  pack = d_torsion_pack(tpack, 1);
  //const int i = pack.i0;
  const int j_index = pack.i1;
  const int k_index = pack.i2;
  const int l_index = pack.i3;

  const int itype = type(i);
  const X_FLOAT xtmp = x(i,0);
  const X_FLOAT ytmp = x(i,1);
  const X_FLOAT ztmp = x(i,2);
  Delta_i = d_Delta_boc[i];

  const int jtype = type(j);

  bo_ij = d_BO(i,j_index);

  delij[0] = x(j,0) - xtmp;
  delij[1] = x(j,1) - ytmp;
  delij[2] = x(j,2) - ztmp;
  const F_FLOAT rsqij = delij[0]*delij[0] + delij[1]*delij[1] + delij[2]*delij[2];
  const F_FLOAT rij = sqrt(rsqij);

  BOA_ij = bo_ij - thb_cut;
  Delta_j = d_Delta_boc[j];
  exp_tor2_ij = exp(-p_tor2 * BOA_ij);
  exp_cot2_ij = exp(-p_cot2 * SQR(BOA_ij - 1.5));
  exp_tor3_DiDj = exp(-p_tor3 * (Delta_i + Delta_j));
  exp_tor4_DiDj = exp(p_tor4  * (Delta_i + Delta_j));
  exp_tor34_inv = 1.0 / (1.0 + exp_tor3_DiDj + exp_tor4_DiDj);
  f11_DiDj = (2.0 + exp_tor3_DiDj) * exp_tor34_inv;

  const int ktype = type(k);

  bo_ik = d_BO(i,k_index);

  BOA_ik = bo_ik - thb_cut;
  for (int d = 0; d < 3; d ++) delik[d] = x(k,d) - x(i,d);
  const F_FLOAT rsqik = delik[0]*delik[0] + delik[1]*delik[1] + delik[2]*delik[2];
  const F_FLOAT rik = sqrt(rsqik);

  cos_ijk = (delij[0]*delik[0]+delij[1]*delik[1]+delij[2]*delik[2])/(rij*rik);
  if (cos_ijk > 1.0) cos_ijk  = 1.0;
  if (cos_ijk < -1.0) cos_ijk  = -1.0;
  theta_ijk = acos(cos_ijk);

  // dcos_ijk
  const F_FLOAT inv_dists = 1.0 / (rij * rik);
  const F_FLOAT cos_ijk_tmp = cos_ijk / ((rij*rik)*(rij*rik));

  for (int d = 0; d < 3; d++) {
    dcos_ijk_di[d] = -(delik[d] + delij[d]) * inv_dists + cos_ijk_tmp * (rsqik * delij[d] + rsqij * delik[d]);
    dcos_ijk_dj[d] = delik[d] * inv_dists - cos_ijk_tmp * rsqik * delij[d];
    dcos_ijk_dk[d] = delij[d] * inv_dists - cos_ijk_tmp * rsqij * delik[d];
  }

  sin_ijk = sin(theta_ijk);
  if (sin_ijk >= 0 && sin_ijk <= MIN_SINE)
    tan_ijk_i = cos_ijk / MIN_SINE;
  else if (sin_ijk <= 0 && sin_ijk >= -MIN_SINE)
    tan_ijk_i = -cos_ijk / MIN_SINE;
  else tan_ijk_i = cos_ijk / sin_ijk;

  exp_tor2_ik = exp(-p_tor2 * BOA_ik);
  exp_cot2_ik = exp(-p_cot2 * SQR(BOA_ik -1.5));

  const int ltype = type(l);

  bo_jl = d_BO(j,l_index);

  for (int d = 0; d < 3; d ++) deljl[d] = x(l,d) - x(j,d);
  const F_FLOAT rsqjl = deljl[0]*deljl[0] + deljl[1]*deljl[1] + deljl[2]*deljl[2];
  const F_FLOAT rjl = sqrt(rsqjl);
  BOA_jl = bo_jl - thb_cut;

  cos_jil = -(delij[0]*deljl[0]+delij[1]*deljl[1]+delij[2]*deljl[2])/(rij*rjl);
  if (cos_jil > 1.0) cos_jil  = 1.0;
  if (cos_jil < -1.0) cos_jil  = -1.0;
  theta_jil = acos(cos_jil);

  // dcos_jil
  const F_FLOAT inv_distjl = 1.0 / (rij * rjl);
  const F_FLOAT cos_jil_tmp = cos_jil / ((rij*rjl)*(rij*rjl));

  for (int d = 0; d < 3; d++) {
    dcos_jil_di[d] = deljl[d] * inv_distjl - cos_jil_tmp * rsqjl * -delij[d];
    dcos_jil_dj[d] = (-deljl[d] + delij[d]) * inv_distjl - cos_jil_tmp * (rsqjl * delij[d] + rsqij * -deljl[d]);
    dcos_jil_dk[d] = -delij[d] * inv_distjl - cos_jil_tmp * rsqij * deljl[d];
  }

  sin_jil = sin(theta_jil);
  if (sin_jil >= 0 && sin_jil <= MIN_SINE)
    tan_jil_i = cos_jil / MIN_SINE;
  else if (sin_jil <= 0 && sin_jil >= -MIN_SINE)
    tan_jil_i = -cos_jil / MIN_SINE;
  else tan_jil_i = cos_jil / sin_jil;

  for (int d = 0; d < 3; d ++) dellk[d] = x(k,d) - x(l,d);
  const F_FLOAT rsqlk = dellk[0]*dellk[0] + dellk[1]*dellk[1] + dellk[2]*dellk[2];
  const F_FLOAT rlk = sqrt(rsqlk);

  // non-Kokkos ReaxFF has a separate function for computing omega, which
  //  limits the scope of the MIN_SINE statements below

  F_FLOAT sin_ijk_rnd = sin_ijk;
  F_FLOAT sin_jil_rnd = sin_jil;

  if (sin_ijk >= 0 && sin_ijk <= MIN_SINE) sin_ijk_rnd = MIN_SINE;
  else if (sin_ijk <= 0 && sin_ijk >= -MIN_SINE) sin_ijk_rnd = -MIN_SINE;
  if (sin_jil >= 0 && sin_jil <= MIN_SINE) sin_jil_rnd = MIN_SINE;
  else if (sin_jil <= 0 && sin_jil >= -MIN_SINE) sin_jil_rnd = -MIN_SINE;

  F_FLOAT unnorm_cos_omega, unnorm_sin_omega, omega;
  F_FLOAT htra, htrb, htrc, hthd, hthe, hnra, hnrc, hnhd, hnhe;
  F_FLOAT arg, poem, tel;
  F_FLOAT cross_ij_jl[3];

  // omega

  F_FLOAT dot_ij_jk = -(delij[0]*delik[0]+delij[1]*delik[1]+delij[2]*delik[2]);
  F_FLOAT dot_ij_lj = delij[0]*deljl[0]+delij[1]*deljl[1]+delij[2]*deljl[2];
  F_FLOAT dot_ik_jl = delik[0]*deljl[0]+delik[1]*deljl[1]+delik[2]*deljl[2];
  unnorm_cos_omega = dot_ij_jk * dot_ij_lj + rsqij * dot_ik_jl;

  cross_ij_jl[0] = delij[1]*deljl[2] - delij[2]*deljl[1];
  cross_ij_jl[1] = delij[2]*deljl[0] - delij[0]*deljl[2];
  cross_ij_jl[2] = delij[0]*deljl[1] - delij[1]*deljl[0];

  unnorm_sin_omega = -rij*(delik[0]*cross_ij_jl[0]+delik[1]*cross_ij_jl[1]+delik[2]*cross_ij_jl[2]);
  omega = atan2(unnorm_sin_omega, unnorm_cos_omega);

  htra = rik + cos_ijk * (rjl * cos_jil - rij);
  htrb = rij - rik * cos_ijk - rjl * cos_jil;
  htrc = rjl + cos_jil * (rik * cos_ijk - rij);
  hthd = rik * sin_ijk_rnd * (rij - rjl * cos_jil);
  hthe = rjl * sin_jil_rnd * (rij - rik * cos_ijk);
  hnra = rjl * sin_ijk_rnd * sin_jil_rnd;
  hnrc = rik * sin_ijk_rnd * sin_jil_rnd;
  hnhd = rik * rjl * cos_ijk * sin_jil_rnd;
  hnhe = rik * rjl * sin_ijk_rnd * cos_jil;

  tel = SQR(rik) + SQR(rij) + SQR(rjl) - SQR(rlk) -
        2.0 * (rik * rij * cos_ijk - rik * rjl * cos_ijk * cos_jil + rij * rjl * cos_jil);

  poem = 2.0 * rik * rjl * sin_ijk_rnd * sin_jil_rnd;
  F_FLOAT inv_poem = 1.0 / poem;

  arg = tel * inv_poem;
  if (arg >  1.0) arg =  1.0;
  if (arg < -1.0) arg = -1.0;

  cos_omega = cos(omega);
  cos2omega = cos(2. * omega);
  cos3omega = cos(3. * omega);

  // torsion energy

  p_tor1 = paramsfbp(ktype,itype,jtype,ltype).p_tor1;
  p_cot1 = paramsfbp(ktype,itype,jtype,ltype).p_cot1;
  V1 = paramsfbp(ktype,itype,jtype,ltype).V1;
  V2 = paramsfbp(ktype,itype,jtype,ltype).V2;
  V3 = paramsfbp(ktype,itype,jtype,ltype).V3;

  exp_tor1 = exp(p_tor1 * SQR(2.0 - d_BO_pi(i,j_index) - f11_DiDj));
  exp_tor2_jl = exp(-p_tor2 * BOA_jl);
  exp_cot2_jl = exp(-p_cot2 * SQR(BOA_jl - 1.5));
  fn10 = (1.0 - exp_tor2_ik) * (1.0 - exp_tor2_ij) * (1.0 - exp_tor2_jl);

  CV = 0.5 * (V1 * (1.0 + cos_omega) + V2 * exp_tor1 * (1.0 - cos2omega) + V3 * (1.0 + cos3omega));

  e_tor = fn10 * sin_ijk * sin_jil * CV;
  if (eflag) ev.ereax[6] += e_tor;

  dfn11 = (-p_tor3 * exp_tor3_DiDj + (p_tor3 * exp_tor3_DiDj - p_tor4 * exp_tor4_DiDj) *
          (2.0 + exp_tor3_DiDj) * exp_tor34_inv) * exp_tor34_inv;

  CEtors1 = sin_ijk * sin_jil * CV;

  CEtors2 = -fn10 * 2.0 * p_tor1 * V2 * exp_tor1 * (2.0 - d_BO_pi(i,j_index) - f11_DiDj) *
            (1.0 - SQR(cos_omega)) * sin_ijk * sin_jil;
  CEtors3 = CEtors2 * dfn11;

  CEtors4 = CEtors1 * p_tor2 * exp_tor2_ik * (1.0 - exp_tor2_ij) * (1.0 - exp_tor2_jl);
  CEtors5 = CEtors1 * p_tor2 * (1.0 - exp_tor2_ik) * exp_tor2_ij * (1.0 - exp_tor2_jl);
  CEtors6 = CEtors1 * p_tor2 * (1.0 - exp_tor2_ik) * (1.0 - exp_tor2_ij) * exp_tor2_jl;

  cmn = -fn10 * CV;
  CEtors7 = cmn * sin_jil * tan_ijk_i;
  CEtors8 = cmn * sin_ijk * tan_jil_i;

  CEtors9 = fn10 * sin_ijk * sin_jil *
    (0.5 * V1 - 2.0 * V2 * exp_tor1 * cos_omega + 1.5 * V3 * (cos2omega + 2.0 * SQR(cos_omega)));

  // 4-body conjugation energy

  fn12 = exp_cot2_ik * exp_cot2_ij * exp_cot2_jl;
  e_con = p_cot1 * fn12 * (1.0 + (SQR(cos_omega) - 1.0) * sin_ijk * sin_jil);
  if (eflag) ev.ereax[7] += e_con;

  Cconj = -2.0 * fn12 * p_cot1 * p_cot2 * (1.0 + (SQR(cos_omega) - 1.0) * sin_ijk * sin_jil);

  CEconj1 = Cconj * (BOA_ik - 1.5e0);
  CEconj2 = Cconj * (BOA_ij - 1.5e0);
  CEconj3 = Cconj * (BOA_jl - 1.5e0);

  CEconj4 = -p_cot1 * fn12 * (SQR(cos_omega) - 1.0) * sin_jil * tan_ijk_i;
  CEconj5 = -p_cot1 * fn12 * (SQR(cos_omega) - 1.0) * sin_ijk * tan_jil_i;
  CEconj6 = 2.0 * p_cot1 * fn12 * cos_omega * sin_ijk * sin_jil;

  // forces

  // contribution to bond order

  a_Cdbopi(i,j_index) += CEtors2;

  a_CdDelta[j] += CEtors3;
  a_CdDelta[i] += CEtors3;

  a_Cdbo(i,k_index) += CEtors4 + CEconj1;
  a_Cdbo(i,j_index) += CEtors5 + CEconj2;
  a_Cdbo(j,l_index) += CEtors6 + CEconj3;

  const F_FLOAT coeff74 = CEtors7 + CEconj4;
  const F_FLOAT coeff85 = CEtors8 + CEconj5;
  const F_FLOAT coeff96 = CEtors9 + CEconj6;

  const F_FLOAT inv_rij = 1.0 / rij;
  const F_FLOAT inv_rik = 1.0 / rik;
  const F_FLOAT inv_rjl = 1.0 / rjl;
  const F_FLOAT inv_sin_ijk_rnd = 1.0 / sin_ijk_rnd;
  const F_FLOAT inv_sin_jil_rnd = 1.0 / sin_jil_rnd;

#ifdef LMP_KOKKOS_GPU
  #pragma unroll
#endif
  for (int d = 0; d < 3; d++) {
    // dcos_omega_di
    F_FLOAT dcos_omega_dk = ((htra-arg*hnra) * inv_rik) * delik[d] - dellk[d];
    dcos_omega_dk += (hthd-arg*hnhd) * inv_sin_ijk_rnd * -dcos_ijk_dk[d];
    dcos_omega_dk *= 2.0 * inv_poem;

    // dcos_omega_dj
    F_FLOAT dcos_omega_di = -((htra-arg*hnra) * inv_rik) * delik[d] - htrb * inv_rij * delij[d];
    dcos_omega_di += -(hthd-arg*hnhd) * inv_sin_ijk_rnd * dcos_ijk_di[d];
    dcos_omega_di += -(hthe-arg*hnhe) * inv_sin_jil_rnd * dcos_jil_di[d];
    dcos_omega_di *= 2.0 * inv_poem;

    // dcos_omega_dk
    F_FLOAT dcos_omega_dj = -((htrc-arg*hnrc) * inv_rjl) * deljl[d] + htrb * inv_rij * delij[d];
    dcos_omega_dj += -(hthd-arg*hnhd) * inv_sin_ijk_rnd * dcos_ijk_dj[d];
    dcos_omega_dj += -(hthe-arg*hnhe) * inv_sin_jil_rnd * dcos_jil_dj[d];
    dcos_omega_dj *= 2.0 * inv_poem;

    // dcos_omega_dl
    F_FLOAT dcos_omega_dl = ((htrc-arg*hnrc) * inv_rjl) * deljl[d] + dellk[d];
    dcos_omega_dl += (hthe-arg*hnhe) * inv_sin_jil_rnd * -dcos_jil_dk[d];
    dcos_omega_dl *= 2.0 * inv_poem;

    // dcos_theta_ijk
    fi_tmp[d] = (coeff74) * dcos_ijk_di[d];
    fj_tmp[d] = (coeff74) * dcos_ijk_dj[d];
    fk_tmp[d] = (coeff74) * dcos_ijk_dk[d];

    // dcos_theta_jil
    fi_tmp[d] += (coeff85) * dcos_jil_di[d];
    fj_tmp[d] += (coeff85) * dcos_jil_dj[d];
    F_FLOAT fl_tmp =  (coeff85) * dcos_jil_dk[d];

    // dcos_omega
    fi_tmp[d] += (coeff96) * dcos_omega_di;
    fj_tmp[d] += (coeff96) * dcos_omega_dj;
    fk_tmp[d] += (coeff96) * dcos_omega_dk;
    fl_tmp += (coeff96) * dcos_omega_dl;

    // total forces
    a_f(i,d) -= fi_tmp[d];
    a_f(j,d) -= fj_tmp[d];
    a_f(k,d) -= fk_tmp[d];
    a_f(l,d) -= fl_tmp;
  }

  // per-atom energy/virial tally

  if (EVFLAG) {
    eng_tmp = e_tor + e_con;
    //if (eflag_atom) this->template ev_tally<NEIGHFLAG>(ev,i,j,eng_tmp,0.0,0.0,0.0,0.0);
    if (eflag_atom) this->template e_tally<NEIGHFLAG>(ev,i,j,eng_tmp);
    if (vflag_either) {
        for (int d = 0; d < 3; d ++) delil[d] = x(l,d) - x(i,d);
        for (int d = 0; d < 3; d ++) delkl[d] = x(l,d) - x(k,d);
        this->template v_tally4<NEIGHFLAG>(ev,k,i,j,l,fk_tmp,fi_tmp,fj_tmp,delkl,delil,deljl);
    }
  }
}

template<class DeviceType>
template<int NEIGHFLAG, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::operator()(TagPairReaxComputeTorsionPreprocessed<NEIGHFLAG,EVFLAG>, const int &tpack) const {
  EV_FLOAT_REAX ev;
  this->template operator()<NEIGHFLAG,EVFLAG>(TagPairReaxComputeTorsionPreprocessed<NEIGHFLAG,EVFLAG>(), tpack, ev);

}





#endif
