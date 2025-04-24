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

#ifndef LMP_PAIR_REAXFF_OTHER_KOKKOS_HPP
#define LMP_PAIR_REAXFF_OTHER_KOKKOS_HPP

#include "pair_reaxff_kokkos.h"

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::operator()(TagPairReaxComputePolar<NEIGHFLAG>, const int &ii, EV_FLOAT_REAX& ev) const {

  const int i = d_ilist[ii];
  const int itype = type(i);
  const F_FLOAT qi = q(i);
  const F_FLOAT chi = paramssing(itype).chi;
  const F_FLOAT eta = paramssing(itype).eta;

  F_FLOAT epol = KCALpMOL_to_EV*(chi*qi+(eta/2.0)*qi*qi);

  /* energy due to coupling with kinetic energy potential */
  if (acks2_flag)
    epol += KCALpMOL_to_EV*qi*d_s[NN + i];

  if (eflag_global) ev.ecoul += epol;
  //if (eflag_atom) this->template ev_tally<NEIGHFLAG>(ev,i,i,epol,0.0,0.0,0.0,0.0);
  if (eflag_atom) this->template e_tally_single<NEIGHFLAG>(ev,i,epol);
}

template<class DeviceType>
template<int NEIGHFLAG>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::operator()(TagPairReaxComputePolar<NEIGHFLAG>, const int &ii) const {
  EV_FLOAT_REAX ev;
  this->template operator()<NEIGHFLAG>(TagPairReaxComputePolar<NEIGHFLAG>(), ii, ev);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::operator()(TagPairReaxComputeTabulatedLJCoulomb<NEIGHFLAG,EVFLAG>, const int &ii, EV_FLOAT_REAX& ev) const {

  // The f array is duplicated for OpenMP, atomic for GPU, and neither for Serial

  auto v_f = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_f),decltype(ndup_f)>::get(dup_f,ndup_f);
  auto a_f = v_f.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  const int i = d_ilist[ii];
  const X_FLOAT xtmp = x(i,0);
  const X_FLOAT ytmp = x(i,1);
  const X_FLOAT ztmp = x(i,2);
  const F_FLOAT qi = q(i);
  const int itype = type(i);
  const tagint itag = tag(i);
  const int jnum = d_numneigh[i];

  F_FLOAT fxtmp, fytmp, fztmp;
  fxtmp = fytmp = fztmp = 0.0;

  for (int jj = 0; jj < jnum; jj++) {
    int j = d_neighbors(i,jj);
    j &= NEIGHMASK;
    const int jtype = type(j);
    const tagint jtag = tag(j);
    const F_FLOAT qj = q(j);

    // skip half of the interactions
    if (j >= nlocal) {
      if (itag > jtag) {
        if ((itag+jtag) % 2 == 0) continue;
      } else if (itag < jtag) {
        if ((itag+jtag) % 2 == 1) continue;
      } else {
        if (x(j,2) < ztmp) continue;
        if (x(j,2) == ztmp && x(j,1)  < ytmp) continue;
        if (x(j,2) == ztmp && x(j,1) == ytmp && x(j,0) < xtmp) continue;
      }
    }

    const X_FLOAT delx = x(j,0) - xtmp;
    const X_FLOAT dely = x(j,1) - ytmp;
    const X_FLOAT delz = x(j,2) - ztmp;
    const F_FLOAT rsq = delx*delx + dely*dely + delz*delz;

    if (rsq > cut_nbsq) continue;
    const F_FLOAT rij = sqrt(rsq);

    const int tmin  = MIN(itype, jtype);
    const int tmax  = MAX(itype, jtype);
    const LR_lookup_table_kk<DeviceType>& t = k_LR.template view<DeviceType>()(tmin,tmax);


    /* Cubic Spline Interpolation */
    int r = (int)(rij * t.inv_dx);
    if (r == 0)  ++r;
    const F_FLOAT base = (double)(r+1) * t.dx;
    const F_FLOAT dif = rij - base;

    const cubic_spline_coef vdW = t.d_vdW[r];
    const cubic_spline_coef ele = t.d_ele[r];
    const cubic_spline_coef CEvd = t.d_CEvd[r];
    const cubic_spline_coef CEclmb = t.d_CEclmb[r];

    const F_FLOAT evdwl = ((vdW.d*dif + vdW.c)*dif + vdW.b)*dif +
      vdW.a;

    F_FLOAT ecoul = (((ele.d*dif + ele.c)*dif + ele.b)*dif +
      ele.a)*qi*qj;

    const F_FLOAT fvdwl = ((CEvd.d*dif + CEvd.c)*dif + CEvd.b)*dif +
      CEvd.a;

    F_FLOAT fcoul = (((CEclmb.d*dif+CEclmb.c)*dif+CEclmb.b)*dif +
      CEclmb.a)*qi*qj;

    /* contribution to energy and gradients (atoms and cell)
     * due to geometry-dependent terms in the ACKS2
     * kinetic energy */
    if (acks2_flag) {

      /* kinetic energy terms */
      double xcut = 0.5 * (paramssing(itype).bcut_acks2
                          + paramssing(jtype).bcut_acks2);

      if (rij <= xcut) {
        const F_FLOAT d = rij / xcut;
        const F_FLOAT bond_softness = gp[34] * pow( d, 3.0 )
                                    * pow( 1.0 - d, 6.0 );

        if (bond_softness > 0.0) {
          /* Coulombic energy contribution */
          const F_FLOAT effpot_diff = d_s[NN + i]
                                    - d_s[NN + j];
          const F_FLOAT e_ele = -0.5 * KCALpMOL_to_EV * bond_softness
                                     * SQR( effpot_diff );

          ecoul += e_ele;

          /* forces contribution */
          F_FLOAT d_bond_softness;
          d_bond_softness = gp[34]
                          * 3.0 / xcut * pow( d, 2.0 )
                          * pow( 1.0 - d, 5.0 ) * (1.0 - 3.0 * d);
          d_bond_softness = -0.5 * d_bond_softness
                          * SQR( effpot_diff );
          d_bond_softness = KCALpMOL_to_EV * d_bond_softness
                          / rij;

          fcoul += d_bond_softness;
        }
      }
    }

    const F_FLOAT ftotal = fvdwl + fcoul;
    fxtmp += delx*ftotal;
    fytmp += dely*ftotal;
    fztmp += delz*ftotal;
    a_f(j,0) -= delx*ftotal;
    a_f(j,1) -= dely*ftotal;
    a_f(j,2) -= delz*ftotal;

    if (EVFLAG) {
      if (eflag_global) ev.evdwl += evdwl;
      if (eflag_global) ev.ecoul += ecoul;

      if (vflag_either || eflag_atom) this->template ev_tally<NEIGHFLAG>(ev,i,j,evdwl+ecoul,-ftotal,delx,dely,delz);
    }
  }

  a_f(i,0) += fxtmp;
  a_f(i,1) += fytmp;
  a_f(i,2) += fztmp;
}

template<class DeviceType>
template<int NEIGHFLAG, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::operator()(TagPairReaxComputeTabulatedLJCoulomb<NEIGHFLAG,EVFLAG>, const int &ii) const {
  EV_FLOAT_REAX ev;
  this->template operator()<NEIGHFLAG,EVFLAG>(TagPairReaxComputeTabulatedLJCoulomb<NEIGHFLAG,EVFLAG>(), ii, ev);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::operator()(TagPairReaxComputeMulti1, const int &ii) const {

  const int i = d_ilist[ii];
  const int itype = type(i);
  const F_FLOAT imass = paramssing(itype).mass;
  F_FLOAT dfvl;

  if (imass > 21.0) dfvl = 0.0;
  else dfvl = 1.0;

  const int jnum = d_bo_num[i];

  F_FLOAT sum_ovun1 = 0.0;
  F_FLOAT sum_ovun2 = 0.0;

  for (int j_index = 0; j_index < jnum; j_index++) {
    int j = d_bo_list(i, j_index);
    j &= NEIGHMASK;
    const int jtype = type(j);

    sum_ovun1 += paramstwbp(itype,jtype).p_ovun1 * paramstwbp(itype,jtype).De_s * d_BO(i,j_index);
    sum_ovun2 += (d_Delta[j] - dfvl * d_Delta_lp_temp[j]) * (d_BO_pi(i, j_index) + d_BO_pi2(i,j_index));
  }
  d_sum_ovun(i,1) += sum_ovun1;
  d_sum_ovun(i,2) += sum_ovun2;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG, int EFLAG>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::operator()(TagPairReaxComputeMulti2<NEIGHFLAG,EFLAG>, const int &ii, EV_FLOAT_REAX& ev) const {

  auto v_CdDelta = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_CdDelta),decltype(ndup_CdDelta)>::get(dup_CdDelta,ndup_CdDelta);
  auto a_CdDelta = v_CdDelta.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  const int i = d_ilist[ii];
  const int itype = type(i);
  const F_FLOAT imass = paramssing(itype).mass;
  const F_FLOAT val_i = paramssing(itype).valency;

  F_FLOAT dfvl;
  if (imass > 21.0) dfvl = 0.0;
  else dfvl = 1.0;

  F_FLOAT e_lp, e_ov, e_un;
  F_FLOAT CEover1, CEover2, CEover3, CEover4;
  F_FLOAT CEunder1, CEunder2, CEunder3, CEunder4;
  const F_FLOAT p_lp3 = gp[5];
  const F_FLOAT p_ovun2 = paramssing(itype).p_ovun2;
  const F_FLOAT p_ovun3 = gp[32];
  const F_FLOAT p_ovun4 = gp[31];
  const F_FLOAT p_ovun5 = paramssing(itype).p_ovun5;
  const F_FLOAT p_ovun6 = gp[6];
  const F_FLOAT p_ovun7 = gp[8];
  const F_FLOAT p_ovun8 = gp[9];

  // lone pair
  const F_FLOAT p_lp2 = paramssing(itype).p_lp2;
  const F_FLOAT expvd2 = exp(-75 * d_Delta_lp[i]);
  const F_FLOAT inv_expvd2 = 1.0 / (1.0+expvd2);

  int numbonds = d_bo_num[i];

  e_lp = 0.0;
  if (numbonds > 0 || enobondsflag)
    e_lp = p_lp2 * d_Delta_lp[i] * inv_expvd2;
  const F_FLOAT dElp = p_lp2 * inv_expvd2 + 75.0 * p_lp2 * d_Delta_lp[i] * expvd2 * inv_expvd2*inv_expvd2;
  const F_FLOAT CElp = dElp * d_dDelta_lp[i];

  if (numbonds > 0 || enobondsflag)
    a_CdDelta[i] += CElp;

  if (EFLAG && eflag_global) ev.ereax[0] += e_lp;
  //if (vflag_either || eflag_atom) this->template ev_tally<NEIGHFLAG>(ev,i,i,e_lp,0.0,0.0,0.0,0.0);
  //if (eflag_atom) this->template e_tally<NEIGHFLAG>(ev,i,i,e_lp);

  // over coordination
  const F_FLOAT exp_ovun1 = p_ovun3 * exp(p_ovun4 * d_sum_ovun(i,2));
  const F_FLOAT inv_exp_ovun1 = 1.0 / (1 + exp_ovun1);
  const F_FLOAT Delta_lpcorr  = d_Delta[i] - (dfvl * d_Delta_lp_temp[i]) * inv_exp_ovun1;

  const F_FLOAT exp_ovun2 = exp(p_ovun2 * Delta_lpcorr);
  const F_FLOAT inv_exp_ovun2 = 1.0 / (1.0 + exp_ovun2);
  const F_FLOAT DlpVi = 1.0 / (Delta_lpcorr + val_i + 1e-8);

  CEover1 = Delta_lpcorr * DlpVi * inv_exp_ovun2;
  e_ov = d_sum_ovun(i,1) * CEover1;

  if (EFLAG && eflag_global) ev.ereax[1] += e_ov;
  //if (eflag_atom) this->template ev_tally<NEIGHFLAG>(ev,i,i,e_ov,0.0,0.0,0.0,0.0);
  //if (eflag_atom) this->template e_tally<NEIGHFLAG>(ev,i,i,e_ov);

  CEover2 = d_sum_ovun(i,1) * DlpVi * inv_exp_ovun2 *
    (1.0 - Delta_lpcorr * (DlpVi + p_ovun2 * exp_ovun2 * inv_exp_ovun2));
  CEover3 = CEover2 * (1.0 - dfvl * d_dDelta_lp[i] * inv_exp_ovun1);
  CEover4 = CEover2 * (dfvl * d_Delta_lp_temp[i]) * p_ovun4 * exp_ovun1 * SQR(inv_exp_ovun1);

  // under coordination

  const F_FLOAT exp_ovun2n = 1.0 / exp_ovun2;
  const F_FLOAT exp_ovun6 = exp(p_ovun6 * Delta_lpcorr);
  const F_FLOAT exp_ovun8 = p_ovun7 * exp(p_ovun8 * d_sum_ovun(i,2));
  const F_FLOAT inv_exp_ovun2n = 1.0 / (1.0 + exp_ovun2n);
  const F_FLOAT inv_exp_ovun8 = 1.0 / (1.0 + exp_ovun8);

  e_un = 0;
  if (numbonds > 0 || enobondsflag)
    e_un = -p_ovun5 * (1.0 - exp_ovun6) * inv_exp_ovun2n * inv_exp_ovun8;

  if (EFLAG && eflag_global) ev.ereax[2] += e_un;
  //if (eflag_atom) this->template ev_tally<NEIGHFLAG>(ev,i,i,e_un,0.0,0.0,0.0,0.0);
  //if (eflag_atom) this->template e_tally<NEIGHFLAG>(ev,i,i,e_un);

  CEunder1 = inv_exp_ovun2n *
    (p_ovun5 * p_ovun6 * exp_ovun6 * inv_exp_ovun8 + p_ovun2 * e_un * exp_ovun2n);
  CEunder2 = -e_un * p_ovun8 * exp_ovun8 * inv_exp_ovun8;
  CEunder3 = CEunder1 * (1.0 - dfvl * d_dDelta_lp[i] * inv_exp_ovun1);
  CEunder4 = CEunder1 * (dfvl * d_Delta_lp_temp[i]) *
      p_ovun4 * exp_ovun1 * inv_exp_ovun1 * inv_exp_ovun1 + CEunder2;

  const F_FLOAT eng_tmp = e_lp + e_ov + e_un;
  if (eflag_atom) this->template e_tally_single<NEIGHFLAG>(ev,i,eng_tmp);

  // multibody forces

  a_CdDelta[i] += CEover3;
  if (numbonds > 0 || enobondsflag)
    a_CdDelta[i] += CEunder3;

  const int jnum = d_bo_num[i];

  F_FLOAT CdDelta_i = 0.0;
  for (int j_index = 0; j_index < jnum; j_index++) {
    int j = d_bo_list(i, j_index);
    j &= NEIGHMASK;
    const int jtype = type(j);
    const F_FLOAT jmass = paramssing(jtype).mass;
    const F_FLOAT De_s = paramstwbp(itype,jtype).De_s;

    // multibody lone pair: correction for C2
    if (p_lp3 > 0.001 && imass == 12.0 && jmass == 12.0) {
      const F_FLOAT Di = d_Delta[i];
      const F_FLOAT vov3 = d_BO(i,j_index) - Di - 0.040*pow(Di,4.0);
      if (vov3 > 3.0) {
        const F_FLOAT e_lph = p_lp3 * (vov3-3.0)*(vov3-3.0);
        const F_FLOAT deahu2dbo = 2.0 * p_lp3 * (vov3 - 3.0);
        const F_FLOAT deahu2dsbo = 2.0 * p_lp3 * (vov3 - 3.0) * (-1.0 - 0.16 * pow(Di,3.0));
        d_Cdbo(i,j_index) += deahu2dbo;
        CdDelta_i += deahu2dsbo;

        if (EFLAG) {
          if (eflag_global) ev.ereax[0] += e_lph;
          if (eflag_atom) this->template e_tally<NEIGHFLAG>(ev,i,j,e_lph);
        }
      }
    }

    // over/under coordination forces merged together
    const F_FLOAT p_ovun1 = paramstwbp(itype,jtype).p_ovun1;
    a_CdDelta[j] += (CEover4 + CEunder4) * (1.0 - dfvl * d_dDelta_lp[j]) * (d_BO_pi(i,j_index) + d_BO_pi2(i,j_index));
    d_Cdbo(i,j_index) += CEover1 * p_ovun1 * De_s;
    d_Cdbopi(i,j_index) += (CEover4 + CEunder4) * (d_Delta[j] - dfvl*d_Delta_lp_temp[j]);
    d_Cdbopi2(i,j_index) += (CEover4 + CEunder4) * (d_Delta[j] - dfvl*d_Delta_lp_temp[j]);
  }
  a_CdDelta[i] += CdDelta_i;

}

template<class DeviceType>
template<int NEIGHFLAG, int EFLAG>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::operator()(TagPairReaxComputeMulti2<NEIGHFLAG,EFLAG>, const int &ii) const {
  EV_FLOAT_REAX ev;
  this->template operator()<NEIGHFLAG,EFLAG>(TagPairReaxComputeMulti2<NEIGHFLAG,EFLAG>(), ii, ev);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
double PairReaxFFKokkos<DeviceType>::memory_usage()
{
  double bytes = 0.0;

  if (cut_hbsq > 0.0) {
    bytes += (double)nmax*3*sizeof(int);
    bytes += (double)maxhb*nmax*sizeof(int);
  }
  bytes += (double)nmax*2*sizeof(int);
  bytes += (double)maxbo*nmax*sizeof(int);

  bytes += (double)nmax*17*sizeof(F_FLOAT);
  bytes += (double)maxbo*nmax*34*sizeof(F_FLOAT);

  // FixReaxFFSpecies
  if (fixspecies_flag) {
    bytes += (double)MAXSPECBOND*nmax*sizeof(tagint);
    bytes += (double)MAXSPECBOND*nmax*sizeof(F_FLOAT);
  }

  // FixReaxFFBonds
  bytes += (double)maxbo*nmax*sizeof(tagint);
  bytes += (double)maxbo*nmax*sizeof(F_FLOAT);
  bytes += (double)nmax*sizeof(int);

  return bytes;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::ev_tally(EV_FLOAT_REAX &ev, const int &i, const int &j,
      const F_FLOAT &epair, const F_FLOAT &fpair, const F_FLOAT &delx,
                const F_FLOAT &dely, const F_FLOAT &delz) const
{
  // The eatom and vatom arrays are duplicated for OpenMP, atomic for GPU, and neither for Serial

  auto v_eatom = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_eatom),decltype(ndup_eatom)>::get(dup_eatom,ndup_eatom);
  auto a_eatom = v_eatom.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  auto v_vatom = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_vatom),decltype(ndup_vatom)>::get(dup_vatom,ndup_vatom);
  auto a_vatom = v_vatom.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  if (eflag_atom) {
    const E_FLOAT epairhalf = 0.5 * epair;
    a_eatom[i] += epairhalf;
    a_eatom[j] += epairhalf;
  }

  if (vflag_either) {
    const E_FLOAT v0 = delx*delx*fpair;
    const E_FLOAT v1 = dely*dely*fpair;
    const E_FLOAT v2 = delz*delz*fpair;
    const E_FLOAT v3 = delx*dely*fpair;
    const E_FLOAT v4 = delx*delz*fpair;
    const E_FLOAT v5 = dely*delz*fpair;

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
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::e_tally(EV_FLOAT_REAX & /*ev*/, const int &i, const int &j,
      const F_FLOAT &epair) const
{
  // The eatom array is duplicated for OpenMP, atomic for GPU, and neither for Serial

  auto v_eatom = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_eatom),decltype(ndup_eatom)>::get(dup_eatom,ndup_eatom);
  auto a_eatom = v_eatom.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  const E_FLOAT epairhalf = 0.5 * epair;
  a_eatom[i] += epairhalf;
  a_eatom[j] += epairhalf;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::e_tally_single(EV_FLOAT_REAX & /*ev*/, const int &i,
      const F_FLOAT &epair) const
{
  // The eatom array is duplicated for OpenMP, atomic for GPU, and neither for Serial

  auto v_eatom = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_eatom),decltype(ndup_eatom)>::get(dup_eatom,ndup_eatom);
  auto a_eatom = v_eatom.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  a_eatom[i] += epair;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::v_tally(EV_FLOAT_REAX &ev, const int &i,
  F_FLOAT *fi, F_FLOAT *drij) const
{
  F_FLOAT v[6];

  v[0] = 0.5*drij[0]*fi[0];
  v[1] = 0.5*drij[1]*fi[1];
  v[2] = 0.5*drij[2]*fi[2];
  v[3] = 0.5*drij[0]*fi[1];
  v[4] = 0.5*drij[0]*fi[2];
  v[5] = 0.5*drij[1]*fi[2];

  if (vflag_global) {
    ev.v[0] += v[0];
    ev.v[1] += v[1];
    ev.v[2] += v[2];
    ev.v[3] += v[3];
    ev.v[4] += v[4];
    ev.v[5] += v[5];
  }

  if (vflag_atom) {
    auto v_vatom = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_vatom),decltype(ndup_vatom)>::get(dup_vatom,ndup_vatom);
    auto a_vatom = v_vatom.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

    a_vatom(i,0) += v[0]; a_vatom(i,1) += v[1]; a_vatom(i,2) += v[2];
    a_vatom(i,3) += v[3]; a_vatom(i,4) += v[4]; a_vatom(i,5) += v[5];
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::v_tally3(EV_FLOAT_REAX &ev, const int &i, const int &j, const int &k,
  F_FLOAT *fj, F_FLOAT *fk, F_FLOAT *drij, F_FLOAT *drik) const
{
  // The eatom and vatom arrays are duplicated for OpenMP, atomic for GPU, and neither for Serial
  auto v_vatom = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_vatom),decltype(ndup_vatom)>::get(dup_vatom,ndup_vatom);
  auto a_vatom = v_vatom.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  F_FLOAT v[6];

  v[0] = drij[0]*fj[0] + drik[0]*fk[0];
  v[1] = drij[1]*fj[1] + drik[1]*fk[1];
  v[2] = drij[2]*fj[2] + drik[2]*fk[2];
  v[3] = drij[0]*fj[1] + drik[0]*fk[1];
  v[4] = drij[0]*fj[2] + drik[0]*fk[2];
  v[5] = drij[1]*fj[2] + drik[1]*fk[2];

  if (vflag_global) {
    ev.v[0] += v[0];
    ev.v[1] += v[1];
    ev.v[2] += v[2];
    ev.v[3] += v[3];
    ev.v[4] += v[4];
    ev.v[5] += v[5];
  }

  if (vflag_atom) {
    a_vatom(i,0) += THIRD * v[0]; a_vatom(i,1) += THIRD * v[1]; a_vatom(i,2) += THIRD * v[2];
    a_vatom(i,3) += THIRD * v[3]; a_vatom(i,4) += THIRD * v[4]; a_vatom(i,5) += THIRD * v[5];
    a_vatom(j,0) += THIRD * v[0]; a_vatom(j,1) += THIRD * v[1]; a_vatom(j,2) += THIRD * v[2];
    a_vatom(j,3) += THIRD * v[3]; a_vatom(j,4) += THIRD * v[4]; a_vatom(j,5) += THIRD * v[5];
    a_vatom(k,0) += THIRD * v[0]; a_vatom(k,1) += THIRD * v[1]; a_vatom(k,2) += THIRD * v[2];
    a_vatom(k,3) += THIRD * v[3]; a_vatom(k,4) += THIRD * v[4]; a_vatom(k,5) += THIRD * v[5];
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::v_tally4(EV_FLOAT_REAX &ev, const int &i, const int &j, const int &k,
  const int &l, F_FLOAT *fi, F_FLOAT *fj, F_FLOAT *fk, F_FLOAT *dril, F_FLOAT *drjl, F_FLOAT *drkl) const
{
  // The vatom array is duplicated for OpenMP, atomic for GPU, and neither for Serial

  F_FLOAT v[6];

  v[0] = dril[0]*fi[0] + drjl[0]*fj[0] + drkl[0]*fk[0];
  v[1] = dril[1]*fi[1] + drjl[1]*fj[1] + drkl[1]*fk[1];
  v[2] = dril[2]*fi[2] + drjl[2]*fj[2] + drkl[2]*fk[2];
  v[3] = dril[0]*fi[1] + drjl[0]*fj[1] + drkl[0]*fk[1];
  v[4] = dril[0]*fi[2] + drjl[0]*fj[2] + drkl[0]*fk[2];
  v[5] = dril[1]*fi[2] + drjl[1]*fj[2] + drkl[1]*fk[2];

  if (vflag_global) {
    ev.v[0] += v[0];
    ev.v[1] += v[1];
    ev.v[2] += v[2];
    ev.v[3] += v[3];
    ev.v[4] += v[4];
    ev.v[5] += v[5];
  }

  if (vflag_atom) {
    auto v_vatom = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_vatom),decltype(ndup_vatom)>::get(dup_vatom,ndup_vatom);
    auto a_vatom = v_vatom.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

    a_vatom(i,0) += 0.25 * v[0]; a_vatom(i,1) += 0.25 * v[1]; a_vatom(i,2) += 0.25 * v[2];
    a_vatom(i,3) += 0.25 * v[3]; a_vatom(i,4) += 0.25 * v[4]; a_vatom(i,5) += 0.25 * v[5];
    a_vatom(j,0) += 0.25 * v[0]; a_vatom(j,1) += 0.25 * v[1]; a_vatom(j,2) += 0.25 * v[2];
    a_vatom(j,3) += 0.25 * v[3]; a_vatom(j,4) += 0.25 * v[4]; a_vatom(j,5) += 0.25 * v[5];
    a_vatom(k,0) += 0.25 * v[0]; a_vatom(k,1) += 0.25 * v[1]; a_vatom(k,2) += 0.25 * v[2];
    a_vatom(k,3) += 0.25 * v[3]; a_vatom(k,4) += 0.25 * v[4]; a_vatom(k,5) += 0.25 * v[5];
    a_vatom(l,0) += 0.25 * v[0]; a_vatom(l,1) += 0.25 * v[1]; a_vatom(l,2) += 0.25 * v[2];
    a_vatom(l,3) += 0.25 * v[3]; a_vatom(l,4) += 0.25 * v[4]; a_vatom(l,5) += 0.25 * v[5];
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::v_tally3_atom(EV_FLOAT_REAX &ev, const int &i, const int & /*j*/,
                                                const int & /*k*/, F_FLOAT *fj, F_FLOAT *fk,
                                                F_FLOAT *drji, F_FLOAT *drjk) const
{
  F_FLOAT v[6];

  v[0] = THIRD * (drji[0]*fj[0] + drjk[0]*fk[0]);
  v[1] = THIRD * (drji[1]*fj[1] + drjk[1]*fk[1]);
  v[2] = THIRD * (drji[2]*fj[2] + drjk[2]*fk[2]);
  v[3] = THIRD * (drji[0]*fj[1] + drjk[0]*fk[1]);
  v[4] = THIRD * (drji[0]*fj[2] + drjk[0]*fk[2]);
  v[5] = THIRD * (drji[1]*fj[2] + drjk[1]*fk[2]);

  if (vflag_global) {
    ev.v[0] += v[0];
    ev.v[1] += v[1];
    ev.v[2] += v[2];
    ev.v[3] += v[3];
    ev.v[4] += v[4];
    ev.v[5] += v[5];
  }

  if (vflag_atom) {
    d_vatom(i,0) += v[0]; d_vatom(i,1) += v[1]; d_vatom(i,2) += v[2];
    d_vatom(i,3) += v[3]; d_vatom(i,4) += v[4]; d_vatom(i,5) += v[5];
  }
}

/* ----------------------------------------------------------------------
   setup for energy, virial computation
   see integrate::ev_set() for values of eflag and vflag
   see pair::ev_setup() for values of eflag_* and vflag_*
   VIRIAL_CENTROID bitflag is not yet supported by ReaxFF
------------------------------------------------------------------------- */

template<class DeviceType>
void PairReaxFFKokkos<DeviceType>::ev_setup(int eflag, int vflag, int)
{
  int i;

  evflag = 1;

  eflag_either = eflag;
  eflag_global = eflag & ENERGY_GLOBAL;
  eflag_atom = eflag & ENERGY_ATOM;

  vflag_either = vflag;
  vflag_global = vflag & (VIRIAL_PAIR | VIRIAL_FDOTR);
  vflag_atom = vflag & VIRIAL_ATOM;

  // reallocate per-atom arrays if necessary

  if (eflag_atom && atom->nmax > maxeatom) {
    maxeatom = atom->nmax;
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->create_kokkos(k_eatom,eatom,maxeatom,"pair:eatom");
    d_eatom = k_eatom.view<DeviceType>();
  }
  if (vflag_atom && atom->nmax > maxvatom) {
    maxvatom = atom->nmax;
    memoryKK->destroy_kokkos(k_vatom,vatom);
    memoryKK->create_kokkos(k_vatom,vatom,maxvatom,"pair:vatom");
    d_vatom = k_vatom.view<DeviceType>();
  }

  // zero accumulators

  if (eflag_global) eng_vdwl = eng_coul = 0.0;
  if (vflag_global) for (i = 0; i < 6; i++) virial[i] = 0.0;
  if (eflag_atom) Kokkos::deep_copy(d_eatom,0.0);
  if (vflag_atom) Kokkos::deep_copy(d_vatom,0.0);

  // if vflag_global = VIRIAL_FDOTR and pair::compute() calls virial_fdotr_compute()
  // compute global virial via (F dot r) instead of via pairwise summation
  // unset other flags as appropriate

  if (vflag_global == VIRIAL_FDOTR && no_virial_fdotr_compute == 0) {
    vflag_fdotr = 1;
    vflag_global = 0;
    if (vflag_atom == 0) vflag_either = 0;
    if (vflag_either == 0 && eflag_either == 0) evflag = 0;
  } else vflag_fdotr = 0;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairReaxFFKokkos<DeviceType>::FindBond(int &numbonds, int groupbit)
{
  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxFindBondZero>(0,nmax),*this);

  bo_cut_bond = api->control->bg_cut;

  atomKK->sync(execution_space,TAG_MASK|MASK_MASK);
  tag = atomKK->k_tag.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();

  const int inum = list->inum;
  NeighListKokkos<DeviceType>* k_list = static_cast<NeighListKokkos<DeviceType>*>(list);
  d_ilist = k_list->d_ilist;

  numbonds = 0;
  PairReaxKokkosFindBondFunctor<DeviceType> find_bond_functor(this, groupbit);
  Kokkos::parallel_reduce(inum,find_bond_functor,numbonds);
  copymode = 0;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::operator()(TagPairReaxFindBondZero, const int &i) const {
  d_numneigh_bonds[i] = 0;
  for (int j = 0; j < maxbo; j++) {
    d_neighid(i,j) = 0;
    d_abo(i,j) = 0.0;
  }
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::calculate_find_bond_item(int ii, int &numbonds, int groupbit) const
{
  const int i = d_ilist[ii];
  int nj = 0;

  if (mask[i] & groupbit) {
    const int jnum = d_bo_num[i];
    for (int j_index = 0; j_index < jnum; j_index++) {
      int j = d_bo_list(i, j_index);
      j &= NEIGHMASK;
      if (mask[j] & groupbit) {
        const tagint jtag = tag[j];
        double bo_tmp = d_BO(i, j_index);

        if (bo_tmp > bo_cut_bond) {
          d_neighid(i,nj) = jtag;
          d_abo(i,nj) = bo_tmp;
          nj++;
        }
      }
    }
  }
  d_numneigh_bonds[i] = nj;
  if (nj > numbonds) numbonds = nj;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairReaxFFKokkos<DeviceType>::PackBondBuffer(DAT::tdual_ffloat_1d k_buf, int &nbuf_local)
{
  d_buf = k_buf.view<DeviceType>();
  k_params_sing.template sync<DeviceType>();
  atomKK->sync(execution_space,TAG_MASK|TYPE_MASK|Q_MASK|MOLECULE_MASK);

  tag = atomKK->k_tag.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  q = atomKK->k_q.view<DeviceType>();
  if (atom->molecule)
    molecule = atomKK->k_molecule.view<DeviceType>();

  copymode = 1;
  nlocal = atomKK->nlocal;
  PairReaxKokkosPackBondBufferFunctor<DeviceType> pack_bond_buffer_functor(this);
  Kokkos::parallel_scan(nlocal,pack_bond_buffer_functor);
  copymode = 0;

  k_buf.modify<DeviceType>();
  k_nbuf_local.modify<DeviceType>();

  k_buf.sync<LMPHostType>();
  k_nbuf_local.sync<LMPHostType>();
  nbuf_local = k_nbuf_local.h_view();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairReaxFFKokkos<DeviceType>::PackReducedBondBuffer(DAT::tdual_ffloat_1d k_buf, int &nbuf_local, bool store_bonds)
{
  d_buf = k_buf.view<DeviceType>();
  k_params_sing.template sync<DeviceType>();

  copymode = 1;
  nlocal = atomKK->nlocal;
  if (store_bonds) {
    PairReaxKokkosPackReducedBondBufferFunctor<DeviceType, true> pack_bond_buffer_functor(this);
    Kokkos::parallel_scan(nlocal,pack_bond_buffer_functor);
  } else {
    PairReaxKokkosPackReducedBondBufferFunctor<DeviceType, false> pack_bond_buffer_functor(this);
    Kokkos::parallel_scan(nlocal,pack_bond_buffer_functor);
  }

  copymode = 0;

  k_buf.modify<DeviceType>();
  k_nbuf_local.modify<DeviceType>();

  k_buf.sync<LMPHostType>();
  k_nbuf_local.sync<LMPHostType>();
  nbuf_local = k_nbuf_local.h_view();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::pack_bond_buffer_item(int i, int &j, const bool &final) const
{
  if (i == 0)
    j += 2;

  if (final) {
    d_buf[j-1] = tag[i];
    d_buf[j+0] = type[i];
    d_buf[j+1] = d_total_bo[i];
    d_buf[j+2] = paramssing(type[i]).nlp_opt - d_Delta_lp[i];
    d_buf[j+3] = q[i];
    d_buf[j+4] = d_numneigh_bonds[i];
  }
  const int numbonds = d_numneigh_bonds[i];

  if (final) {
    for (int k = 5; k < 5+numbonds; k++) {
      d_buf[j+k] = d_neighid(i,k-5);
    }
  }
  j += (5+numbonds);

  if (final) {
    if (!molecule.data()) d_buf[j] = 0.0;
    else d_buf[j] = molecule[i];
  }
  j++;

  if (final) {
    for (int k = 0; k < numbonds; k++) {
      d_buf[j+k] = d_abo(i,k);
    }
  }
  j += (1+numbonds);

  if (final && i == nlocal-1)
    k_nbuf_local.view<DeviceType>()() = j - 1;
}

template<class DeviceType>
template<bool STORE_BONDS>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::pack_reduced_bond_buffer_item(int i, int &j, const bool &final) const
{
  const int numbonds = d_numneigh_bonds[i];
  if (final) {
    d_buf[j] = d_total_bo[i];
    d_buf[j+1] = paramssing(type[i]).nlp_opt - d_Delta_lp[i];
    d_buf[j+2] = numbonds;
  }

  j += 3;

  if constexpr(STORE_BONDS) {
    if (final) {
      for (int k = 0; k < numbonds; ++k) {
        d_buf[j+k] = d_neighid(i,k);
      }
    }

    j += numbonds;

    if (final) {
      for (int k = 0; k < numbonds; k++) {
        d_buf[j+k] = d_abo(i,k);
      }
    }

    j += numbonds;
  }

  if (final && i == nlocal-1)
    k_nbuf_local.view<DeviceType>()() = j - 1;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairReaxFFKokkos<DeviceType>::FindBondSpecies()
{

  if (nmax > (int)k_tmpid.extent(0)) {
    memoryKK->destroy_kokkos(k_tmpid,tmpid);
    memoryKK->destroy_kokkos(k_tmpbo,tmpbo);
    memoryKK->create_kokkos(k_tmpid,tmpid,nmax,MAXSPECBOND,"pair:tmpid");
    memoryKK->create_kokkos(k_tmpbo,tmpbo,nmax,MAXSPECBOND,"pair:tmpbo");
  }

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxFindBondSpeciesZero>(0,nmax),*this);

  nlocal = atomKK->nlocal;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxFindBondSpecies>(0,nlocal),*this);
  copymode = 0;

  // NOTE: Could improve performance if a Kokkos version of ComputeSpecAtom is added

  k_tmpbo.modify<DeviceType>();
  k_tmpid.modify<DeviceType>();
  k_error_flag.modify<DeviceType>();

  k_tmpbo.sync<LMPHostType>();
  k_tmpid.sync<LMPHostType>();
  k_error_flag.sync<LMPHostType>();

  if (k_error_flag.h_view())
    error->all(FLERR,"Increase MAXSPECBOND in reaxff_defs.h");
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::operator()(TagPairReaxFindBondSpeciesZero, const int &i) const {
  for (int j = 0; j < MAXSPECBOND; j++) {
    k_tmpbo.view<DeviceType>()(i,j) = 0.0;
    k_tmpid.view<DeviceType>()(i,j) = 0;
  }
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::operator()(TagPairReaxFindBondSpecies, const int &i) const {
  int nj = 0;

  const int jnum = d_bo_num[i];
  for (int j_index = 0; j_index < jnum; j_index++) {
    int j = d_bo_list(i, j_index);
    j &= NEIGHMASK;
    if (j < i) continue;

    double bo_tmp = d_BO(i, j_index);

    if (bo_tmp >= 0.10) { // Why is this a hardcoded value?
      k_tmpid.view<DeviceType>()(i,nj) = j;
      k_tmpbo.view<DeviceType>()(i,nj) = bo_tmp;
      nj++;
      if (nj > MAXSPECBOND) k_error_flag.view<DeviceType>()() = 1;
    }
  }
}


#endif
