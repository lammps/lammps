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

#ifndef LMP_PAIR_REAXFF_HYDROGEN_KOKKOS_HPP
#define LMP_PAIR_REAXFF_HYDROGEN_KOKKOS_HPP

#include "pair_reaxff_kokkos.h"

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::operator()(TagPairReaxComputeHydrogen<NEIGHFLAG,EVFLAG>, const int &ii, EV_FLOAT_REAX& ev) const {

  auto v_f = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_f),decltype(ndup_f)>::get(dup_f,ndup_f);
  auto a_f = v_f.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  int hblist[MAX_BONDS];
  F_FLOAT theta, cos_theta, sin_xhz4, cos_xhz1, sin_theta2;
  F_FLOAT e_hb, exp_hb2, exp_hb3, CEhb1, CEhb2, CEhb3;
  F_FLOAT dcos_theta_di[3], dcos_theta_dj[3], dcos_theta_dk[3];

  // tally variables
  F_FLOAT fi_tmp[3], fj_tmp[3], fk_tmp[3], delij[3], delji[3], delik[3], delki[3];
  for (int d = 0; d < 3; d++) fi_tmp[d] = fj_tmp[d] = fk_tmp[d] = 0.0;

  const int i = d_ilist[ii];
  const int itype = type(i);
  if (paramssing(itype).p_hbond != 1) return;

  const X_FLOAT xtmp = x(i,0);
  const X_FLOAT ytmp = x(i,1);
  const X_FLOAT ztmp = x(i,2);

  const int jnum = d_bo_num[i];
  const int knum = d_hb_num[i];

  int top = 0;
  for (int j_index = 0; j_index < jnum; j_index++) {
    int j = d_bo_list(i, j_index);
    j &= NEIGHMASK;
    const int jtype = type(j);
    const F_FLOAT bo_ij = d_BO(i,j_index);

    if (paramssing(jtype).p_hbond == 2 && bo_ij >= HB_THRESHOLD) {
      hblist[top] = j_index;
      top ++;
    }
  }

  F_FLOAT fitmp[3];
  for (int d = 0; d < 3; d++) fitmp[d] = 0.0;

  for (int k_index = 0; k_index < knum; k_index++) {
    int k = d_hb_list(i, k_index);
    k &= NEIGHMASK;
    const tagint ktag = tag(k);
    const int ktype = type(k);

    delik[0] = x(k,0) - xtmp;
    delik[1] = x(k,1) - ytmp;
    delik[2] = x(k,2) - ztmp;
    const F_FLOAT rsqik = delik[0]*delik[0] + delik[1]*delik[1] + delik[2]*delik[2];
    const F_FLOAT rik = sqrt(rsqik);

    for (int itr = 0; itr < top; itr++) {
      const int j_index = hblist[itr];
      int j = d_bo_list(i, j_index);
      j &= NEIGHMASK;
      const tagint jtag = tag(j);
      if (jtag == ktag) continue;

      const int jtype = type(j);
      const F_FLOAT bo_ij = d_BO(i,j_index);

      delij[0] = x(j,0) - xtmp;
      delij[1] = x(j,1) - ytmp;
      delij[2] = x(j,2) - ztmp;
      const F_FLOAT rsqij = delij[0]*delij[0] + delij[1]*delij[1] + delij[2]*delij[2];
      const F_FLOAT rij = sqrt(rsqij);

      // theta and derivatives
      cos_theta = (delij[0]*delik[0]+delij[1]*delik[1]+delij[2]*delik[2])/(rij*rik);
      if (cos_theta > 1.0) cos_theta  = 1.0;
      if (cos_theta < -1.0) cos_theta  = -1.0;
      theta = acos(cos_theta);

      const F_FLOAT inv_dists = 1.0 / (rij * rik);
      const F_FLOAT Cdot_inv3 = cos_theta * inv_dists * inv_dists;

      for (int d = 0; d < 3; d++) {
        dcos_theta_di[d] = -(delik[d] + delij[d]) * inv_dists + Cdot_inv3 * (rsqik * delij[d] + rsqij * delik[d]);
        dcos_theta_dj[d] = delik[d] * inv_dists - Cdot_inv3 * rsqik * delij[d];
        dcos_theta_dk[d] = delij[d] * inv_dists - Cdot_inv3 * rsqij * delik[d];
      }

      // hydrogen bond energy
      const F_FLOAT p_hb1 = paramshbp(jtype,itype,ktype).p_hb1;
      const F_FLOAT p_hb2 = paramshbp(jtype,itype,ktype).p_hb2;
      const F_FLOAT p_hb3 = paramshbp(jtype,itype,ktype).p_hb3;
      const F_FLOAT r0_hb = paramshbp(jtype,itype,ktype).r0_hb;

      sin_theta2 = sin(theta/2.0);
      sin_xhz4 = SQR(sin_theta2);
      sin_xhz4 *= sin_xhz4;
      cos_xhz1 = (1.0 - cos_theta);
      exp_hb2 = exp(-p_hb2 * bo_ij);
      exp_hb3 = exp(-p_hb3 * (r0_hb/rik + rik/r0_hb - 2.0));

      e_hb = p_hb1 * (1.0 - exp_hb2) * exp_hb3 * sin_xhz4;
      if (EVFLAG && eflag_global) ev.ereax[8] += e_hb;

      // hydrogen bond forces
      CEhb1 = p_hb1 * p_hb2 * exp_hb2 * exp_hb3 * sin_xhz4;
      CEhb2 = -p_hb1/2.0 * (1.0 - exp_hb2) * exp_hb3 * cos_xhz1;
      CEhb3 = -p_hb3 * (-r0_hb/SQR(rik) + 1.0/r0_hb) * e_hb;

      d_Cdbo(i,j_index) += CEhb1; // dbo term

      // dcos terms
      for (int d = 0; d < 3; d++) fi_tmp[d] = CEhb2 * dcos_theta_di[d];
      for (int d = 0; d < 3; d++) fj_tmp[d] = CEhb2 * dcos_theta_dj[d];
      for (int d = 0; d < 3; d++) fk_tmp[d] = CEhb2 * dcos_theta_dk[d];

      // dr terms
      for (int d = 0; d < 3; d++) fi_tmp[d] -= CEhb3/rik * delik[d];
      for (int d = 0; d < 3; d++) fk_tmp[d] += CEhb3/rik * delik[d];

      for (int d = 0; d < 3; d++) fitmp[d] -= fi_tmp[d];
      for (int d = 0; d < 3; d++) a_f(j,d) -= fj_tmp[d];
      for (int d = 0; d < 3; d++) a_f(k,d) -= fk_tmp[d];

      for (int d = 0; d < 3; d++) delki[d] = -1.0 * delik[d];
      for (int d = 0; d < 3; d++) delji[d] = -1.0 * delij[d];

      if (EVFLAG) {
        if (eflag_atom) this->template e_tally<NEIGHFLAG>(ev,i,j,e_hb);
        if (vflag_either) this->template v_tally3<NEIGHFLAG>(ev,i,j,k,fj_tmp,fk_tmp,delji,delki);
      }
    }
  }
  for (int d = 0; d < 3; d++) a_f(i,d) += fitmp[d];
}

template<class DeviceType>
template<int NEIGHFLAG, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::operator()(TagPairReaxComputeHydrogen<NEIGHFLAG,EVFLAG>, const int &ii) const {
  EV_FLOAT_REAX ev;
  this->template operator()<NEIGHFLAG,EVFLAG>(TagPairReaxComputeHydrogen<NEIGHFLAG,EVFLAG>(), ii, ev);

}





#endif
