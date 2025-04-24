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

#ifndef LMP_PAIR_REAXFF_BONDORDER_KOKKOS_HPP
#define LMP_PAIR_REAXFF_BONDORDER_KOKKOS_HPP

#include "pair_reaxff_kokkos.h"

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::operator()(TagPairReaxBondOrder1, const int &ii) const {

  const int i = d_ilist[ii];
  const int itype = type(i);

  const F_FLOAT val_i = paramssing(itype).valency;
  d_Deltap[i] = d_total_bo[i] - val_i;
  d_Deltap_boc[i] = d_total_bo[i] - paramssing(itype).valency_val;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int EREAXFF_FLAG>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::operator()(TagPairReaxBondOrder2<EREAXFF_FLAG>, const int &ii) const {

  F_FLOAT exp_p1i, exp_p2i, exp_p1j, exp_p2j, f1, f2, f3, u1_ij, u1_ji, Cf1A_ij, Cf1B_ij, Cf1_ij, Cf1_ji;
  F_FLOAT f4, f5, exp_f4, exp_f5, f4f5, Cf45_ij, Cf45_ji;
  F_FLOAT A0_ij, A1_ij, A2_ij, A3_ij, A2_ji, A3_ji;

  const int i = d_ilist[ii];
  const int itype = type(i);
  const int jnum = d_bo_num[i];

  const F_FLOAT val_i = paramssing(itype).valency;

  d_total_bo[i] = 0.0;
  F_FLOAT total_bo = 0.0;

  for (int j_index = 0; j_index < jnum; j_index++) {
    int j = d_bo_list(i, j_index);
    j &= NEIGHMASK;
    const int jtype = type(j);

    // calculate corrected BO and total bond order
    const F_FLOAT val_j = paramssing(jtype).valency;
    const F_FLOAT ovc = paramstwbp(itype,jtype).ovc;
    const F_FLOAT v13cor = paramstwbp(itype,jtype).v13cor;
    const F_FLOAT p_boc3 = paramstwbp(itype,jtype).p_boc3;
    const F_FLOAT p_boc4 = paramstwbp(itype,jtype).p_boc4;
    const F_FLOAT p_boc5 = paramstwbp(itype,jtype).p_boc5;

    if (ovc < 0.001 && v13cor < 0.001) {
      d_C1dbo(i,j_index) = 1.0;
      d_C2dbo(i,j_index) = 0.0;
      d_C3dbo(i,j_index) = 0.0;
      d_C1dbopi(i,j_index) = 1.0;
      d_C2dbopi(i,j_index) = 0.0;
      d_C3dbopi(i,j_index) = 0.0;
      d_C4dbopi(i,j_index) = 0.0;
      d_C1dbopi2(i,j_index) = 1.0;
      d_C2dbopi2(i,j_index) = 0.0;
      d_C3dbopi2(i,j_index) = 0.0;
      d_C4dbopi2(i,j_index) = 0.0;
    } else {
      if (ovc >= 0.001) {
        exp_p1i = exp(-p_boc1 * d_Deltap[i]);
        exp_p2i = exp(-p_boc2 * d_Deltap[i]);
        exp_p1j = exp(-p_boc1 * d_Deltap[j]);
        exp_p2j = exp(-p_boc2 * d_Deltap[j]);

        f2 = exp_p1i + exp_p1j;
        f3 = -1.0/p_boc2*log(0.5*(exp_p2i+exp_p2j));
        f1 = 0.5 * ((val_i + f2)/(val_i + f2 + f3) + (val_j + f2)/(val_j + f2 + f3));
        u1_ij = val_i + f2 + f3;
        u1_ji = val_j + f2 + f3;
        Cf1A_ij = 0.5 * f3 * (1.0/(u1_ij*u1_ij)+1.0/(u1_ji*u1_ji));
        Cf1B_ij = -0.5 * ((u1_ij - f3)/(u1_ij*u1_ij)+(u1_ji - f3)/(u1_ji*u1_ji));
        Cf1_ij = 0.5 * (-p_boc1 * exp_p1i / u1_ij - ((val_i+f2) / (u1_ij*u1_ij)) *
                       (-p_boc1 * exp_p1i + exp_p2i / (exp_p2i + exp_p2j)) +
                        -p_boc1 * exp_p1i / u1_ji - ((val_j+f2) / (u1_ji*u1_ji)) *
                       (-p_boc1 * exp_p1i + exp_p2i / (exp_p2i + exp_p2j)));
        Cf1_ji = -Cf1A_ij * p_boc1 * exp_p1j + Cf1B_ij * exp_p2j / (exp_p2i + exp_p2j);
      } else {
        f1 = 1.0;
        Cf1_ij = Cf1_ji = 0.0;
      }

      if (v13cor >= 0.001) {
        exp_f4 =exp(-(p_boc4*(d_BO(i,j_index)*d_BO(i,j_index))-d_Deltap_boc[i])*p_boc3+p_boc5);
        exp_f5 =exp(-(p_boc4*(d_BO(i,j_index)*d_BO(i,j_index))-d_Deltap_boc[j])*p_boc3+p_boc5);
        f4 = 1. / (1. + exp_f4);
        f5 = 1. / (1. + exp_f5);
        f4f5 = f4 * f5;

        Cf45_ij = -f4 * exp_f4;
        Cf45_ji = -f5 * exp_f5;
      } else {
        f4 = f5 = f4f5 = 1.0;
        Cf45_ij = Cf45_ji = 0.0;
      }

      A0_ij = f1 * f4f5;
      A1_ij = -2 * p_boc3 * p_boc4 * d_BO(i,j_index) * (Cf45_ij + Cf45_ji);
      A2_ij = Cf1_ij / f1 + p_boc3 * Cf45_ij;
      A2_ji = Cf1_ji / f1 + p_boc3 * Cf45_ji;
      A3_ij = A2_ij + Cf1_ij / f1;
      A3_ji = A2_ji + Cf1_ji / f1;

      d_BO(i,j_index) = d_BO(i,j_index) * A0_ij;
      d_BO_pi(i,j_index) = d_BO_pi(i,j_index) * A0_ij * f1;
      d_BO_pi2(i,j_index) = d_BO_pi2(i,j_index) * A0_ij * f1;
      d_BO_s(i,j_index) = d_BO(i,j_index)-(d_BO_pi(i,j_index)+d_BO_pi2(i,j_index));

      d_C1dbo(i,j_index) = A0_ij + d_BO(i,j_index) * A1_ij;
      d_C2dbo(i,j_index) = d_BO(i,j_index) * A2_ij;
      d_C3dbo(i,j_index) = d_BO(i,j_index) * A2_ji;

      d_C1dbopi(i,j_index) = f1*f1*f4*f5;
      d_C2dbopi(i,j_index) = d_BO_pi(i,j_index) * A1_ij;
      d_C3dbopi(i,j_index) = d_BO_pi(i,j_index) * A3_ij;
      d_C4dbopi(i,j_index) = d_BO_pi(i,j_index) * A3_ji;

      d_C1dbopi2(i,j_index) = f1*f1*f4*f5;
      d_C2dbopi2(i,j_index) = d_BO_pi2(i,j_index) * A1_ij;
      d_C3dbopi2(i,j_index) = d_BO_pi2(i,j_index) * A3_ij;
      d_C4dbopi2(i,j_index) = d_BO_pi2(i,j_index) * A3_ji;
    }

    if (d_BO(i,j_index) < 1e-10) d_BO(i,j_index) = 0.0;
    if (d_BO_s(i,j_index) < 1e-10) d_BO_s(i,j_index) = 0.0;
    if (d_BO_pi(i,j_index) < 1e-10) d_BO_pi(i,j_index) = 0.0;
    if (d_BO_pi2(i,j_index) < 1e-10) d_BO_pi2(i,j_index) = 0.0;

    total_bo += d_BO(i,j_index);

    d_Cdbo(i,j_index) = 0.0;
    d_Cdbopi(i,j_index) = 0.0;
    d_Cdbopi2(i,j_index) = 0.0;
    d_CdDelta[j] = 0.0;
  }
  d_CdDelta[i] = 0.0;
  d_total_bo[i] += total_bo;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int EREAXFF_FLAG>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::operator()(TagPairReaxBondOrder3<EREAXFF_FLAG>, const int &ii) const {
  // bot part of BO()

  const int i = d_ilist[ii];
  const int itype = type(i);
  F_FLOAT nlp_temp;

  if(EREAXFF_FLAG) {

    d_Delta[i] = d_total_bo[i] - paramssing(itype).valency;

  } else
    d_Delta[i] = d_total_bo[i] - paramssing(itype).valency;

  const F_FLOAT Delta_e = d_total_bo[i] - paramssing(itype).valency_e;
  d_Delta_boc[i] = d_total_bo[i] - paramssing(itype).valency_boc;

  const F_FLOAT vlpex = Delta_e - 2.0 * (int)(Delta_e/2.0);
  const F_FLOAT explp1 = exp(-gp[15] * SQR(2.0 + vlpex));
  const F_FLOAT nlp = explp1 - (int)(Delta_e / 2.0);
  d_Delta_lp[i] = paramssing(itype).nlp_opt - nlp;
  const F_FLOAT Clp = 2.0 * gp[15] * explp1 * (2.0 + vlpex);
  d_dDelta_lp[i] = Clp;

  if (paramssing(itype).mass > 21.0) {
    nlp_temp = 0.5 * (paramssing(itype).valency_e - paramssing(itype).valency);
    d_Delta_lp_temp[i] = paramssing(itype).nlp_opt - nlp_temp;
  } else {
    nlp_temp = nlp;
    d_Delta_lp_temp[i] = paramssing(itype).nlp_opt - nlp_temp;
  }

  d_sum_ovun(i,1) = 0.0;
  d_sum_ovun(i,2) = 0.0;
}

#endif
