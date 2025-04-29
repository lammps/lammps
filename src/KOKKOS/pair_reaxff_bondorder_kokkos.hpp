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

// clang-format off
/* ----------------------------------------------------------------------
   Unified bond-order kernel with optimisations 1-7 applied
------------------------------------------------------------------------- */

/*

template<int EREAXFF_FLAG>
struct TagPairReaxBondOrderUnified {};

template<class DeviceType>
template<int EREAXFF_FLAG>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::operator()(
        TagPairReaxBondOrderUnified<EREAXFF_FLAG>,
        const TeamPolicy<DeviceType>::member_type &team) const
{
  const int ii = team.league_rank();
  const int i  = d_ilist[ii];
  const int itype = type(i);
  const int jnum  = d_bo_num[i];

  // scratch block
  const int lvl = 0;
  F_FLOAT *s_bo   = (F_FLOAT*) team.team_shmem().get_shmem(jnum*4*sizeof(F_FLOAT));
  int     *s_j    = (int*)    team.team_shmem().get_shmem(jnum*sizeof(int));

  // preload neighbours into scratch – AoSoA: [bo,bo_pi,bo_pi2,bo_s | j]
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team,jnum), [&](int jj){
      const int jidx = d_bo_list(i,jj);
      s_j[jj]      = jidx & NEIGHMASK;
      s_bo[4*jj]   = d_BO(i,jj);
      s_bo[4*jj+1] = d_BO_pi(i,jj);
      s_bo[4*jj+2] = d_BO_pi2(i,jj);
      s_bo[4*jj+3] = 0.0;      // placeholder for BO_s
  });
  team.team_barrier();

  if(team.team_rank()==0){
    d_Deltap[i]     = d_total_bo[i] - c_paramssing(itype*PSSING_STRIDE+VALENCY);
    d_Deltap_boc[i] = d_total_bo[i] - c_paramssing(itype*PSSING_STRIDE+VALENCY_VAL);
    d_total_bo[i]   = 0.0;
  }
  team.team_barrier();

  // hierarchical j-loop
  F_FLOAT total_bo_team = 0.0;

  Kokkos::parallel_reduce(Kokkos::TeamVectorRange(team,jnum),
  [&](const int jj, F_FLOAT &bo_acc){

    F_FLOAT &bo   = s_bo[4*jj];
    if(bo < 1e-12) return;                      // early exit (#6)

    const F_FLOAT &bo_pi  = s_bo[4*jj+1];
    const F_FLOAT &bo_pi2 = s_bo[4*jj+2];

    const int j      = s_j[jj];
    const int jtype  = type(j);

    const F_FLOAT val_i = c_paramssing(itype*PSSING_STRIDE+VALENCY);
    const F_FLOAT val_j = c_paramssing(jtype*PSSING_STRIDE+VALENCY);

    const int twbp_off = (itype*n_species + jtype)*PWB_STRIDE;
    const F_FLOAT ovc    = c_paramstwbp(twbp_off+OVC);
    const F_FLOAT v13cor = c_paramstwbp(twbp_off+V13COR);
    const F_FLOAT p_boc1 = c_paramstwbp(twbp_off+PBOC1);
    const F_FLOAT p_boc2 = c_paramstwbp(twbp_off+PBOC2);
    const F_FLOAT p_boc3 = c_paramstwbp(twbp_off+PBOC3);
    const F_FLOAT p_boc4 = c_paramstwbp(twbp_off+PBOC4);
    const F_FLOAT p_boc5 = c_paramstwbp(twbp_off+PBOC5);

    const F_FLOAT Dpi = d_Deltap[i], Dpj = d_Deltap[j];
    const F_FLOAT Dboci = d_Deltap_boc[i], Dbocj = d_Deltap_boc[j];

    F_FLOAT f1=1,f2=0,f3=0,Cf1_i=0,Cf1_j=0;
    if(ovc>=0.001){
      const F_FLOAT exp_p1i=exp(-p_boc1*Dpi);
      const F_FLOAT exp_p1j=exp(-p_boc1*Dpj);
      const F_FLOAT exp_p2i=exp(-p_boc2*Dpi);
      const F_FLOAT exp_p2j=exp(-p_boc2*Dpj);
      f2 = exp_p1i+exp_p1j;
      f3 = -1.0/p_boc2*log(0.5*(exp_p2i+exp_p2j));
      f1 = 0.5*((val_i+f2)/(val_i+f2+f3) + (val_j+f2)/(val_j+f2+f3));
      const F_FLOAT u1_i = val_i+f2+f3;
      const F_FLOAT u1_j = val_j+f2+f3;
      Cf1_i = -0.5*p_boc1*exp_p1i/u1_i -0.5*(val_i+f2)/(u1_i*u1_i)*
              (-p_boc1*exp_p1i+exp_p2i/(exp_p2i+exp_p2j));
      Cf1_j = -0.5*p_boc1*exp_p1j/u1_j -0.5*(val_j+f2)/(u1_j*u1_j)*
              (-p_boc1*exp_p1j+exp_p2j/(exp_p2i+exp_p2j));
    }

    F_FLOAT f4=1,f5=1,Cf4=0,Cf5=0;
    if(v13cor>=0.001){
      const F_FLOAT bo2 = bo*bo;
      const F_FLOAT exp_f4 = exp(-(p_boc4*bo2-Dboci)*p_boc3+p_boc5);
      const F_FLOAT exp_f5 = exp(-(p_boc4*bo2-Dbocj)*p_boc3+p_boc5);
      f4  = 1.0/(1.0+exp_f4);
      f5  = 1.0/(1.0+exp_f5);
      Cf4 = -f4*exp_f4;
      Cf5 = -f5*exp_f5;
    }

    const F_FLOAT A0 = f1*f4*f5;
    const F_FLOAT A1 = -2*p_boc3*p_boc4*bo*(Cf4+Cf5);
    const F_FLOAT A2 = Cf1_i/f1 + p_boc3*Cf4;
    const F_FLOAT A3 = Cf1_j/f1 + p_boc3*Cf5;

    // fused updates (#5)
    bo           = std::fma(bo,          A0, 0.0);
    const F_FLOAT bo_pi_up  = std::fma(bo_pi ,A0*f1,0.0);
    const F_FLOAT bo_pi2_up = std::fma(bo_pi2,A0*f1,0.0);
    const F_FLOAT bo_s_up   = bo - bo_pi_up - bo_pi2_up;

    // store back
    d_BO(i,jj)       = bo;
    d_BO_pi(i,jj)    = bo_pi_up;
    d_BO_pi2(i,jj)   = bo_pi2_up;
    d_BO_s(i,jj)     = bo_s_up;

    d_C1dbo(i,jj)    = A0 + bo*A1;
    d_C2dbo(i,jj)    = bo*A2;
    d_C3dbo(i,jj)    = bo*A3;

    d_C1dbopi(i,jj)  = f1*f1*f4*f5;
    d_C2dbopi(i,jj)  = bo_pi_up*A1;
    d_C3dbopi(i,jj)  = bo_pi_up*(A2+A3);

    d_C1dbopi2(i,jj) = f1*f1*f4*f5;
    d_C2dbopi2(i,jj) = bo_pi2_up*A1;
    d_C3dbopi2(i,jj) = bo_pi2_up*(A2+A3);

    bo_acc += bo;
  }, total_bo_team);

  // warp-aggregated atomic add (#7)
  KokkosKernels::Impl::warp_atomic_add(&d_total_bo[i], total_bo_team);

  team.team_barrier();

  if(team.team_rank()==0){
    d_Delta[i]     = d_total_bo[i] - c_paramssing(itype*PSSING_STRIDE+VALENCY);
    d_Delta_boc[i] = d_total_bo[i] - c_paramssing(itype*PSSING_STRIDE+VALENCY_BOC);
  }
}

template<class DeviceType, int EREAXFF_FLAG>
using PairReaxFFBondOrderUnifiedPolicy =
  Kokkos::TeamPolicy<DeviceType,
    TagPairReaxBondOrderUnified<EREAXFF_FLAG>,
    32, 8>; // launch-bounds (#1, #9)

*/

#endif
