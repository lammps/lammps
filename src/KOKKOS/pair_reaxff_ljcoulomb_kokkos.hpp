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

#ifndef LMP_PAIR_REAXFF_LJCOULOMB_KOKKOS_HPP
#define LMP_PAIR_REAXFF_LJCOULOMB_KOKKOS_HPP

#include "pair_reaxff_kokkos.h"

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairReaxFFKokkos<DeviceType>::operator()(TagPairReaxComputeLJCoulomb<NEIGHFLAG,EVFLAG>, const int &ii, EV_FLOAT_REAX& ev) const {

  // The f array is duplicated for OpenMP, atomic for GPU, and neither for Serial

  auto v_f = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_f),decltype(ndup_f)>::get(dup_f,ndup_f);
  auto a_f = v_f.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  F_FLOAT powr_vdw, powgi_vdw, fn13, dfn13, exp1, exp2, etmp;
  F_FLOAT evdwl, fvdwl;
  evdwl = fvdwl = 0.0;

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

    // LJ energy/force
    F_FLOAT Tap = d_tap[7] * rij + d_tap[6];
    Tap = Tap * rij + d_tap[5];
    Tap = Tap * rij + d_tap[4];
    Tap = Tap * rij + d_tap[3];
    Tap = Tap * rij + d_tap[2];
    Tap = Tap * rij + d_tap[1];
    Tap = Tap * rij + d_tap[0];

    F_FLOAT dTap = 7*d_tap[7] * rij + 6*d_tap[6];
    dTap = dTap * rij + 5*d_tap[5];
    dTap = dTap * rij + 4*d_tap[4];
    dTap = dTap * rij + 3*d_tap[3];
    dTap = dTap * rij + 2*d_tap[2];
    dTap += d_tap[1]/rij;

    const F_FLOAT gamma_w = paramstwbp(itype,jtype).gamma_w;
    const F_FLOAT alpha = paramstwbp(itype,jtype).alpha;
    const F_FLOAT r_vdw = paramstwbp(itype,jtype).r_vdw;
    const F_FLOAT epsilon = paramstwbp(itype,jtype).epsilon;

    // shielding
    if (vdwflag == 1 || vdwflag == 3) {
      F_FLOAT tmp_var;
      tmp_var = pow(rij,gp[28]-2.0);
      powr_vdw = tmp_var*rij*rij;
      powgi_vdw = pow(1.0/gamma_w,gp[28]);
      dfn13 = pow(powr_vdw+powgi_vdw,1.0/gp[28]-1.0);
      fn13  = dfn13*(powr_vdw+powgi_vdw);
      dfn13 = dfn13*tmp_var;

      exp2 = exp(0.5*alpha*(1.0-fn13/r_vdw));
      exp1 = exp2*exp2;
      etmp = epsilon*(exp1-2.0*exp2);
      evdwl = Tap*etmp;
      fvdwl = dTap*etmp-Tap*epsilon*(alpha/r_vdw)*(exp1-exp2)*dfn13;
    } else {
      exp2 = exp(0.5*alpha*(1.0-rij/r_vdw));
      exp1 = exp2*exp2;
      etmp = epsilon*(exp1-2.0*exp2);
      evdwl = Tap*etmp;
      fvdwl = dTap*etmp-Tap*epsilon*(alpha/r_vdw)*(exp1-exp2)*rij;
    }
    // inner wall
    if (vdwflag == 2 || vdwflag == 3) {
      const F_FLOAT ecore = paramstwbp(itype,jtype).ecore;
      const F_FLOAT acore = paramstwbp(itype,jtype).acore;
      const F_FLOAT rcore = paramstwbp(itype,jtype).rcore;
      const F_FLOAT e_core = ecore*exp(acore*(1.0-(rij/rcore)));
      const F_FLOAT de_core = -(acore/rcore)*e_core;
      evdwl += Tap*e_core;
      fvdwl += dTap*e_core+Tap*de_core/rij;

      if (lgflag) {
        const F_FLOAT lgre = paramstwbp(itype,jtype).lgre;
        const F_FLOAT lgcij = paramstwbp(itype,jtype).lgcij;
        const F_FLOAT rij5 = rsq*rsq*rij;
        const F_FLOAT rij6 = rij5*rij;
        const F_FLOAT re6 = lgre*lgre*lgre*lgre*lgre*lgre;
        const F_FLOAT elg = -lgcij/(rij6+re6);
        const F_FLOAT delg = -6.0*elg*rij5/(rij6+re6);
        evdwl += Tap*elg;
        fvdwl += dTap*elg+Tap*delg/rij;
      }
    }

    // Coulomb energy/force
    const F_FLOAT shld = paramstwbp(itype,jtype).gamma;
    const F_FLOAT denom1 = rij * rij * rij + shld;
    const F_FLOAT denom3 = cbrt(denom1);
    F_FLOAT ecoul = C_ele * qi*qj*Tap/denom3;
    F_FLOAT fcoul = C_ele * qi*qj*(dTap-Tap*rij/denom1)/denom3;

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
    a_f(j,0) -= delx*ftotal;
    fytmp += dely*ftotal;
    a_f(j,1) -= dely*ftotal;
    fztmp += delz*ftotal;
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
void PairReaxFFKokkos<DeviceType>::operator()(TagPairReaxComputeLJCoulomb<NEIGHFLAG,EVFLAG>, const int &ii) const {
  EV_FLOAT_REAX ev;
  this->template operator()<NEIGHFLAG,EVFLAG>(TagPairReaxComputeLJCoulomb<NEIGHFLAG,EVFLAG>(), ii, ev);
}

#endif
