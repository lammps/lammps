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

#include "pair_lj_cut_tip4p_cut_kokkos.h"

#include "atom_kokkos.h"

#include <cmath>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairLJCutTIP4PCutKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
  const int inum = this->prepare(eflag_in, vflag_in);

  // copy per-type-pair LJ coefficients to device (small, refreshed each call)

  const int ntp1 = this->atom->ntypes + 1;
  if ((int)d_lj1.extent(0) != ntp1) {
    d_lj1 = typename AT::t_kkfloat_2d("lj/cut/tip4p/cut/kk:lj1",ntp1,ntp1);
    d_lj2 = typename AT::t_kkfloat_2d("lj/cut/tip4p/cut/kk:lj2",ntp1,ntp1);
    d_lj3 = typename AT::t_kkfloat_2d("lj/cut/tip4p/cut/kk:lj3",ntp1,ntp1);
    d_lj4 = typename AT::t_kkfloat_2d("lj/cut/tip4p/cut/kk:lj4",ntp1,ntp1);
    d_offset = typename AT::t_kkfloat_2d("lj/cut/tip4p/cut/kk:offset",ntp1,ntp1);
    d_cut_ljsq = typename AT::t_kkfloat_2d("lj/cut/tip4p/cut/kk:cut_ljsq",ntp1,ntp1);
  }
  {
    auto h_lj1 = Kokkos::create_mirror_view(d_lj1);
    auto h_lj2 = Kokkos::create_mirror_view(d_lj2);
    auto h_lj3 = Kokkos::create_mirror_view(d_lj3);
    auto h_lj4 = Kokkos::create_mirror_view(d_lj4);
    auto h_offset = Kokkos::create_mirror_view(d_offset);
    auto h_cut_ljsq = Kokkos::create_mirror_view(d_cut_ljsq);
    for (int i = 1; i < ntp1; i++)
      for (int j = 1; j < ntp1; j++) {
        h_lj1(i,j) = this->lj1[i][j]; h_lj2(i,j) = this->lj2[i][j];
        h_lj3(i,j) = this->lj3[i][j]; h_lj4(i,j) = this->lj4[i][j];
        h_offset(i,j) = this->offset[i][j]; h_cut_ljsq(i,j) = this->cut_ljsq[i][j];
      }
    Kokkos::deep_copy(d_lj1,h_lj1); Kokkos::deep_copy(d_lj2,h_lj2);
    Kokkos::deep_copy(d_lj3,h_lj3); Kokkos::deep_copy(d_lj4,h_lj4);
    Kokkos::deep_copy(d_offset,h_offset); Kokkos::deep_copy(d_cut_ljsq,h_cut_ljsq);
  }

  this->copymode = 1;

  // M (virtual charge) site of every O atom (local + ghost)
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagPairTIP4PNewsite>(0,this->nall), *this);

  EV_FLOAT ev;
  if (this->eflag || this->vflag)
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType,TagPairLJCutTIP4PCutCompute<1>>(0,inum), *this, ev);
  else
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType,TagPairLJCutTIP4PCutCompute<0>>(0,inum), *this, ev);

  this->copymode = 0;

  if (this->eflag_global) this->eng_vdwl += ev.evdwl;
  this->finalize(ev);
}

/* ----------------------------------------------------------------------
   standard pairwise (LJ) energy/virial tally for a half neighbor list
------------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairLJCutTIP4PCutKokkos<DeviceType>::ev_tally(EV_FLOAT &ev, const int &i, const int &j,
    const KK_FLOAT &evdwl, const KK_FLOAT &fpair,
    const KK_FLOAT &delx, const KK_FLOAT &dely, const KK_FLOAT &delz) const
{
  if (this->eflag_global) ev.evdwl += evdwl;
  if (this->eflag_atom) {
    Kokkos::atomic_add(&d_eatom[i], (KK_ACC_FLOAT)(0.5*evdwl));
    Kokkos::atomic_add(&d_eatom[j], (KK_ACC_FLOAT)(0.5*evdwl));
  }
  if (this->vflag_global || this->vflag_atom) {
    const KK_FLOAT v0 = delx*delx*fpair;
    const KK_FLOAT v1 = dely*dely*fpair;
    const KK_FLOAT v2 = delz*delz*fpair;
    const KK_FLOAT v3 = delx*dely*fpair;
    const KK_FLOAT v4 = delx*delz*fpair;
    const KK_FLOAT v5 = dely*delz*fpair;
    if (this->vflag_global) {
      ev.v[0] += v0; ev.v[1] += v1; ev.v[2] += v2;
      ev.v[3] += v3; ev.v[4] += v4; ev.v[5] += v5;
    }
    if (this->vflag_atom) {
      Kokkos::atomic_add(&d_vatom(i,0), (KK_ACC_FLOAT)(0.5*v0));
      Kokkos::atomic_add(&d_vatom(i,1), (KK_ACC_FLOAT)(0.5*v1));
      Kokkos::atomic_add(&d_vatom(i,2), (KK_ACC_FLOAT)(0.5*v2));
      Kokkos::atomic_add(&d_vatom(i,3), (KK_ACC_FLOAT)(0.5*v3));
      Kokkos::atomic_add(&d_vatom(i,4), (KK_ACC_FLOAT)(0.5*v4));
      Kokkos::atomic_add(&d_vatom(i,5), (KK_ACC_FLOAT)(0.5*v5));
      Kokkos::atomic_add(&d_vatom(j,0), (KK_ACC_FLOAT)(0.5*v0));
      Kokkos::atomic_add(&d_vatom(j,1), (KK_ACC_FLOAT)(0.5*v1));
      Kokkos::atomic_add(&d_vatom(j,2), (KK_ACC_FLOAT)(0.5*v2));
      Kokkos::atomic_add(&d_vatom(j,3), (KK_ACC_FLOAT)(0.5*v3));
      Kokkos::atomic_add(&d_vatom(j,4), (KK_ACC_FLOAT)(0.5*v4));
      Kokkos::atomic_add(&d_vatom(j,5), (KK_ACC_FLOAT)(0.5*v5));
    }
  }
}

/* ----------------------------------------------------------------------
   main interaction kernel: LJ on atom positions + Coulomb on the M sites
------------------------------------------------------------------------- */

template<class DeviceType>
template<int EVFLAG>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairLJCutTIP4PCutKokkos<DeviceType>::operator()(TagPairLJCutTIP4PCutCompute<EVFLAG>,
                                                     const int &ii, EV_FLOAT &ev) const
{
  const int i = d_ilist(ii);
  const int itype = type(i);
  const KK_FLOAT qtmp = q(i);
  const KK_FLOAT xtmp = x(i,0);
  const KK_FLOAT ytmp = x(i,1);
  const KK_FLOAT ztmp = x(i,2);

  KK_FLOAT x1[3];
  int iH1 = -1, iH2 = -1;
  if (itype == m_typeO) {
    iH1 = d_hneigh(i,0);
    iH2 = d_hneigh(i,1);
    x1[0] = d_newsite(i,0); x1[1] = d_newsite(i,1); x1[2] = d_newsite(i,2);
  } else {
    x1[0] = xtmp; x1[1] = ytmp; x1[2] = ztmp;
  }

  const int jnum = d_numneigh(i);

  for (int jj = 0; jj < jnum; jj++) {
    int j = d_neighbors(i,jj);
    const KK_FLOAT factor_lj = special_lj[sbmask(j)];
    const KK_FLOAT factor_coul = special_coul[sbmask(j)];
    j &= NEIGHMASK;
    const int jtype = type(j);

    // atom-based separation (used for LJ and the coulomb pre-screen)
    const KK_FLOAT adelx = xtmp - x(j,0);
    const KK_FLOAT adely = ytmp - x(j,1);
    const KK_FLOAT adelz = ztmp - x(j,2);
    const KK_FLOAT arsq = adelx*adelx + adely*adely + adelz*adelz;

    // ---- Lennard-Jones (on real atom positions) ----
    if (arsq < d_cut_ljsq(itype,jtype)) {
      const KK_FLOAT r2inv = (KK_FLOAT)1.0 / arsq;
      const KK_FLOAT r6inv = r2inv*r2inv*r2inv;
      const KK_FLOAT forcelj = factor_lj * r6inv *
          (d_lj1(itype,jtype)*r6inv - d_lj2(itype,jtype)) * r2inv;
      Kokkos::atomic_add(&f(i,0), (KK_ACC_FLOAT)(adelx*forcelj));
      Kokkos::atomic_add(&f(i,1), (KK_ACC_FLOAT)(adely*forcelj));
      Kokkos::atomic_add(&f(i,2), (KK_ACC_FLOAT)(adelz*forcelj));
      Kokkos::atomic_add(&f(j,0), (KK_ACC_FLOAT)(-adelx*forcelj));
      Kokkos::atomic_add(&f(j,1), (KK_ACC_FLOAT)(-adely*forcelj));
      Kokkos::atomic_add(&f(j,2), (KK_ACC_FLOAT)(-adelz*forcelj));
      if (EVFLAG) {
        KK_FLOAT evdwl = 0.0;
        if (this->eflag)
          evdwl = factor_lj * (r6inv*(d_lj3(itype,jtype)*r6inv - d_lj4(itype,jtype))
                               - d_offset(itype,jtype));
        ev_tally(ev,i,j,evdwl,forcelj,adelx,adely,adelz);
      }
    }

    // ---- Coulomb (on the M sites) ----
    if (arsq >= m_cut_coulsqplus) continue;

    KK_FLOAT delx = adelx, dely = adely, delz = adelz, rsq = arsq;
    int jH1 = -1, jH2 = -1;
    if (itype == m_typeO || jtype == m_typeO) {
      KK_FLOAT x2[3];
      if (jtype == m_typeO) {
        jH1 = d_hneigh(j,0);
        jH2 = d_hneigh(j,1);
        x2[0] = d_newsite(j,0); x2[1] = d_newsite(j,1); x2[2] = d_newsite(j,2);
      } else {
        x2[0] = x(j,0); x2[1] = x(j,1); x2[2] = x(j,2);
      }
      delx = x1[0] - x2[0];
      dely = x1[1] - x2[1];
      delz = x1[2] - x2[2];
      rsq = delx*delx + dely*dely + delz*delz;
    }

    if (rsq >= m_cut_coulsq) continue;

    const KK_FLOAT r2inv = (KK_FLOAT)1.0 / rsq;
    const KK_FLOAT forcecoul = qqrd2e * qtmp * q(j) * Kokkos::sqrt(r2inv);
    const KK_FLOAT cforce = factor_coul * forcecoul * r2inv;

    int n = 0, key = 0, vlist[6];
    KK_FLOAT v[6] = {0,0,0,0,0,0};

    if (itype != m_typeO) {
      Kokkos::atomic_add(&f(i,0), (KK_ACC_FLOAT)(delx*cforce));
      Kokkos::atomic_add(&f(i,1), (KK_ACC_FLOAT)(dely*cforce));
      Kokkos::atomic_add(&f(i,2), (KK_ACC_FLOAT)(delz*cforce));
      if (EVFLAG && this->vflag) {
        v[0] = x(i,0)*delx*cforce; v[1] = x(i,1)*dely*cforce; v[2] = x(i,2)*delz*cforce;
        v[3] = x(i,0)*dely*cforce; v[4] = x(i,0)*delz*cforce; v[5] = x(i,1)*delz*cforce;
      }
      vlist[n++] = i;
    } else {
      key++;
      const KK_FLOAT fdx = delx*cforce, fdy = dely*cforce, fdz = delz*cforce;
      const KK_FLOAT fOx = fdx*(1.0-m_alpha), fOy = fdy*(1.0-m_alpha), fOz = fdz*(1.0-m_alpha);
      const KK_FLOAT fHx = 0.5*m_alpha*fdx, fHy = 0.5*m_alpha*fdy, fHz = 0.5*m_alpha*fdz;
      Kokkos::atomic_add(&f(i,0), (KK_ACC_FLOAT)fOx);
      Kokkos::atomic_add(&f(i,1), (KK_ACC_FLOAT)fOy);
      Kokkos::atomic_add(&f(i,2), (KK_ACC_FLOAT)fOz);
      Kokkos::atomic_add(&f(iH1,0), (KK_ACC_FLOAT)fHx);
      Kokkos::atomic_add(&f(iH1,1), (KK_ACC_FLOAT)fHy);
      Kokkos::atomic_add(&f(iH1,2), (KK_ACC_FLOAT)fHz);
      Kokkos::atomic_add(&f(iH2,0), (KK_ACC_FLOAT)fHx);
      Kokkos::atomic_add(&f(iH2,1), (KK_ACC_FLOAT)fHy);
      Kokkos::atomic_add(&f(iH2,2), (KK_ACC_FLOAT)fHz);
      if (EVFLAG && this->vflag) {
        v[0] = x(i,0)*fOx + x(iH1,0)*fHx + x(iH2,0)*fHx;
        v[1] = x(i,1)*fOy + x(iH1,1)*fHy + x(iH2,1)*fHy;
        v[2] = x(i,2)*fOz + x(iH1,2)*fHz + x(iH2,2)*fHz;
        v[3] = x(i,0)*fOy + x(iH1,0)*fHy + x(iH2,0)*fHy;
        v[4] = x(i,0)*fOz + x(iH1,0)*fHz + x(iH2,0)*fHz;
        v[5] = x(i,1)*fOz + x(iH1,1)*fHz + x(iH2,1)*fHz;
      }
      vlist[n++] = i; vlist[n++] = iH1; vlist[n++] = iH2;
    }

    if (jtype != m_typeO) {
      Kokkos::atomic_add(&f(j,0), (KK_ACC_FLOAT)(-delx*cforce));
      Kokkos::atomic_add(&f(j,1), (KK_ACC_FLOAT)(-dely*cforce));
      Kokkos::atomic_add(&f(j,2), (KK_ACC_FLOAT)(-delz*cforce));
      if (EVFLAG && this->vflag) {
        v[0] -= x(j,0)*delx*cforce; v[1] -= x(j,1)*dely*cforce; v[2] -= x(j,2)*delz*cforce;
        v[3] -= x(j,0)*dely*cforce; v[4] -= x(j,0)*delz*cforce; v[5] -= x(j,1)*delz*cforce;
      }
      vlist[n++] = j;
    } else {
      key += 2;
      const KK_FLOAT fdx = -delx*cforce, fdy = -dely*cforce, fdz = -delz*cforce;
      const KK_FLOAT fOx = fdx*(1.0-m_alpha), fOy = fdy*(1.0-m_alpha), fOz = fdz*(1.0-m_alpha);
      const KK_FLOAT fHx = 0.5*m_alpha*fdx, fHy = 0.5*m_alpha*fdy, fHz = 0.5*m_alpha*fdz;
      Kokkos::atomic_add(&f(j,0), (KK_ACC_FLOAT)fOx);
      Kokkos::atomic_add(&f(j,1), (KK_ACC_FLOAT)fOy);
      Kokkos::atomic_add(&f(j,2), (KK_ACC_FLOAT)fOz);
      Kokkos::atomic_add(&f(jH1,0), (KK_ACC_FLOAT)fHx);
      Kokkos::atomic_add(&f(jH1,1), (KK_ACC_FLOAT)fHy);
      Kokkos::atomic_add(&f(jH1,2), (KK_ACC_FLOAT)fHz);
      Kokkos::atomic_add(&f(jH2,0), (KK_ACC_FLOAT)fHx);
      Kokkos::atomic_add(&f(jH2,1), (KK_ACC_FLOAT)fHy);
      Kokkos::atomic_add(&f(jH2,2), (KK_ACC_FLOAT)fHz);
      if (EVFLAG && this->vflag) {
        v[0] += x(j,0)*fOx + x(jH1,0)*fHx + x(jH2,0)*fHx;
        v[1] += x(j,1)*fOy + x(jH1,1)*fHy + x(jH2,1)*fHy;
        v[2] += x(j,2)*fOz + x(jH1,2)*fHz + x(jH2,2)*fHz;
        v[3] += x(j,0)*fOy + x(jH1,0)*fHy + x(jH2,0)*fHy;
        v[4] += x(j,0)*fOz + x(jH1,0)*fHz + x(jH2,0)*fHz;
        v[5] += x(j,1)*fOz + x(jH1,1)*fHz + x(jH2,1)*fHz;
      }
      vlist[n++] = j; vlist[n++] = jH1; vlist[n++] = jH2;
    }

    KK_FLOAT ecoul = 0.0;
    if (EVFLAG && this->eflag) ecoul = factor_coul * qqrd2e * qtmp * q(j) * Kokkos::sqrt(r2inv);

    if (EVFLAG) ev_tally_tip4p(ev,key,vlist,v,ecoul);
  }
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class PairLJCutTIP4PCutKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairLJCutTIP4PCutKokkos<LMPHostType>;
#endif
}
