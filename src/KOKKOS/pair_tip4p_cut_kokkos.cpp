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

#include "pair_tip4p_cut_kokkos.h"

#include <cmath>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairTIP4PCutKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
  const int inum = this->prepare(eflag_in, vflag_in);

  this->copymode = 1;

  // M (virtual charge) site of every O atom (local + ghost)
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagPairTIP4PNewsite>(0,this->nall), *this);

  EV_FLOAT ev;
  if (this->eflag || this->vflag)
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType,TagPairTIP4PCutCompute<1>>(0,inum), *this, ev);
  else
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType,TagPairTIP4PCutCompute<0>>(0,inum), *this, ev);

  this->copymode = 0;

  this->finalize(ev);
}

/* ----------------------------------------------------------------------
   1/r cut Coulomb on the M sites, force redistributed onto O and its two H
------------------------------------------------------------------------- */

template<class DeviceType>
template<int EVFLAG>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairTIP4PCutKokkos<DeviceType>::operator()(TagPairTIP4PCutCompute<EVFLAG>,
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
    const KK_FLOAT factor_coul = special_coul[sbmask(j)];
    j &= NEIGHMASK;
    const int jtype = type(j);

    KK_FLOAT delx = xtmp - x(j,0);
    KK_FLOAT dely = ytmp - x(j,1);
    KK_FLOAT delz = ztmp - x(j,2);
    KK_FLOAT rsq = delx*delx + dely*dely + delz*delz;

    if (rsq >= m_cut_coulsqplus) continue;

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
template class PairTIP4PCutKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairTIP4PCutKokkos<LMPHostType>;
#endif
}
