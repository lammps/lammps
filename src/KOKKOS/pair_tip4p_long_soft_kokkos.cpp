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

/* ----------------------------------------------------------------------
   Contributing author: Claude Opus 5.5 (Anthropic), under the direction of
   Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "pair_tip4p_long_soft_kokkos.h"

#include "atom_kokkos.h"
#include "kspace.h"

#include <cmath>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairTIP4PLongSoftKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
  const int inum = this->prepare(eflag_in, vflag_in);
  this->prepare_coul_long();

  // copy the per-type-pair soft-core coefficients to the device when they
  // have changed, i.e. after init_one() calls from setup or fix adapt

  if (!soft_coeffs_synced) {
    const int ntp1 = this->atom->ntypes + 1;
    if ((int)d_lam1.extent(0) != ntp1) {
      d_lam1 = typename AT::t_kkfloat_2d("tip4p/long/soft/kk:lam1",ntp1,ntp1);
      d_lam2 = typename AT::t_kkfloat_2d("tip4p/long/soft/kk:lam2",ntp1,ntp1);
    }
    auto h_lam1 = Kokkos::create_mirror_view(d_lam1);
    auto h_lam2 = Kokkos::create_mirror_view(d_lam2);
    for (int i = 1; i < ntp1; i++)
      for (int j = 1; j < ntp1; j++) {
        h_lam1(i,j) = static_cast<KK_FLOAT>(this->lam1[i][j]);
        h_lam2(i,j) = static_cast<KK_FLOAT>(this->lam2[i][j]);
      }
    Kokkos::deep_copy(d_lam1,h_lam1);
    Kokkos::deep_copy(d_lam2,h_lam2);
    soft_coeffs_synced = 1;
  }

  this->copymode = 1;

  // M (virtual charge) site of every O atom (local + ghost)
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagPairTIP4PNewsite>(0,this->nall), *this);

  EV_FLOAT ev;
  if (this->eflag || this->vflag)
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType,TagPairTIP4PLongSoftCompute<1>>(0,inum), *this, ev);
  else
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagPairTIP4PLongSoftCompute<0>>(0,inum), *this);

  this->copymode = 0;

  this->finalize(ev);
}

/* ----------------------------------------------------------------------
   mark the device copy of the soft-core coefficients as stale
------------------------------------------------------------------------- */

template<class DeviceType>
double PairTIP4PLongSoftKokkos<DeviceType>::init_one(int i, int j)
{
  const double cutone = PairTIP4PLongSoft::init_one(i,j);
  soft_coeffs_synced = 0;
  return cutone;
}

/* ----------------------------------------------------------------------
   main interaction kernel: soft-core long-range (Ewald) Coulomb on the M sites
------------------------------------------------------------------------- */

template<class DeviceType>
template<int EVFLAG>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairTIP4PLongSoftKokkos<DeviceType>::operator()(TagPairTIP4PLongSoftCompute<EVFLAG>,
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
    if (iH1 < 0) { d_h_missing() = 1; return; }
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
        if (jH1 < 0) { d_h_missing() = 1; continue; }
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

    // soft-core long-range (Ewald) Coulomb force and energy on the M-site distance
    KK_FLOAT ecoul = 0.0;
    const KK_FLOAT cforce = coul_long_soft(rsq, qtmp, q(j), factor_coul, d_lam1(itype,jtype),
                                           d_lam2(itype,jtype), EVFLAG && this->eflag, ecoul);

    int n = 0, key = 0, vlist[6];
    KK_ACC_FLOAT v[6] = {0,0,0,0,0,0};
    const bool do_virial = EVFLAG && this->vflag;

    apply_site_force(i,iH1,iH2,itype == m_typeO,delx,dely,delz, cforce,do_virial,1,n,key,vlist,v);
    apply_site_force(j,jH1,jH2,jtype == m_typeO,delx,dely,delz,-cforce,do_virial,2,n,key,vlist,v);

    if (EVFLAG) ev_tally_tip4p(ev,key,vlist,v,ecoul);
  }
}

template<class DeviceType>
template<int EVFLAG>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairTIP4PLongSoftKokkos<DeviceType>::operator()(TagPairTIP4PLongSoftCompute<EVFLAG>,
                                                     const int &ii) const
{
  EV_FLOAT ev;
  this->operator()(TagPairTIP4PLongSoftCompute<EVFLAG>(), ii, ev);
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class PairTIP4PLongSoftKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairTIP4PLongSoftKokkos<LMPHostType>;
#endif
}
