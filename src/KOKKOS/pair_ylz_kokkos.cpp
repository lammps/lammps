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
   Contributing author:  Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "pair_ylz_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "atom_vec_ellipsoid_kokkos.h"
#include "error.h"
#include "force.h"
#include "kokkos.h"
#include "math_const.h"
#include "math_extra_kokkos.h"
#include "memory_kokkos.h"
#include "neigh_list_kokkos.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "update.h"

#include <cmath>

using namespace LAMMPS_NS;
using MathConst::MY_PI2;
using MathConst::MY_PI;
using MathConst::MY_TWOBYSIXTH;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairYLZKokkos<DeviceType>::PairYLZKokkos(LAMMPS *lmp) : PairYLZ(lmp), avecKK(nullptr)
{
  respa_enable = 0;

  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | F_MASK | TORQUE_MASK | TYPE_MASK | ELLIPSOID_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | TORQUE_MASK | ENERGY_MASK | VIRIAL_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairYLZKokkos<DeviceType>::~PairYLZKokkos()
{
  if (copymode) return;

  if (allocated) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->destroy_kokkos(k_vatom,vatom);
    memoryKK->destroy_kokkos(k_cutsq,cutsq);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairYLZKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
  eflag = eflag_in;
  vflag = vflag_in;

  if (neighflag == FULL) no_virial_fdotr_compute = 1;

  ev_init(eflag,vflag,0);

  // reallocate per-atom arrays if necessary

  if (eflag_atom) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->create_kokkos(k_eatom,eatom,maxeatom,"pair:eatom");
    d_eatom = k_eatom.view<DeviceType>();
  }
  if (vflag_atom) {
    memoryKK->destroy_kokkos(k_vatom,vatom);
    memoryKK->create_kokkos(k_vatom,vatom,maxvatom,"pair:vatom");
    d_vatom = k_vatom.view<DeviceType>();
  }

  atomKK->sync(execution_space,datamask_read);
  k_cutsq.template sync<DeviceType>();
  k_params.template sync<DeviceType>();
  if (eflag || vflag) atomKK->modified(execution_space,datamask_modify);
  else atomKK->modified(execution_space,F_MASK | TORQUE_MASK);

  x           = atomKK->k_x.view<DeviceType>();
  f           = atomKK->k_f.view<DeviceType>();
  torque      = atomKK->k_torque.view<DeviceType>();
  type        = atomKK->k_type.view<DeviceType>();
  d_ellipsoid = atomKK->k_ellipsoid.view<DeviceType>();
  bonus       = avecKK->k_bonus.view<DeviceType>();
  nlocal      = atom->nlocal;
  nall        = atom->nlocal + atom->nghost;
  newton_pair = force->newton_pair;
  special_lj[0] = static_cast<KK_FLOAT>(force->special_lj[0]);
  special_lj[1] = static_cast<KK_FLOAT>(force->special_lj[1]);
  special_lj[2] = static_cast<KK_FLOAT>(force->special_lj[2]);
  special_lj[3] = static_cast<KK_FLOAT>(force->special_lj[3]);

  NeighListKokkos<DeviceType>* k_list = static_cast<NeighListKokkos<DeviceType>*>(list);
  d_numneigh  = k_list->d_numneigh;
  d_neighbors = k_list->d_neighbors;
  d_ilist     = k_list->d_ilist;
  int inum    = list->inum;

  copymode = 1;

  EV_FLOAT ev;

  if (atom->ntypes > MAX_TYPES_STACKPARAMS) {
    if (evflag) {
      if (neighflag == HALF) {
        if (newton_pair) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType,TagPairYLZKernel<HALF,1,1,false>>(0,inum),*this,ev);
        else             Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType,TagPairYLZKernel<HALF,0,1,false>>(0,inum),*this,ev);
      } else if (neighflag == HALFTHREAD) {
        if (newton_pair) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType,TagPairYLZKernel<HALFTHREAD,1,1,false>>(0,inum),*this,ev);
        else             Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType,TagPairYLZKernel<HALFTHREAD,0,1,false>>(0,inum),*this,ev);
      } else if (neighflag == FULL) {
        if (newton_pair) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType,TagPairYLZKernel<FULL,1,1,false>>(0,inum),*this,ev);
        else             Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType,TagPairYLZKernel<FULL,0,1,false>>(0,inum),*this,ev);
      }
    } else {
      if (neighflag == HALF) {
        if (newton_pair) Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagPairYLZKernel<HALF,1,0,false>>(0,inum),*this);
        else             Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagPairYLZKernel<HALF,0,0,false>>(0,inum),*this);
      } else if (neighflag == HALFTHREAD) {
        if (newton_pair) Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagPairYLZKernel<HALFTHREAD,1,0,false>>(0,inum),*this);
        else             Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagPairYLZKernel<HALFTHREAD,0,0,false>>(0,inum),*this);
      } else if (neighflag == FULL) {
        if (newton_pair) Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagPairYLZKernel<FULL,1,0,false>>(0,inum),*this);
        else             Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagPairYLZKernel<FULL,0,0,false>>(0,inum),*this);
      }
    }
  } else {
    if (evflag) {
      if (neighflag == HALF) {
        if (newton_pair) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType,TagPairYLZKernel<HALF,1,1,true>>(0,inum),*this,ev);
        else             Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType,TagPairYLZKernel<HALF,0,1,true>>(0,inum),*this,ev);
      } else if (neighflag == HALFTHREAD) {
        if (newton_pair) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType,TagPairYLZKernel<HALFTHREAD,1,1,true>>(0,inum),*this,ev);
        else             Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType,TagPairYLZKernel<HALFTHREAD,0,1,true>>(0,inum),*this,ev);
      } else if (neighflag == FULL) {
        if (newton_pair) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType,TagPairYLZKernel<FULL,1,1,true>>(0,inum),*this,ev);
        else             Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType,TagPairYLZKernel<FULL,0,1,true>>(0,inum),*this,ev);
      }
    } else {
      if (neighflag == HALF) {
        if (newton_pair) Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagPairYLZKernel<HALF,1,0,true>>(0,inum),*this);
        else             Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagPairYLZKernel<HALF,0,0,true>>(0,inum),*this);
      } else if (neighflag == HALFTHREAD) {
        if (newton_pair) Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagPairYLZKernel<HALFTHREAD,1,0,true>>(0,inum),*this);
        else             Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagPairYLZKernel<HALFTHREAD,0,0,true>>(0,inum),*this);
      } else if (neighflag == FULL) {
        if (newton_pair) Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagPairYLZKernel<FULL,1,0,true>>(0,inum),*this);
        else             Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagPairYLZKernel<FULL,0,0,true>>(0,inum),*this);
      }
    }
  }

  if (eflag_global) eng_vdwl += static_cast<double>(ev.evdwl);
  if (vflag_global) {
    virial[0] += static_cast<double>(ev.v[0]);
    virial[1] += static_cast<double>(ev.v[1]);
    virial[2] += static_cast<double>(ev.v[2]);
    virial[3] += static_cast<double>(ev.v[3]);
    virial[4] += static_cast<double>(ev.v[4]);
    virial[5] += static_cast<double>(ev.v[5]);
  }

  if (eflag_atom) {
    k_eatom.template modify<DeviceType>();
    k_eatom.sync_host();
  }
  if (vflag_atom) {
    k_vatom.template modify<DeviceType>();
    k_vatom.sync_host();
  }

  if (vflag_fdotr) pair_virial_fdotr_compute(this);

  copymode = 0;
}

/* ----------------------------------------------------------------------
   Pairwise YLZ kernel: compute forces and torques for one central atom i
------------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG, bool STACKPARAMS>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairYLZKokkos<DeviceType>::operator()(TagPairYLZKernel<NEIGHFLAG,NEWTON_PAIR,EVFLAG,STACKPARAMS>, const int ii, EV_FLOAT &ev) const
{
  // atomic views for half-neighbor accumulation
  Kokkos::View<KK_ACC_FLOAT*[3],typename DAT::t_kkacc_1d_3::array_layout,
               typename KKDevice<DeviceType>::value,
               Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value>> a_f = f;
  Kokkos::View<KK_ACC_FLOAT*[3],typename DAT::t_kkacc_1d_3::array_layout,
               typename KKDevice<DeviceType>::value,
               Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value>> a_torque = torque;

  const int i = d_ilist[ii];
  const int itype = type(i);
  const int ibnum = d_ellipsoid(i);
  if (ibnum < 0) return;   // skip non-ellipsoid atoms

  const KK_FLOAT xtmp = x(i,0);
  const KK_FLOAT ytmp = x(i,1);
  const KK_FLOAT ztmp = x(i,2);

  // rotation matrix for atom i: quat_to_mat gives R, we extract
  // the first row of R^T = first column of R, i.e., ni1 = R[:,0]
  const auto *iquat = bonus(ibnum).quat;
  KK_FLOAT ai[3][3];
  MathExtraKokkos::quat_to_mat(iquat, ai);
  // ni1[k] = ai[k][0]  (first column of rotation matrix)

  const int jnum = d_numneigh[i];

  KK_ACC_FLOAT fx_i    = 0.0;
  KK_ACC_FLOAT fy_i    = 0.0;
  KK_ACC_FLOAT fz_i    = 0.0;
  KK_ACC_FLOAT torx_i  = 0.0;
  KK_ACC_FLOAT tory_i  = 0.0;
  KK_ACC_FLOAT torz_i  = 0.0;

  for (int jj = 0; jj < jnum; jj++) {
    int j = d_neighbors(i,jj);
    const KK_FLOAT factor_lj = special_lj[sbmask(j)];
    j &= NEIGHMASK;

    const int jbnum = d_ellipsoid(j);
    if (jbnum < 0) continue;   // skip non-ellipsoid atoms

    const KK_FLOAT r12x = x(j,0) - xtmp;
    const KK_FLOAT r12y = x(j,1) - ytmp;
    const KK_FLOAT r12z = x(j,2) - ztmp;
    const KK_FLOAT rsq  = r12x*r12x + r12y*r12y + r12z*r12z;
    const int jtype = type(j);

    KK_FLOAT cutsq_ij = STACKPARAMS ? m_cutsq[itype][jtype] : d_cutsq(itype,jtype);
    if (rsq >= cutsq_ij) continue;

    // pair parameters
    const KK_FLOAT energy_well = STACKPARAMS ? m_params[itype][jtype].epsilon : d_params(itype,jtype).epsilon;
    const KK_FLOAT sig         = STACKPARAMS ? m_params[itype][jtype].sigma   : d_params(itype,jtype).sigma;
    const KK_FLOAT zt          = STACKPARAMS ? m_params[itype][jtype].zeta    : d_params(itype,jtype).zeta;
    const KK_FLOAT muu         = STACKPARAMS ? m_params[itype][jtype].mu      : d_params(itype,jtype).mu;
    const KK_FLOAT sint        = STACKPARAMS ? m_params[itype][jtype].beta    : d_params(itype,jtype).beta;
    const KK_FLOAT rmin        = MY_TWOBYSIXTH * sig;
    const KK_FLOAT rcut        = Kokkos::sqrt(cutsq_ij);

    // rotation matrix for atom j
    const auto *jquat = bonus(jbnum).quat;
    KK_FLOAT aj[3][3];
    MathExtraKokkos::quat_to_mat(jquat, aj);

    // orientation axes: ni1 = first column of ai, nj1 = first column of aj
    const KK_FLOAT ni1x = ai[0][0];
    const KK_FLOAT ni1y = ai[1][0];
    const KK_FLOAT ni1z = ai[2][0];
    const KK_FLOAT nj1x = aj[0][0];
    const KK_FLOAT nj1y = aj[1][0];
    const KK_FLOAT nj1z = aj[2][0];

    const KK_FLOAT r = Kokkos::sqrt(rsq);
    // r12hat = unit vector from i to j
    const KK_FLOAT inv_r = static_cast<KK_FLOAT>(1.0)/r;
    const KK_FLOAT r12hx = r12x*inv_r;
    const KK_FLOAT r12hy = r12y*inv_r;
    const KK_FLOAT r12hz = r12z*inv_r;

    const KK_FLOAT ninj     = ni1x*nj1x + ni1y*nj1y + ni1z*nj1z;
    const KK_FLOAT ni1rhat  = ni1x*r12hx + ni1y*r12hy + ni1z*r12hz;
    const KK_FLOAT nj1rhat  = nj1x*r12hx + nj1y*r12hy + nj1z*r12hz;

    const KK_FLOAT a_val = ninj + (sint - ni1rhat)*(sint + nj1rhat)
                           - static_cast<KK_FLOAT>(2.0)*sint*sint;
    const KK_FLOAT phi   = static_cast<KK_FLOAT>(1.0) + (a_val - static_cast<KK_FLOAT>(1.0))*muu;

    // dphi/drhat
    const KK_FLOAT dphi_drhat_x = muu*((sint - ni1rhat)*nj1x - ni1x*(sint + nj1rhat));
    const KK_FLOAT dphi_drhat_y = muu*((sint - ni1rhat)*nj1y - ni1y*(sint + nj1rhat));
    const KK_FLOAT dphi_drhat_z = muu*((sint - ni1rhat)*nj1z - ni1z*(sint + nj1rhat));

    // dphi/dni1
    const KK_FLOAT spnj = sint + nj1rhat;
    const KK_FLOAT dphi_dni1x = muu*(nj1x - r12hx*spnj);
    const KK_FLOAT dphi_dni1y = muu*(nj1y - r12hy*spnj);
    const KK_FLOAT dphi_dni1z = muu*(nj1z - r12hz*spnj);

    // dphi/dnj1
    const KK_FLOAT spni = sint - ni1rhat;
    const KK_FLOAT dphi_dnj1x = muu*(ni1x + r12hx*spni);
    const KK_FLOAT dphi_dnj1y = muu*(ni1y + r12hy*spni);
    const KK_FLOAT dphi_dnj1z = muu*(ni1z + r12hz*spni);

    // radial and orientation-derivative parts of the energy/force
    KK_FLOAT U, dUdr, dUdphi;
    if (r < rmin) {
      const KK_FLOAT t = rmin/r;
      const KK_FLOAT t2 = t*t;
      const KK_FLOAT t4 = t2*t2;
      const KK_FLOAT uR = (t4 - static_cast<KK_FLOAT>(2.0)*t2)*energy_well;
      U      = uR + (static_cast<KK_FLOAT>(1.0) - phi)*energy_well;
      dUdr   = static_cast<KK_FLOAT>(4.0)*(t2 - t4)/r*energy_well;
      dUdphi = -energy_well;
    } else {
      const KK_FLOAT t   = MY_PI2*(r - rmin)/(rcut - rmin);
      const KK_FLOAT cos_t = Kokkos::cos(t);
      // t1 = cos_t^(2*zt-1)
      KK_FLOAT t1 = cos_t;
      for (int k = 1; k <= static_cast<int>(static_cast<KK_FLOAT>(2.0)*zt) - 2; k++) t1 *= cos_t;
      const KK_FLOAT uA = -energy_well*t1*cos_t;
      U      = uA*phi;
      dUdr   = MY_PI*zt/(rcut - rmin)*t1*Kokkos::sin(t)*phi*energy_well;
      dUdphi = uA;
    }

    // dU/drhat
    const KK_FLOAT dUdrhat_x = dUdphi*dphi_drhat_x;
    const KK_FLOAT dUdrhat_y = dUdphi*dphi_drhat_y;
    const KK_FLOAT dUdrhat_z = dUdphi*dphi_drhat_z;
    const KK_FLOAT dUdrhatrhat = dUdrhat_x*r12hx + dUdrhat_y*r12hy + dUdrhat_z*r12hz;

    // force on i from j (note sign: fforce points from i toward j in ylz_analytic)
    const KK_FLOAT ff_x = factor_lj*(dUdr*r12hx + (dUdrhat_x - dUdrhatrhat*r12hx)*inv_r);
    const KK_FLOAT ff_y = factor_lj*(dUdr*r12hy + (dUdrhat_y - dUdrhatrhat*r12hy)*inv_r);
    const KK_FLOAT ff_z = factor_lj*(dUdr*r12hz + (dUdrhat_z - dUdrhatrhat*r12hz)*inv_r);

    // torque on i: cross(dU/dni1, ni1)  [sign swap: cross(ni1, dU/dni1) with sign flip]
    const KK_FLOAT dUdni1x = factor_lj*dUdphi*dphi_dni1x;
    const KK_FLOAT dUdni1y = factor_lj*dUdphi*dphi_dni1y;
    const KK_FLOAT dUdni1z = factor_lj*dUdphi*dphi_dni1z;
    // ttor = cross(dUdni1, ni1) = dUdni1 x ni1
    const KK_FLOAT ttor_x = dUdni1y*ni1z - dUdni1z*ni1y;
    const KK_FLOAT ttor_y = dUdni1z*ni1x - dUdni1x*ni1z;
    const KK_FLOAT ttor_z = dUdni1x*ni1y - dUdni1y*ni1x;

    fx_i   += ff_x;
    fy_i   += ff_y;
    fz_i   += ff_z;
    torx_i += ttor_x;
    tory_i += ttor_y;
    torz_i += ttor_z;

    if ((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD) && (NEWTON_PAIR || j < nlocal)) {
      // torque on j: cross(dU/dnj1, nj1)
      const KK_FLOAT dUdnj1x = factor_lj*dUdphi*dphi_dnj1x;
      const KK_FLOAT dUdnj1y = factor_lj*dUdphi*dphi_dnj1y;
      const KK_FLOAT dUdnj1z = factor_lj*dUdphi*dphi_dnj1z;
      const KK_FLOAT rtor_x = dUdnj1y*nj1z - dUdnj1z*nj1y;
      const KK_FLOAT rtor_y = dUdnj1z*nj1x - dUdnj1x*nj1z;
      const KK_FLOAT rtor_z = dUdnj1x*nj1y - dUdnj1y*nj1x;

      a_f(j,0) -= ff_x;
      a_f(j,1) -= ff_y;
      a_f(j,2) -= ff_z;
      a_torque(j,0) += rtor_x;
      a_torque(j,1) += rtor_y;
      a_torque(j,2) += rtor_z;
    }

    if (EVFLAG) {
      const KK_FLOAT evdwl = factor_lj*U;
      ev_tally_xyz<NEIGHFLAG,NEWTON_PAIR>(ev, i, j, evdwl, ff_x, ff_y, ff_z, -r12x, -r12y, -r12z);
    }
  }

  a_f(i,0) += fx_i;
  a_f(i,1) += fy_i;
  a_f(i,2) += fz_i;
  a_torque(i,0) += torx_i;
  a_torque(i,1) += tory_i;
  a_torque(i,2) += torz_i;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG, bool STACKPARAMS>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairYLZKokkos<DeviceType>::operator()(TagPairYLZKernel<NEIGHFLAG,NEWTON_PAIR,EVFLAG,STACKPARAMS>, const int ii) const
{
  EV_FLOAT ev;
  this->template operator()<NEIGHFLAG,NEWTON_PAIR,EVFLAG,STACKPARAMS>(TagPairYLZKernel<NEIGHFLAG,NEWTON_PAIR,EVFLAG,STACKPARAMS>(), ii, ev);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairYLZKokkos<DeviceType>::ev_tally_xyz(EV_FLOAT &ev, int i, int j, const KK_FLOAT &epair,
                                             KK_FLOAT fx, KK_FLOAT fy, KK_FLOAT fz,
                                             KK_FLOAT delx, KK_FLOAT dely, KK_FLOAT delz) const
{
  const KK_FLOAT half = static_cast<KK_FLOAT>(0.5);
  const KK_FLOAT efactor = ((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD)&&(NEWTON_PAIR||(j<nlocal))) ?
      static_cast<KK_FLOAT>(1.0) : half;

  if (eflag_atom) {
    const KK_FLOAT epairhalf = half*epair;
    Kokkos::atomic_add(&d_eatom[i], epairhalf);
    if ((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD) && (NEWTON_PAIR || j < nlocal))
      Kokkos::atomic_add(&d_eatom[j], epairhalf);
  }

  if (vflag_either) {
    const KK_FLOAT v0 = delx*fx;
    const KK_FLOAT v1 = dely*fy;
    const KK_FLOAT v2 = delz*fz;
    const KK_FLOAT v3 = delx*fy;
    const KK_FLOAT v4 = delx*fz;
    const KK_FLOAT v5 = dely*fz;

    if (vflag_global) {
      if ((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD) && (NEWTON_PAIR || j < nlocal)) {
        ev.v[0] += v0; ev.v[1] += v1; ev.v[2] += v2;
        ev.v[3] += v3; ev.v[4] += v4; ev.v[5] += v5;
      } else {
        ev.v[0] += half*v0; ev.v[1] += half*v1; ev.v[2] += half*v2;
        ev.v[3] += half*v3; ev.v[4] += half*v4; ev.v[5] += half*v5;
      }
    }

    if (vflag_atom) {
      Kokkos::atomic_add(&d_vatom(i,0), half*v0);
      Kokkos::atomic_add(&d_vatom(i,1), half*v1);
      Kokkos::atomic_add(&d_vatom(i,2), half*v2);
      Kokkos::atomic_add(&d_vatom(i,3), half*v3);
      Kokkos::atomic_add(&d_vatom(i,4), half*v4);
      Kokkos::atomic_add(&d_vatom(i,5), half*v5);
      if ((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD) && (NEWTON_PAIR || j < nlocal)) {
        Kokkos::atomic_add(&d_vatom(j,0), half*v0);
        Kokkos::atomic_add(&d_vatom(j,1), half*v1);
        Kokkos::atomic_add(&d_vatom(j,2), half*v2);
        Kokkos::atomic_add(&d_vatom(j,3), half*v3);
        Kokkos::atomic_add(&d_vatom(j,4), half*v4);
        Kokkos::atomic_add(&d_vatom(j,5), half*v5);
      }
    }
  }

  if (eflag_global)
    ev.evdwl += efactor*epair;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
int PairYLZKokkos<DeviceType>::sbmask(const int &j) const
{
  return j >> SBBITS & 3;
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

template<class DeviceType>
void PairYLZKokkos<DeviceType>::allocate()
{
  PairYLZ::allocate();

  int n = atom->ntypes;
  memory->destroy(cutsq);
  memoryKK->create_kokkos(k_cutsq,cutsq,n+1,n+1,"pair:cutsq");
  d_cutsq = k_cutsq.template view<DeviceType>();
  k_params = Kokkos::DualView<params_ylz**,Kokkos::LayoutRight,DeviceType>("PairYLZ::params",n+1,n+1);
  d_params = k_params.template view<DeviceType>();
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

template<class DeviceType>
void PairYLZKokkos<DeviceType>::init_style()
{
  PairYLZ::init_style();

  avecKK = dynamic_cast<AtomVecEllipsoidKokkos *>(atom->style_match("ellipsoid"));
  if (!avecKK) error->all(FLERR,"Pair style ylz/kk requires atom style ellipsoid/kk");

  // adjust neighbor list request for KOKKOS

  neighflag = lmp->kokkos->neighflag;
  auto request = neighbor->find_request(this);
  request->set_kokkos_host(std::is_same_v<DeviceType,LMPHostType> &&
                           !std::is_same_v<DeviceType,LMPDeviceType>);
  request->set_kokkos_device(std::is_same_v<DeviceType,LMPDeviceType>);
  if (neighflag == FULL) request->enable_full();
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

template<class DeviceType>
double PairYLZKokkos<DeviceType>::init_one(int i, int j)
{
  double cutone = PairYLZ::init_one(i,j);

  k_params.view_host()(i,j).epsilon = static_cast<KK_FLOAT>(epsilon[i][j]);
  k_params.view_host()(i,j).sigma   = static_cast<KK_FLOAT>(sigma[i][j]);
  k_params.view_host()(i,j).zeta    = static_cast<KK_FLOAT>(zeta[i][j]);
  k_params.view_host()(i,j).mu      = static_cast<KK_FLOAT>(mu[i][j]);
  k_params.view_host()(i,j).beta    = static_cast<KK_FLOAT>(beta[i][j]);
  k_params.view_host()(i,j).cutsq   = static_cast<KK_FLOAT>(cutone*cutone);
  k_params.view_host()(j,i) = k_params.view_host()(i,j);
  if (i<MAX_TYPES_STACKPARAMS+1 && j<MAX_TYPES_STACKPARAMS+1) {
    m_params[i][j] = m_params[j][i] = k_params.view_host()(i,j);
    m_cutsq[j][i] = m_cutsq[i][j] = static_cast<KK_FLOAT>(cutone*cutone);
  }

  k_cutsq.view_host()(i,j) = k_cutsq.view_host()(j,i) = cutone*cutone;
  k_cutsq.modify_host();
  k_params.modify_host();

  return cutone;
}

namespace LAMMPS_NS {
template class PairYLZKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairYLZKokkos<LMPHostType>;
#endif
}
