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

#include "pair_mesomem_dipole_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "error.h"
#include "force.h"
#include "kokkos.h"
#include "memory_kokkos.h"
#include "neigh_list_kokkos.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "respa.h"
#include "update.h"

#include <cmath>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairMesomemDipoleKokkos<DeviceType>::PairMesomemDipoleKokkos(LAMMPS *lmp) : PairMesomemDipole(lmp)
{
  respa_enable = 0;

  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | F_MASK | TORQUE_MASK | TYPE_MASK | MU_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | TORQUE_MASK | ENERGY_MASK | VIRIAL_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairMesomemDipoleKokkos<DeviceType>::~PairMesomemDipoleKokkos()
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
void PairMesomemDipoleKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
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

  x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  torque = atomKK->k_torque.view<DeviceType>();
  mu = atomKK->k_mu.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  nlocal = atom->nlocal;
  nall = atom->nlocal + atom->nghost;
  for (int m = 0; m < 4; ++m) special_lj[m] = static_cast<KK_FLOAT>(force->special_lj[m]);
  newton_pair = force->newton_pair;

  // get the neighbor list and neighbors used in operator()

  NeighListKokkos<DeviceType>* k_list = static_cast<NeighListKokkos<DeviceType>*>(list);
  d_numneigh = k_list->d_numneigh;
  d_neighbors = k_list->d_neighbors;
  d_ilist = k_list->d_ilist;
  inum = list->inum;

  // loop over neighbors of my atoms

  copymode = 1;

  EV_FLOAT ev;
  if (neighflag == HALF) {
    if (newton_pair) launch<HALF,1>(ev);
    else launch<HALF,0>(ev);
  } else if (neighflag == HALFTHREAD) {
    if (newton_pair) launch<HALFTHREAD,1>(ev);
    else launch<HALFTHREAD,0>(ev);
  } else if (neighflag == FULL) {
    if (newton_pair) launch<FULL,1>(ev);
    else launch<FULL,0>(ev);
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

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR>
void PairMesomemDipoleKokkos<DeviceType>::launch(EV_FLOAT &ev)
{
  if (evflag)
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairMesomemDipoleKernel<NEIGHFLAG,NEWTON_PAIR,1> >(0,inum),*this,ev);
  else
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairMesomemDipoleKernel<NEIGHFLAG,NEWTON_PAIR,0> >(0,inum),*this);
}

/* ----------------------------------------------------------------------
   same model as PairMesomemDipole::mesomem_analytic(): isotropic attraction
   plus tilt and splay terms between the two directors, switched off with a
   distance dependent weight function
------------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairMesomemDipoleKokkos<DeviceType>::operator()(TagPairMesomemDipoleKernel<NEIGHFLAG,NEWTON_PAIR,EVFLAG>, const int ii, EV_FLOAT &ev) const {

  // The f and torque arrays are atomic for Half/Thread neighbor style
  Kokkos::View<KK_ACC_FLOAT*[3], typename DAT::t_kkacc_1d_3::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > a_f = f;
  Kokkos::View<KK_ACC_FLOAT*[3], typename DAT::t_kkacc_1d_3::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > a_torque = torque;

  const int i = d_ilist[ii];
  const KK_FLOAT xtmp = x(i,0);
  const KK_FLOAT ytmp = x(i,1);
  const KK_FLOAT ztmp = x(i,2);
  const int itype = type(i);
  const int jnum = d_numneigh[i];

  // particles with a zero dipole moment only have the isotropic interaction

  const bool has_ni = mu(i,3) > static_cast<KK_FLOAT>(0.0);
  KK_FLOAT ni[3] = {0.0, 0.0, 0.0};
  if (has_ni) {
    const KK_FLOAT inv_mag = static_cast<KK_FLOAT>(1.0) / mu(i,3);
    ni[0] = mu(i,0) * inv_mag;
    ni[1] = mu(i,1) * inv_mag;
    ni[2] = mu(i,2) * inv_mag;
  }

  KK_ACC_FLOAT fx_i = 0.0;
  KK_ACC_FLOAT fy_i = 0.0;
  KK_ACC_FLOAT fz_i = 0.0;
  KK_ACC_FLOAT torquex_i = 0.0;
  KK_ACC_FLOAT torquey_i = 0.0;
  KK_ACC_FLOAT torquez_i = 0.0;

  for (int jj = 0; jj < jnum; jj++) {
    int j = d_neighbors(i,jj);
    const KK_FLOAT factor_lj = special_lj[sbmask(j)];
    j &= NEIGHMASK;

    const KK_FLOAT delx = xtmp - x(j,0);
    const KK_FLOAT dely = ytmp - x(j,1);
    const KK_FLOAT delz = ztmp - x(j,2);
    const int jtype = type(j);
    const KK_FLOAT rsq = delx*delx + dely*dely + delz*delz;

    if (rsq >= d_cutsq(itype,jtype)) continue;

    const params_mesomem &p = params(itype,jtype);
    const KK_FLOAT r = Kokkos::sqrt(rsq);
    const KK_FLOAT rinv = static_cast<KK_FLOAT>(1.0) / r;
    const KK_FLOAT rhat[3] = {delx * rinv, dely * rinv, delz * rinv};

    // isotropic part:
    // r < sigma: eps * [(sigma/r)^4 - 2 (sigma/r)^2]
    // r >= sigma: -eps * cos^(2 zeta)[pi/2 (r - sigma) / (r_c - sigma)]

    KK_FLOAT energy, fpair;
    if (r < p.sigma) {
      const KK_FLOAT sr2 = p.sigma * p.sigma / rsq;
      const KK_FLOAT sr4 = sr2 * sr2;
      energy = p.eps * (sr4 - static_cast<KK_FLOAT>(2.0) * sr2);
      fpair = static_cast<KK_FLOAT>(4.0) * p.eps * rinv * (sr4 - sr2);
    } else {
      const KK_FLOAT g = p.gscale * (r - p.sigma);
      const KK_FLOAT cosg = Kokkos::cos(g);
      const KK_FLOAT sing = Kokkos::sin(g);
      KK_FLOAT cospow;
      if (p.zpow < 0) {
        cospow = Kokkos::pow(cosg, static_cast<KK_FLOAT>(2.0) * p.zeta - static_cast<KK_FLOAT>(1.0));
      } else {
        // integer power by repeated squaring, as MathSpecial::powint()
        cospow = static_cast<KK_FLOAT>(1.0);
        if (cosg == static_cast<KK_FLOAT>(0.0)) {
          if (p.zpow > 0) cospow = static_cast<KK_FLOAT>(0.0);
        } else {
          KK_FLOAT ww = cosg;
          for (int nn = p.zpow; nn != 0; nn >>= 1, ww *= ww)
            if (nn & 1) cospow *= ww;
        }
      }
      energy = -p.eps * cospow * cosg;
      fpair = static_cast<KK_FLOAT>(-2.0) * p.zeta * p.eps * cospow * sing * p.gscale;
    }

    KK_FLOAT fi[3] = {fpair * rhat[0], fpair * rhat[1], fpair * rhat[2]};
    KK_FLOAT ti[3] = {0.0, 0.0, 0.0};
    KK_FLOAT tj[3] = {0.0, 0.0, 0.0};

    // orientation dependent part, requires a dipole on both particles and r < w_c

    if (has_ni && (mu(j,3) > static_cast<KK_FLOAT>(0.0)) && (r < p.wc)) {

      // weight function w = exp[-r^2 / ((w_c/2)^2 (1 - (r/w_c)^4))]
      // it underflows to zero close to w_c, which also avoids a division by zero

      const KK_FLOAT rw = r * p.wc_inv;
      const KK_FLOAT rw4 = rw * rw * rw * rw;
      const KK_FLOAT dw = rw4 - static_cast<KK_FLOAT>(1.0);
      const KK_FLOAT wfac = rsq * p.wc_half2inv;
      const KK_FLOAT w = (dw < static_cast<KK_FLOAT>(0.0)) ? Kokkos::exp(wfac / dw) : static_cast<KK_FLOAT>(0.0);

      if (w > static_cast<KK_FLOAT>(0.0)) {
        const KK_FLOAT inv_mag = static_cast<KK_FLOAT>(1.0) / mu(j,3);
        const KK_FLOAT nj[3] = {mu(j,0) * inv_mag, mu(j,1) * inv_mag, mu(j,2) * inv_mag};

        const KK_FLOAT nirhat = ni[0]*rhat[0] + ni[1]*rhat[1] + ni[2]*rhat[2];
        const KK_FLOAT njrhat = nj[0]*rhat[0] + nj[1]*rhat[1] + nj[2]*rhat[2];
        const KK_FLOAT ninj = ni[0]*nj[0] + ni[1]*nj[1] + ni[2]*nj[2];

        // tilt: 1/2 k_tilt (d_i^2 + d_j^2) with d_i = n_i.rhat + s, d_j = n_j.rhat - s
        // splay: 1/2 k_splay (n_i.n_j - 1 + 2 s^2)^2
        // with the spontaneous curvature shift s = r C_0 / 2

        const KK_FLOAT s = static_cast<KK_FLOAT>(0.5) * r * p.c0;
        const KK_FLOAT di = nirhat + s;
        const KK_FLOAT dj = njrhat - s;
        const KK_FLOAT splay = ninj - static_cast<KK_FLOAT>(1.0) + static_cast<KK_FLOAT>(2.0) * s * s;
        const KK_FLOAT uang = static_cast<KK_FLOAT>(0.5) * p.ktilt * (di * di + dj * dj) +
          static_cast<KK_FLOAT>(0.5) * p.ksplay * splay * splay;

        // force from the angular dependence of the tilt term, with
        // d(n.rhat)/dr = (n - (n.rhat) rhat) / r

        const KK_FLOAT ftilt = -p.ktilt * rinv;
        for (int k = 0; k < 3; ++k)
          fi[k] += w * ftilt * (di * (ni[k] - nirhat * rhat[k]) + dj * (nj[k] - njrhat * rhat[k]));

        // radial forces from the weight function and from the r dependence of s

        const KK_FLOAT fradial = uang * static_cast<KK_FLOAT>(2.0) * wfac * (rw4 + static_cast<KK_FLOAT>(1.0)) * rinv / (dw * dw) +
          static_cast<KK_FLOAT>(0.5) * p.ktilt * p.c0 * (dj - di) - p.ksplay * splay * p.c0 * p.c0 * r;
        fi[0] += w * fradial * rhat[0];
        fi[1] += w * fradial * rhat[1];
        fi[2] += w * fradial * rhat[2];

        // torques: t = n x (-dU/dn)

        const KK_FLOAT rh_x_ni[3] = {rhat[1]*ni[2] - rhat[2]*ni[1], rhat[2]*ni[0] - rhat[0]*ni[2],
                                     rhat[0]*ni[1] - rhat[1]*ni[0]};
        const KK_FLOAT rh_x_nj[3] = {rhat[1]*nj[2] - rhat[2]*nj[1], rhat[2]*nj[0] - rhat[0]*nj[2],
                                     rhat[0]*nj[1] - rhat[1]*nj[0]};
        const KK_FLOAT ni_x_nj[3] = {ni[1]*nj[2] - ni[2]*nj[1], ni[2]*nj[0] - ni[0]*nj[2],
                                     ni[0]*nj[1] - ni[1]*nj[0]};

        for (int k = 0; k < 3; ++k) {
          ti[k] = w * (p.ktilt * di * rh_x_ni[k] - p.ksplay * splay * ni_x_nj[k]);
          tj[k] = w * (p.ktilt * dj * rh_x_nj[k] + p.ksplay * splay * ni_x_nj[k]);
        }

        energy += w * uang;
      }
    }

    if (factor_lj != static_cast<KK_FLOAT>(1.0)) {
      energy *= factor_lj;
      for (int k = 0; k < 3; ++k) {
        fi[k] *= factor_lj;
        ti[k] *= factor_lj;
        tj[k] *= factor_lj;
      }
    }

    fx_i += static_cast<KK_ACC_FLOAT>(fi[0]);
    fy_i += static_cast<KK_ACC_FLOAT>(fi[1]);
    fz_i += static_cast<KK_ACC_FLOAT>(fi[2]);
    torquex_i += static_cast<KK_ACC_FLOAT>(ti[0]);
    torquey_i += static_cast<KK_ACC_FLOAT>(ti[1]);
    torquez_i += static_cast<KK_ACC_FLOAT>(ti[2]);

    if ((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD) && (NEWTON_PAIR || j < nlocal)) {
      a_f(j,0) -= static_cast<KK_ACC_FLOAT>(fi[0]);
      a_f(j,1) -= static_cast<KK_ACC_FLOAT>(fi[1]);
      a_f(j,2) -= static_cast<KK_ACC_FLOAT>(fi[2]);
      a_torque(j,0) += static_cast<KK_ACC_FLOAT>(tj[0]);
      a_torque(j,1) += static_cast<KK_ACC_FLOAT>(tj[1]);
      a_torque(j,2) += static_cast<KK_ACC_FLOAT>(tj[2]);
    }

    if (EVFLAG) {
      if (eflag_global)
        ev.evdwl += (((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD)&&(NEWTON_PAIR||(j<nlocal)))?static_cast<KK_ACC_FLOAT>(1.0):static_cast<KK_ACC_FLOAT>(0.5))*static_cast<KK_ACC_FLOAT>(energy);

      // forces are not central, so the virial needs the force vector
      if (eflag_atom || vflag_either)
        ev_tally_xyz<NEIGHFLAG, NEWTON_PAIR>(ev, i, j, energy, fi[0], fi[1], fi[2], delx, dely, delz);
    }
  }

  a_f(i,0) += fx_i;
  a_f(i,1) += fy_i;
  a_f(i,2) += fz_i;
  a_torque(i,0) += torquex_i;
  a_torque(i,1) += torquey_i;
  a_torque(i,2) += torquez_i;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairMesomemDipoleKokkos<DeviceType>::operator()(TagPairMesomemDipoleKernel<NEIGHFLAG,NEWTON_PAIR,EVFLAG>, const int ii) const {
  EV_FLOAT ev;
  this->template operator()<NEIGHFLAG,NEWTON_PAIR,EVFLAG>(TagPairMesomemDipoleKernel<NEIGHFLAG,NEWTON_PAIR,EVFLAG>(), ii, ev);
}

template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairMesomemDipoleKokkos<DeviceType>::ev_tally_xyz(EV_FLOAT & ev, int i, int j, const KK_FLOAT &epair,
                                                        KK_FLOAT fx, KK_FLOAT fy, KK_FLOAT fz,
                                                        KK_FLOAT delx, KK_FLOAT dely, KK_FLOAT delz) const
{
  Kokkos::View<KK_ACC_FLOAT*, typename DAT::t_kkacc_1d::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > v_eatom = d_eatom;
  Kokkos::View<KK_ACC_FLOAT*[6], typename DAT::t_kkacc_1d_6::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > v_vatom = d_vatom;

  if (eflag_atom) {
    const KK_ACC_FLOAT epairhalf = static_cast<KK_ACC_FLOAT>(0.5) * static_cast<KK_ACC_FLOAT>(epair);
    if (NEIGHFLAG == FULL || newton_pair || i < nlocal)
      v_eatom[i] += epairhalf;
    if (NEIGHFLAG != FULL && (newton_pair || j < nlocal))
      v_eatom[j] += epairhalf;
  }

  if (vflag_either) {
    const KK_ACC_FLOAT v0 = static_cast<KK_ACC_FLOAT>(delx*fx);
    const KK_ACC_FLOAT v1 = static_cast<KK_ACC_FLOAT>(dely*fy);
    const KK_ACC_FLOAT v2 = static_cast<KK_ACC_FLOAT>(delz*fz);
    const KK_ACC_FLOAT v3 = static_cast<KK_ACC_FLOAT>(delx*fy);
    const KK_ACC_FLOAT v4 = static_cast<KK_ACC_FLOAT>(delx*fz);
    const KK_ACC_FLOAT v5 = static_cast<KK_ACC_FLOAT>(dely*fz);

    if (vflag_global) {
      if (NEIGHFLAG != FULL) {
        if (NEWTON_PAIR) { // neigh half, newton on
          ev.v[0] += v0;
          ev.v[1] += v1;
          ev.v[2] += v2;
          ev.v[3] += v3;
          ev.v[4] += v4;
          ev.v[5] += v5;
        } else { // neigh half, newton off
          if (i < nlocal) {
            ev.v[0] += static_cast<KK_ACC_FLOAT>(0.5)*v0;
            ev.v[1] += static_cast<KK_ACC_FLOAT>(0.5)*v1;
            ev.v[2] += static_cast<KK_ACC_FLOAT>(0.5)*v2;
            ev.v[3] += static_cast<KK_ACC_FLOAT>(0.5)*v3;
            ev.v[4] += static_cast<KK_ACC_FLOAT>(0.5)*v4;
            ev.v[5] += static_cast<KK_ACC_FLOAT>(0.5)*v5;
          }
          if (j < nlocal) {
            ev.v[0] += static_cast<KK_ACC_FLOAT>(0.5)*v0;
            ev.v[1] += static_cast<KK_ACC_FLOAT>(0.5)*v1;
            ev.v[2] += static_cast<KK_ACC_FLOAT>(0.5)*v2;
            ev.v[3] += static_cast<KK_ACC_FLOAT>(0.5)*v3;
            ev.v[4] += static_cast<KK_ACC_FLOAT>(0.5)*v4;
            ev.v[5] += static_cast<KK_ACC_FLOAT>(0.5)*v5;
          }
        }
      } else { //neigh full
        ev.v[0] += static_cast<KK_ACC_FLOAT>(0.5)*v0;
        ev.v[1] += static_cast<KK_ACC_FLOAT>(0.5)*v1;
        ev.v[2] += static_cast<KK_ACC_FLOAT>(0.5)*v2;
        ev.v[3] += static_cast<KK_ACC_FLOAT>(0.5)*v3;
        ev.v[4] += static_cast<KK_ACC_FLOAT>(0.5)*v4;
        ev.v[5] += static_cast<KK_ACC_FLOAT>(0.5)*v5;
      }
    }

    if (vflag_atom) {

      if (NEIGHFLAG == FULL || NEWTON_PAIR || i < nlocal) {
        v_vatom(i,0) += static_cast<KK_ACC_FLOAT>(0.5)*v0;
        v_vatom(i,1) += static_cast<KK_ACC_FLOAT>(0.5)*v1;
        v_vatom(i,2) += static_cast<KK_ACC_FLOAT>(0.5)*v2;
        v_vatom(i,3) += static_cast<KK_ACC_FLOAT>(0.5)*v3;
        v_vatom(i,4) += static_cast<KK_ACC_FLOAT>(0.5)*v4;
        v_vatom(i,5) += static_cast<KK_ACC_FLOAT>(0.5)*v5;
      }
      if (NEIGHFLAG != FULL && (NEWTON_PAIR || j < nlocal)) {
        v_vatom(j,0) += static_cast<KK_ACC_FLOAT>(0.5)*v0;
        v_vatom(j,1) += static_cast<KK_ACC_FLOAT>(0.5)*v1;
        v_vatom(j,2) += static_cast<KK_ACC_FLOAT>(0.5)*v2;
        v_vatom(j,3) += static_cast<KK_ACC_FLOAT>(0.5)*v3;
        v_vatom(j,4) += static_cast<KK_ACC_FLOAT>(0.5)*v4;
        v_vatom(j,5) += static_cast<KK_ACC_FLOAT>(0.5)*v5;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

template<class DeviceType>
void PairMesomemDipoleKokkos<DeviceType>::allocate()
{
  PairMesomemDipole::allocate();

  int n = atom->ntypes;
  memory->destroy(cutsq);
  memoryKK->create_kokkos(k_cutsq,cutsq,n+1,n+1,"pair:cutsq");
  d_cutsq = k_cutsq.template view<DeviceType>();
  k_params = Kokkos::DualView<params_mesomem**,Kokkos::LayoutRight,DeviceType>("PairMesomemDipole::params",n+1,n+1);
  params = k_params.template view<DeviceType>();
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

template<class DeviceType>
void PairMesomemDipoleKokkos<DeviceType>::init_style()
{
  PairMesomemDipole::init_style();

  // error if rRESPA with inner levels

  if (update->whichflag == 1 && utils::strmatch(update->integrate_style,"^respa")) {
    int respa = 0;
    if (((Respa *) update->integrate)->level_inner >= 0) respa = 1;
    if (((Respa *) update->integrate)->level_middle >= 0) respa = 2;
    if (respa)
      error->all(FLERR,"Cannot use Kokkos pair style with rRESPA inner/middle");
  }

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
double PairMesomemDipoleKokkos<DeviceType>::init_one(int i, int j)
{
  double cutone = PairMesomemDipole::init_one(i,j);

  k_params.view_host()(i,j).sigma = static_cast<KK_FLOAT>(sigma[i][j]);
  k_params.view_host()(i,j).eps = static_cast<KK_FLOAT>(eps[i][j]);
  k_params.view_host()(i,j).zeta = static_cast<KK_FLOAT>(zeta[i][j]);
  k_params.view_host()(i,j).gscale = static_cast<KK_FLOAT>(gscale[i][j]);
  k_params.view_host()(i,j).ktilt = static_cast<KK_FLOAT>(ktilt[i][j]);
  k_params.view_host()(i,j).ksplay = static_cast<KK_FLOAT>(ksplay[i][j]);
  k_params.view_host()(i,j).c0 = static_cast<KK_FLOAT>(c0[i][j]);
  k_params.view_host()(i,j).wc = static_cast<KK_FLOAT>(weight_rcut[i][j]);
  k_params.view_host()(i,j).wc_inv = static_cast<KK_FLOAT>(wc_inv[i][j]);
  k_params.view_host()(i,j).wc_half2inv = static_cast<KK_FLOAT>(wc_half2inv[i][j]);
  k_params.view_host()(i,j).zpow = zpow[i][j];
  k_params.view_host()(j,i) = k_params.view_host()(i,j);
  k_params.modify_host();

  k_cutsq.view_host()(i,j) = k_cutsq.view_host()(j,i) = cutone*cutone;
  k_cutsq.modify_host();

  return cutone;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
int PairMesomemDipoleKokkos<DeviceType>::sbmask(const int& j) const {
  return j >> SBBITS & 3;
}

namespace LAMMPS_NS {
template class PairMesomemDipoleKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairMesomemDipoleKokkos<LMPHostType>;
#endif
}
