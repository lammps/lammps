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
   Contributing author: Stan Moore (SNL)
------------------------------------------------------------------------- */

#include "pair_colloid_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "error.h"
#include "force.h"
#include "kokkos.h"
#include "math_special_kokkos.h"
#include "memory_kokkos.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "respa.h"
#include "update.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace MathSpecialKokkos;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairColloidKokkos<DeviceType>::PairColloidKokkos(LAMMPS *lmp) : PairColloid(lmp)
{
  respa_enable = 0;

  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | F_MASK | TYPE_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | ENERGY_MASK | VIRIAL_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairColloidKokkos<DeviceType>::~PairColloidKokkos()
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
void PairColloidKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
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
  else atomKK->modified(execution_space,F_MASK);

  x = atomKK->k_x.view<DeviceType>();
  c_x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  nlocal = atom->nlocal;
  nall = atom->nlocal + atom->nghost;
  newton_pair = force->newton_pair;
  special_lj[0] = force->special_lj[0];
  special_lj[1] = force->special_lj[1];
  special_lj[2] = force->special_lj[2];
  special_lj[3] = force->special_lj[3];

  // loop over neighbors of my atoms

  copymode = 1;
  EV_FLOAT ev = pair_compute<PairColloidKokkos<DeviceType>,void >(this,(NeighListKokkos<DeviceType>*)list);

  if (eflag_global) eng_vdwl += ev.evdwl;
  if (vflag_global) {
    virial[0] += ev.v[0];
    virial[1] += ev.v[1];
    virial[2] += ev.v[2];
    virial[3] += ev.v[3];
    virial[4] += ev.v[4];
    virial[5] += ev.v[5];
  }

  if (vflag_fdotr) pair_virial_fdotr_compute(this);

  if (eflag_atom) {
    k_eatom.template modify<DeviceType>();
    k_eatom.sync_host();
  }

  if (vflag_atom) {
    k_vatom.template modify<DeviceType>();
    k_vatom.sync_host();
  }

  copymode = 0;
}

template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairColloidKokkos<DeviceType>::
compute_fpair(const KK_FLOAT &rsq, const int &, const int &, const int &itype, const int &jtype) const {
  const int    iForm  = STACKPARAMS ? m_params[itype][jtype].form   : params(itype,jtype).form;
  const KK_FLOAT a12   = STACKPARAMS ? m_params[itype][jtype].a12    : params(itype,jtype).a12;
  const KK_FLOAT sigma3 = STACKPARAMS ? m_params[itype][jtype].sigma3 : params(itype,jtype).sigma3;
  const KK_FLOAT sigma6 = STACKPARAMS ? m_params[itype][jtype].sigma6 : params(itype,jtype).sigma6;
  const KK_FLOAT lj1   = STACKPARAMS ? m_params[itype][jtype].lj1   : params(itype,jtype).lj1;
  const KK_FLOAT lj2   = STACKPARAMS ? m_params[itype][jtype].lj2   : params(itype,jtype).lj2;
  const KK_FLOAT ppA1  = STACKPARAMS ? m_params[itype][jtype].a1    : params(itype,jtype).a1;
  const KK_FLOAT ppA2  = STACKPARAMS ? m_params[itype][jtype].a2    : params(itype,jtype).a2;

  if (iForm == SMALL_SMALL) {
    const KK_FLOAT r2inv = static_cast<KK_FLOAT>(1.0) / rsq;
    const KK_FLOAT r6inv = r2inv*r2inv*r2inv;
    return r6inv*(lj1*r6inv - lj2)*r2inv;

  } else if (iForm == SMALL_LARGE) {
    const KK_FLOAT c2  = ppA2;
    const KK_FLOAT K1  = c2*c2;
    const KK_FLOAT K2  = rsq;
    const KK_FLOAT K0  = K1 - rsq;
    const KK_FLOAT K4  = rsq*rsq;
    KK_FLOAT K3 = K1 - K2;
    K3 *= K3*K3;
    const KK_FLOAT K6  = K3*K3;
    const KK_FLOAT fR  = sigma3*a12*c2*K1/K3;
    return static_cast<KK_FLOAT>(4.0/15.0)*fR
            * (static_cast<KK_FLOAT>(2.0)*(K1+K2)
               * (K1*(static_cast<KK_FLOAT>(5.0)*K1+static_cast<KK_FLOAT>(22.0)*K2)
                  +static_cast<KK_FLOAT>(5.0)*K4)
               * sigma6/K6 - static_cast<KK_FLOAT>(5.0)) / K0;

  } else { // LARGE_LARGE
    const KK_FLOAT r  = sqrt(rsq);
    const KK_FLOAT c1 = ppA1;
    const KK_FLOAT c2 = ppA2;
    const KK_FLOAT K0 = c1*c2;
    const KK_FLOAT K1 = c1+c2;
    const KK_FLOAT K2 = c1-c2;
    const KK_FLOAT K3 = K1+r;
    const KK_FLOAT K4 = K1-r;
    const KK_FLOAT K5 = K2+r;
    const KK_FLOAT K6 = K2-r;
    const KK_FLOAT K7 = static_cast<KK_FLOAT>(1.0)/(K3*K4);
    const KK_FLOAT K8 = static_cast<KK_FLOAT>(1.0)/(K5*K6);
    const KK_FLOAT g0 = powint(K3,-7);
    const KK_FLOAT g1 = powint(K4,-7);
    const KK_FLOAT g2 = powint(K5,-7);
    const KK_FLOAT g3 = powint(K6,-7);
    const KK_FLOAT h0 = ((K3+static_cast<KK_FLOAT>(5.0)*K1)*K3+static_cast<KK_FLOAT>(30.0)*K0)*g0;
    const KK_FLOAT h1 = ((K4+static_cast<KK_FLOAT>(5.0)*K1)*K4+static_cast<KK_FLOAT>(30.0)*K0)*g1;
    const KK_FLOAT h2 = ((K5+static_cast<KK_FLOAT>(5.0)*K2)*K5-static_cast<KK_FLOAT>(30.0)*K0)*g2;
    const KK_FLOAT h3 = ((K6+static_cast<KK_FLOAT>(5.0)*K2)*K6-static_cast<KK_FLOAT>(30.0)*K0)*g3;
    const KK_FLOAT g0m = g0*(static_cast<KK_FLOAT>(42.0)*K0/K3+static_cast<KK_FLOAT>(6.0)*K1+K3);
    const KK_FLOAT g1m = g1*(static_cast<KK_FLOAT>(42.0)*K0/K4+static_cast<KK_FLOAT>(6.0)*K1+K4);
    const KK_FLOAT g2m = g2*(-static_cast<KK_FLOAT>(42.0)*K0/K5+static_cast<KK_FLOAT>(6.0)*K2+K5);
    const KK_FLOAT g3m = g3*(-static_cast<KK_FLOAT>(42.0)*K0/K6+static_cast<KK_FLOAT>(6.0)*K2+K6);
    const KK_FLOAT fR  = a12*sigma6/r/static_cast<KK_FLOAT>(37800.0);
    const KK_FLOAT evdwl_partial = fR*(h0-h1-h2+h3);
    const KK_FLOAT dUR = evdwl_partial/r + static_cast<KK_FLOAT>(5.0)*fR*(g0m+g1m-g2m-g3m);
    const KK_FLOAT dUA = -a12/static_cast<KK_FLOAT>(3.0)*r
                         *((static_cast<KK_FLOAT>(2.0)*K0*K7+static_cast<KK_FLOAT>(1.0))*K7
                           +(static_cast<KK_FLOAT>(2.0)*K0*K8-static_cast<KK_FLOAT>(1.0))*K8);
    return (dUR+dUA)/r;
  }
}

template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairColloidKokkos<DeviceType>::
compute_evdwl(const KK_FLOAT &rsq, const int &, const int &, const int &itype, const int &jtype) const {
  const int    iForm  = STACKPARAMS ? m_params[itype][jtype].form   : params(itype,jtype).form;
  const KK_FLOAT a12   = STACKPARAMS ? m_params[itype][jtype].a12    : params(itype,jtype).a12;
  const KK_FLOAT sigma3 = STACKPARAMS ? m_params[itype][jtype].sigma3 : params(itype,jtype).sigma3;
  const KK_FLOAT sigma6 = STACKPARAMS ? m_params[itype][jtype].sigma6 : params(itype,jtype).sigma6;
  const KK_FLOAT lj3   = STACKPARAMS ? m_params[itype][jtype].lj3   : params(itype,jtype).lj3;
  const KK_FLOAT lj4   = STACKPARAMS ? m_params[itype][jtype].lj4   : params(itype,jtype).lj4;
  const KK_FLOAT ppA1  = STACKPARAMS ? m_params[itype][jtype].a1    : params(itype,jtype).a1;
  const KK_FLOAT ppA2  = STACKPARAMS ? m_params[itype][jtype].a2    : params(itype,jtype).a2;
  const KK_FLOAT offset = STACKPARAMS ? m_params[itype][jtype].offset : params(itype,jtype).offset;

  if (iForm == SMALL_SMALL) {
    const KK_FLOAT r2inv = static_cast<KK_FLOAT>(1.0) / rsq;
    const KK_FLOAT r6inv = r2inv*r2inv*r2inv;
    return r6inv*(r6inv*lj3 - lj4) - offset;

  } else if (iForm == SMALL_LARGE) {
    const KK_FLOAT c2  = ppA2;
    const KK_FLOAT K1  = c2*c2;
    const KK_FLOAT K2  = rsq;
    const KK_FLOAT K4  = rsq*rsq;
    KK_FLOAT K3 = K1 - K2;
    K3 *= K3*K3;
    const KK_FLOAT K6  = K3*K3;
    const KK_FLOAT fR  = sigma3*a12*c2*K1/K3;
    return static_cast<KK_FLOAT>(2.0/9.0)*fR
           * (static_cast<KK_FLOAT>(1.0)
              - (K1*(K1*(K1/static_cast<KK_FLOAT>(3.0)+static_cast<KK_FLOAT>(3.0)*K2)
                     +static_cast<KK_FLOAT>(4.2)*K4)+K2*K4)
              * sigma6/K6) - offset;

  } else { // LARGE_LARGE
    const KK_FLOAT r  = sqrt(rsq);
    const KK_FLOAT c1 = ppA1;
    const KK_FLOAT c2 = ppA2;
    const KK_FLOAT K0 = c1*c2;
    const KK_FLOAT K1 = c1+c2;
    const KK_FLOAT K2 = c1-c2;
    const KK_FLOAT K3 = K1+r;
    const KK_FLOAT K4 = K1-r;
    const KK_FLOAT K5 = K2+r;
    const KK_FLOAT K6 = K2-r;
    const KK_FLOAT K7 = static_cast<KK_FLOAT>(1.0)/(K3*K4);
    const KK_FLOAT K8 = static_cast<KK_FLOAT>(1.0)/(K5*K6);
    const KK_FLOAT g0 = powint(K3,-7);
    const KK_FLOAT g1 = powint(K4,-7);
    const KK_FLOAT g2 = powint(K5,-7);
    const KK_FLOAT g3 = powint(K6,-7);
    const KK_FLOAT h0 = ((K3+static_cast<KK_FLOAT>(5.0)*K1)*K3+static_cast<KK_FLOAT>(30.0)*K0)*g0;
    const KK_FLOAT h1 = ((K4+static_cast<KK_FLOAT>(5.0)*K1)*K4+static_cast<KK_FLOAT>(30.0)*K0)*g1;
    const KK_FLOAT h2 = ((K5+static_cast<KK_FLOAT>(5.0)*K2)*K5-static_cast<KK_FLOAT>(30.0)*K0)*g2;
    const KK_FLOAT h3 = ((K6+static_cast<KK_FLOAT>(5.0)*K2)*K6-static_cast<KK_FLOAT>(30.0)*K0)*g3;
    const KK_FLOAT fR  = a12*sigma6/r/static_cast<KK_FLOAT>(37800.0);
    return fR*(h0-h1-h2+h3)
           + a12/static_cast<KK_FLOAT>(6.0)
             *(static_cast<KK_FLOAT>(2.0)*K0*(K7+K8)-log(K8/K7)) - offset;
  }
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

template<class DeviceType>
void PairColloidKokkos<DeviceType>::allocate()
{
  PairColloid::allocate();

  int n = atom->ntypes;
  memory->destroy(cutsq);
  memoryKK->create_kokkos(k_cutsq,cutsq,n+1,n+1,"pair:cutsq");
  d_cutsq = k_cutsq.template view<DeviceType>();
  k_params = Kokkos::DualView<params_colloid**,Kokkos::LayoutRight,DeviceType>("PairColloid::params",n+1,n+1);
  params = k_params.template view<DeviceType>();
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

template<class DeviceType>
void PairColloidKokkos<DeviceType>::init_style()
{
  PairColloid::init_style();

  // error if rRESPA with inner levels

  if (update->whichflag == 1 && utils::strmatch(update->integrate_style,"^respa")) {
    int respa = 0;
    if (((Respa *) update->integrate)->level_inner >= 0) respa = 1;
    if (((Respa *) update->integrate)->level_middle >= 0) respa = 2;
    if (respa) error->all(FLERR,"Cannot use Kokkos pair style with rRESPA inner/middle");
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
double PairColloidKokkos<DeviceType>::init_one(int i, int j)
{
  double cutone = PairColloid::init_one(i,j);

  k_params.view_host()(i,j).form   = form[i][j];
  k_params.view_host()(i,j).a12    = static_cast<KK_FLOAT>(a12[i][j]);
  k_params.view_host()(i,j).sigma3 = static_cast<KK_FLOAT>(sigma3[i][j]);
  k_params.view_host()(i,j).sigma6 = static_cast<KK_FLOAT>(sigma6[i][j]);
  k_params.view_host()(i,j).lj1    = static_cast<KK_FLOAT>(lj1[i][j]);
  k_params.view_host()(i,j).lj2    = static_cast<KK_FLOAT>(lj2[i][j]);
  k_params.view_host()(i,j).lj3    = static_cast<KK_FLOAT>(lj3[i][j]);
  k_params.view_host()(i,j).lj4    = static_cast<KK_FLOAT>(lj4[i][j]);
  k_params.view_host()(i,j).a1     = static_cast<KK_FLOAT>(a1[i][j]);
  k_params.view_host()(i,j).a2     = static_cast<KK_FLOAT>(a2[i][j]);
  k_params.view_host()(i,j).offset = static_cast<KK_FLOAT>(offset[i][j]);
  k_params.view_host()(i,j).cutsq  = static_cast<KK_FLOAT>(cutone*cutone);
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
template class PairColloidKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairColloidKokkos<LMPHostType>;
#endif
}
