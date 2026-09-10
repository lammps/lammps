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

#include "pair_buck_mdf_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "error.h"
#include "force.h"
#include "kokkos.h"
#include "memory_kokkos.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "update.h"

#include <cmath>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairBuckMDFKokkos<DeviceType>::PairBuckMDFKokkos(LAMMPS *lmp) : PairBuckMDF(lmp)
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
PairBuckMDFKokkos<DeviceType>::~PairBuckMDFKokkos()
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
void PairBuckMDFKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
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
  special_lj[0] = static_cast<KK_FLOAT>(force->special_lj[0]);
  special_lj[1] = static_cast<KK_FLOAT>(force->special_lj[1]);
  special_lj[2] = static_cast<KK_FLOAT>(force->special_lj[2]);
  special_lj[3] = static_cast<KK_FLOAT>(force->special_lj[3]);

  // loop over neighbors of my atoms

  copymode = 1;

  EV_FLOAT ev = pair_compute<PairBuckMDFKokkos<DeviceType>,void >(this,(NeighListKokkos<DeviceType>*)list);

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
   MDF (Mei-Davenport-Fernando) tapering functions
   tt = taper function value (goes 1 -> 0 from cut_inner to cut)
   dt = -r * d(tt)/dr (switching derivative times r)
------------------------------------------------------------------------- */

template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairBuckMDFKokkos<DeviceType>::
compute_fpair(const KK_FLOAT &rsq, const int &, const int &, const int &itype, const int &jtype) const {
  const KK_FLOAT r = Kokkos::sqrt(rsq);
  const KK_FLOAT r2inv = static_cast<KK_FLOAT>(1.0)/rsq;
  const KK_FLOAT r6inv = r2inv*r2inv*r2inv;
  const KK_FLOAT rhoinv_val    = (STACKPARAMS?m_params[itype][jtype].rhoinv:params(itype,jtype).rhoinv);
  const KK_FLOAT buck1_val     = (STACKPARAMS?m_params[itype][jtype].buck1:params(itype,jtype).buck1);
  const KK_FLOAT buck2_val     = (STACKPARAMS?m_params[itype][jtype].buck2:params(itype,jtype).buck2);
  const KK_FLOAT a_val         = (STACKPARAMS?m_params[itype][jtype].a:params(itype,jtype).a);
  const KK_FLOAT c_val         = (STACKPARAMS?m_params[itype][jtype].c:params(itype,jtype).c);
  const KK_FLOAT cut_inner_val = (STACKPARAMS?m_params[itype][jtype].cut_inner:params(itype,jtype).cut_inner);
  const KK_FLOAT cut_inner_sq  = (STACKPARAMS?m_params[itype][jtype].cut_inner_sq:params(itype,jtype).cut_inner_sq);
  const KK_FLOAT cutsq_val     = (STACKPARAMS?m_params[itype][jtype].cutsq:params(itype,jtype).cutsq);
  const KK_FLOAT rexp = Kokkos::exp(-r*rhoinv_val);
  KK_FLOAT forcebuck = buck1_val*r*rexp - buck2_val*r6inv;

  if (rsq > cut_inner_sq) {
    const KK_FLOAT cut_val = Kokkos::sqrt(cutsq_val);
    const KK_FLOAT dp = cut_val - cut_inner_val;
    const KK_FLOAT d = (r - cut_inner_val)/dp;
    const KK_FLOAT dd = static_cast<KK_FLOAT>(1.0) - d;
    // tapering function
    const KK_FLOAT tt = (static_cast<KK_FLOAT>(1.0) + static_cast<KK_FLOAT>(3.0)*d
                         + static_cast<KK_FLOAT>(6.0)*d*d)*dd*dd*dd;
    // minus derivative of tapering function times r
    const KK_FLOAT dt = static_cast<KK_FLOAT>(30.0)*d*d*dd*dd*r/dp;
    const KK_FLOAT phibuck = a_val*rexp - c_val*r6inv;
    forcebuck = forcebuck*tt + phibuck*dt;
  }

  return forcebuck*r2inv;
}

template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairBuckMDFKokkos<DeviceType>::
compute_evdwl(const KK_FLOAT &rsq, const int &, const int &, const int &itype, const int &jtype) const {
  const KK_FLOAT r = Kokkos::sqrt(rsq);
  const KK_FLOAT r2inv = static_cast<KK_FLOAT>(1.0)/rsq;
  const KK_FLOAT r6inv = r2inv*r2inv*r2inv;
  const KK_FLOAT rhoinv_val    = (STACKPARAMS?m_params[itype][jtype].rhoinv:params(itype,jtype).rhoinv);
  const KK_FLOAT a_val         = (STACKPARAMS?m_params[itype][jtype].a:params(itype,jtype).a);
  const KK_FLOAT c_val         = (STACKPARAMS?m_params[itype][jtype].c:params(itype,jtype).c);
  const KK_FLOAT cut_inner_val = (STACKPARAMS?m_params[itype][jtype].cut_inner:params(itype,jtype).cut_inner);
  const KK_FLOAT cut_inner_sq  = (STACKPARAMS?m_params[itype][jtype].cut_inner_sq:params(itype,jtype).cut_inner_sq);
  const KK_FLOAT cutsq_val     = (STACKPARAMS?m_params[itype][jtype].cutsq:params(itype,jtype).cutsq);
  const KK_FLOAT rexp = Kokkos::exp(-r*rhoinv_val);
  KK_FLOAT evdwl = a_val*rexp - c_val*r6inv;

  if (rsq > cut_inner_sq) {
    const KK_FLOAT cut_val = Kokkos::sqrt(cutsq_val);
    const KK_FLOAT dp = cut_val - cut_inner_val;
    const KK_FLOAT d = (r - cut_inner_val)/dp;
    const KK_FLOAT dd = static_cast<KK_FLOAT>(1.0) - d;
    const KK_FLOAT tt = (static_cast<KK_FLOAT>(1.0) + static_cast<KK_FLOAT>(3.0)*d
                         + static_cast<KK_FLOAT>(6.0)*d*d)*dd*dd*dd;
    evdwl *= tt;
  }

  return evdwl;
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

template<class DeviceType>
void PairBuckMDFKokkos<DeviceType>::allocate()
{
  PairBuckMDF::allocate();

  int n = atom->ntypes;
  memory->destroy(cutsq);
  memoryKK->create_kokkos(k_cutsq,cutsq,n+1,n+1,"pair:cutsq");
  d_cutsq = k_cutsq.template view<DeviceType>();
  k_params = Kokkos::DualView<params_buck_mdf**,Kokkos::LayoutRight,DeviceType>("PairBuckMDF::params",n+1,n+1);
  params = k_params.template view<DeviceType>();
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

template<class DeviceType>
void PairBuckMDFKokkos<DeviceType>::init_style()
{
  PairBuckMDF::init_style();

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
double PairBuckMDFKokkos<DeviceType>::init_one(int i, int j)
{
  double cutone = PairBuckMDF::init_one(i,j);

  k_params.view_host()(i,j).a           = static_cast<KK_FLOAT>(a[i][j]);
  k_params.view_host()(i,j).rhoinv      = static_cast<KK_FLOAT>(rhoinv[i][j]);
  k_params.view_host()(i,j).c           = static_cast<KK_FLOAT>(c[i][j]);
  k_params.view_host()(i,j).buck1       = static_cast<KK_FLOAT>(buck1[i][j]);
  k_params.view_host()(i,j).buck2       = static_cast<KK_FLOAT>(buck2[i][j]);
  k_params.view_host()(i,j).cut_inner   = static_cast<KK_FLOAT>(cut_inner[i][j]);
  k_params.view_host()(i,j).cut_inner_sq= static_cast<KK_FLOAT>(cut_inner_sq[i][j]);
  k_params.view_host()(i,j).cutsq       = static_cast<KK_FLOAT>(cutone*cutone);
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
template class PairBuckMDFKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairBuckMDFKokkos<LMPHostType>;
#endif
}
