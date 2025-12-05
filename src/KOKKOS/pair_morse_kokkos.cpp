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
   Contributing authors: Stefan Paquay (Eindhoven University of Technology)
------------------------------------------------------------------------- */

#include "pair_morse_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "error.h"
#include "force.h"
#include "kokkos.h"
#include "math_const.h"
#include "memory_kokkos.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "respa.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairMorseKokkos<DeviceType>::PairMorseKokkos(LAMMPS *lmp) : PairMorse(lmp)
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
PairMorseKokkos<DeviceType>::~PairMorseKokkos()
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
void PairMorseKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
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

  EV_FLOAT ev = pair_compute<PairMorseKokkos<DeviceType>,void >(this,(NeighListKokkos<DeviceType>*)list);

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
}

template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairMorseKokkos<DeviceType>::
compute_fpair(const KK_FLOAT &rsq, const int &, const int &, const int &itype, const int &jtype) const {
  const KK_FLOAT rr = sqrt(rsq);
  const KK_FLOAT r0 = STACKPARAMS ? m_params[itype][jtype].r0 : params(itype,jtype).r0;
  const KK_FLOAT d0 = STACKPARAMS ? m_params[itype][jtype].d0 : params(itype,jtype).d0;
  const KK_FLOAT aa = STACKPARAMS ? m_params[itype][jtype].alpha : params(itype,jtype).alpha;
  const KK_FLOAT dr = rr - r0;

  // U  =  d0 * [ exp( -2*a*(x-r0)) - 2*exp(-a*(x-r0)) ]
  // f  = -2*a*d0*[ -exp( -2*a*(x-r0) ) + exp( -a*(x-r0) ) ] * grad(r)
  //    = +2*a*d0*[  exp( -2*a*(x-r0) ) - exp( -a*(x-r0) ) ] * grad(r)
  const KK_FLOAT dexp    = exp( -aa*dr );
  const KK_FLOAT forcelj = 2*aa*d0*dexp*(dexp-1.0);

  return forcelj / rr;
}

template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairMorseKokkos<DeviceType>::
compute_evdwl(const KK_FLOAT &rsq, const int &, const int &, const int &itype, const int &jtype) const {
  const KK_FLOAT rr = sqrt(rsq);
  const KK_FLOAT r0 = STACKPARAMS ? m_params[itype][jtype].r0 : params(itype,jtype).r0;
  const KK_FLOAT d0 = STACKPARAMS ? m_params[itype][jtype].d0 : params(itype,jtype).d0;
  const KK_FLOAT aa = STACKPARAMS ? m_params[itype][jtype].alpha : params(itype,jtype).alpha;
  const KK_FLOAT dr = rr - r0;

  // U  =  d0 * [ exp( -2*a*(x-r0)) - 2*exp(-a*(x-r0)) ]
  // f  = -2*a*d0*[ -exp( -2*a*(x-r0) ) + exp( -a*(x-r0) ) ] * grad(r)
  //    = +2*a*d0*[  exp( -2*a*(x-r0) ) - exp( -a*(x-r0) ) ] * grad(r)
  const KK_FLOAT dexp    = exp( -aa*dr );

  return d0 * dexp * ( dexp - 2.0 );
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

template<class DeviceType>
void PairMorseKokkos<DeviceType>::allocate()
{
  PairMorse::allocate();

  int n = atom->ntypes;
  memory->destroy(cutsq);
  memoryKK->create_kokkos(k_cutsq,cutsq,n+1,n+1,"pair:cutsq");
  d_cutsq = k_cutsq.template view<DeviceType>();
  k_params = Kokkos::DualView<params_morse**,Kokkos::LayoutRight,DeviceType>("PairMorse::params",n+1,n+1);
  params = k_params.template view<DeviceType>();
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

template<class DeviceType>
void PairMorseKokkos<DeviceType>::init_style()
{
  PairMorse::init_style();

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
// Rewrite this.
template<class DeviceType>
double PairMorseKokkos<DeviceType>::init_one(int i, int j)
{
  double cutone = PairMorse::init_one(i,j);

  k_params.view_host()(i,j).d0     = d0[i][j];
  k_params.view_host()(i,j).alpha  = alpha[i][j];
  k_params.view_host()(i,j).r0     = r0[i][j];
  k_params.view_host()(i,j).offset = offset[i][j];
  k_params.view_host()(i,j).cutsq  = cutone*cutone;
  k_params.view_host()(j,i)        = k_params.view_host()(i,j);

  if (i<MAX_TYPES_STACKPARAMS+1 && j<MAX_TYPES_STACKPARAMS+1) {
    m_params[i][j] = m_params[j][i] = k_params.view_host()(i,j);
    m_cutsq[j][i] = m_cutsq[i][j] = cutone*cutone;
  }

  k_cutsq.view_host()(i,j) = k_cutsq.view_host()(j,i) = cutone*cutone;
  k_cutsq.modify_host();
  k_params.modify_host();

  return cutone;
}



namespace LAMMPS_NS {
template class PairMorseKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairMorseKokkos<LMPHostType>;
#endif
}

