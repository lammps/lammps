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
   Contributing author: Trung Nguyen (U Chicago)
------------------------------------------------------------------------- */

#include "pair_yukawa_colloid_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "error.h"
#include "force.h"
#include "kokkos.h"
#include "memory_kokkos.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "respa.h"
#include "update.h"

#include <cmath>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairYukawaColloidKokkos<DeviceType>::PairYukawaColloidKokkos(LAMMPS *lmp) : PairYukawaColloid(lmp)
{
  respa_enable = 0;

  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | F_MASK | TYPE_MASK | ENERGY_MASK | VIRIAL_MASK | RADIUS_MASK;
  datamask_modify = F_MASK | ENERGY_MASK | VIRIAL_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairYukawaColloidKokkos<DeviceType>::~PairYukawaColloidKokkos()
{
  if (copymode) return;

  if (allocated) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->destroy_kokkos(k_vatom,vatom);
    memoryKK->destroy_kokkos(k_cutsq,cutsq);
  }
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

template<class DeviceType>
void PairYukawaColloidKokkos<DeviceType>::allocate()
{
  PairYukawaColloid::allocate();

  int n = atom->ntypes;
  memory->destroy(cutsq);
  memoryKK->create_kokkos(k_cutsq,cutsq,n+1,n+1,"pair:cutsq");
  d_cutsq = k_cutsq.template view<DeviceType>();
  k_params = Kokkos::DualView<params_yukawa**,
                              Kokkos::LayoutRight,DeviceType>(
                              "PairYukawaColloid::params",n+1,n+1);

  params = k_params.template view<DeviceType>();
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

template<class DeviceType>
void PairYukawaColloidKokkos<DeviceType>::init_style()
{
  PairYukawaColloid::init_style();

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
  request->set_kokkos_host(std::is_same<DeviceType,LMPHostType>::value &&
                           !std::is_same<DeviceType,LMPDeviceType>::value);
  request->set_kokkos_device(std::is_same<DeviceType,LMPDeviceType>::value);
  if (neighflag == FULL) request->enable_full();
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */
// Rewrite this.
template<class DeviceType>
double PairYukawaColloidKokkos<DeviceType>::init_one(int i, int j)
{
  double cutone = PairYukawaColloid::init_one(i,j);

  k_params.h_view(i,j).a      = a[i][j];
  k_params.h_view(i,j).offset = offset[i][j];
  k_params.h_view(i,j).cutsq  = cutone*cutone;
  k_params.h_view(j,i)        = k_params.h_view(i,j);

  if (i<MAX_TYPES_STACKPARAMS+1 && j<MAX_TYPES_STACKPARAMS+1) {
    m_params[i][j] = m_params[j][i] = k_params.h_view(i,j);
    m_cutsq[j][i] = m_cutsq[i][j] = cutone*cutone;
  }

  k_cutsq.h_view(i,j) = k_cutsq.h_view(j,i) = cutone*cutone;
  k_cutsq.template modify<LMPHostType>();
  k_params.template modify<LMPHostType>();

  return cutone;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairYukawaColloidKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
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
  radius = atomKK->k_radius.view<DeviceType>();
  nlocal = atom->nlocal;
  nall = atom->nlocal + atom->nghost;
  newton_pair = force->newton_pair;
  special_lj[0] = force->special_lj[0];
  special_lj[1] = force->special_lj[1];
  special_lj[2] = force->special_lj[2];
  special_lj[3] = force->special_lj[3];

  // loop over neighbors of my atoms

  EV_FLOAT ev = pair_compute<PairYukawaColloidKokkos<DeviceType>,void >(
    this,(NeighListKokkos<DeviceType>*)list);

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
    k_eatom.template sync<LMPHostType>();
  }

  if (vflag_atom) {
    k_vatom.template modify<DeviceType>();
    k_vatom.template sync<LMPHostType>();
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
KOKKOS_INLINE_FUNCTION
F_FLOAT PairYukawaColloidKokkos<DeviceType>::
compute_fpair(const F_FLOAT& rsq, const int& i, const int&j,
              const int& itype, const int& jtype) const {
  (void) i;
  (void) j;
  const F_FLOAT radi   = radius[i];
  const F_FLOAT radj   = radius[j];
  const F_FLOAT rr     = sqrt(rsq);
  // Fetch the params either off the stack or from some mapped memory?
  const F_FLOAT aa     = STACKPARAMS ? m_params[itype][jtype].a
                                     : params(itype,jtype).a;

  // U   = a * exp(-kappa*(r-(radi+radj))) / kappa
  // f   = -dU/dr = a * exp(-kappa*r)
  // f/r = a * exp(-kappa*r) / r
  const F_FLOAT rinv = 1.0 / rr;
  const F_FLOAT screening = exp(-kappa*(rr-(radi+radj)));
  const F_FLOAT forceyukawa = aa * screening;
  const F_FLOAT fpair = forceyukawa * rinv;

  return fpair;
}

template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
KOKKOS_INLINE_FUNCTION
F_FLOAT PairYukawaColloidKokkos<DeviceType>::
compute_evdwl(const F_FLOAT& rsq, const int& i, const int&j,
              const int& itype, const int& jtype) const {
  (void) i;
  (void) j;
  const F_FLOAT radi   = radius[i];
  const F_FLOAT radj   = radius[j];
  const F_FLOAT rr     = sqrt(rsq);
  const F_FLOAT aa     = STACKPARAMS ? m_params[itype][jtype].a
                                     : params(itype,jtype).a;
  const F_FLOAT offset = STACKPARAMS ? m_params[itype][jtype].offset
                                     : params(itype,jtype).offset;

  // U   = a * exp(-kappa*(r-(radi+radj))) / kappa
  const F_FLOAT screening = exp(-kappa*(rr-(radi+radj)));

  return aa / kappa * screening - offset;
}


namespace LAMMPS_NS {
template class PairYukawaColloidKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairYukawaColloidKokkos<LMPHostType>;
#endif
}
