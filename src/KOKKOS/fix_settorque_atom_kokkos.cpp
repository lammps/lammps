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

#include "fix_settorque_atom_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "error.h"
#include "input.h"
#include "kokkos_base.h"
#include "memory_kokkos.h"
#include "modify.h"
#include "region.h"
#include "update.h"
#include "variable.h"

using namespace LAMMPS_NS;

enum { NONE, CONSTANT, EQUAL, ATOM };

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixSetTorqueAtomKokkos<DeviceType>::FixSetTorqueAtomKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixSetTorqueAtom(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = TORQUE_MASK | MASK_MASK;
  datamask_modify = TORQUE_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixSetTorqueAtomKokkos<DeviceType>::~FixSetTorqueAtomKokkos()
{
  if (copymode) return;
  memoryKK->destroy_kokkos(k_storque,storque);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixSetTorqueAtomKokkos<DeviceType>::init()
{
  FixSetTorqueAtom::init();

  if (utils::strmatch(update->integrate_style,"^respa"))
    error->all(FLERR,"Cannot (yet) use respa with fix settorque/atom/kk");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixSetTorqueAtomKokkos<DeviceType>::post_force(int /*vflag*/)
{
  atomKK->sync(execution_space,datamask_read);

  d_torque = atomKK->k_torque.view<DeviceType>();
  d_mask = atomKK->k_mask.view<DeviceType>();

  const int nlocal = atom->nlocal;

  // update region if necessary

  if (region) {
    auto *regionKKBase = dynamic_cast<KokkosBase *>(region);
    if (!regionKKBase)
      error->all(FLERR,"Cannot (yet) use {}-style region with fix settorque/atom/kk",
                 region->style);
    region->prematch();
    DAT::tdual_int_1d k_match = DAT::tdual_int_1d("settorque/atom:k_match",nlocal);
    regionKKBase->match_all_kokkos(groupbit,k_match);
    k_match.template sync<DeviceType>();
    d_match = k_match.template view<DeviceType>();
  }

  // reallocate storque array if necessary

  if (varflag == ATOM && atom->nmax > maxatom) {
    maxatom = atom->nmax;
    memoryKK->destroy_kokkos(k_storque,storque);
    memoryKK->create_kokkos(k_storque,storque,maxatom,3,"settorque/atom:storque");
    d_storque = k_storque.view<DeviceType>();
  }

  double_3 toriginal_kk;
  torque_flag = 0;

  if (varflag == CONSTANT) {
    copymode = 1;
    Kokkos::parallel_reduce(
      Kokkos::RangePolicy<DeviceType,TagFixSetTorqueAtomConstant>(0,nlocal),*this,toriginal_kk);
    copymode = 0;

  // variable torques, wrap with clear/add

  } else {

    modify->clearstep_compute();

    if (xstyle == EQUAL) xvalue = input->variable->compute_equal(xvar);
    else if (xstyle == ATOM)
      input->variable->compute_atom(xvar,igroup,&storque[0][0],3,0);
    if (ystyle == EQUAL) yvalue = input->variable->compute_equal(yvar);
    else if (ystyle == ATOM)
      input->variable->compute_atom(yvar,igroup,&storque[0][1],3,0);
    if (zstyle == EQUAL) zvalue = input->variable->compute_equal(zvar);
    else if (zstyle == ATOM)
      input->variable->compute_atom(zvar,igroup,&storque[0][2],3,0);

    modify->addstep_compute(update->ntimestep + 1);

    // atom-style variables are evaluated on the host, so the result has to be
    // copied to the device for the kernel below

    if (varflag == ATOM) {
      k_storque.modify_host();
      k_storque.sync<DeviceType>();
    }

    copymode = 1;
    Kokkos::parallel_reduce(
      Kokkos::RangePolicy<DeviceType,TagFixSetTorqueAtomNonConstant>(0,nlocal),*this,toriginal_kk);
    copymode = 0;
  }

  atomKK->modified(execution_space,datamask_modify);

  toriginal[0] = toriginal_kk.d0;
  toriginal[1] = toriginal_kk.d1;
  toriginal[2] = toriginal_kk.d2;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixSetTorqueAtomKokkos<DeviceType>::operator()(TagFixSetTorqueAtomConstant, const int &i,
                                                    double_3 &toriginal_kk) const
{
  const KK_FLOAT xvalue_kk = static_cast<KK_FLOAT>(xvalue);
  const KK_FLOAT yvalue_kk = static_cast<KK_FLOAT>(yvalue);
  const KK_FLOAT zvalue_kk = static_cast<KK_FLOAT>(zvalue);
  if (d_mask[i] & groupbit) {
    if (region && !d_match[i]) return;
    toriginal_kk.d0 += static_cast<double>(d_torque(i,0));
    toriginal_kk.d1 += static_cast<double>(d_torque(i,1));
    toriginal_kk.d2 += static_cast<double>(d_torque(i,2));
    if (xstyle) d_torque(i,0) = static_cast<KK_ACC_FLOAT>(xvalue_kk);
    if (ystyle) d_torque(i,1) = static_cast<KK_ACC_FLOAT>(yvalue_kk);
    if (zstyle) d_torque(i,2) = static_cast<KK_ACC_FLOAT>(zvalue_kk);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixSetTorqueAtomKokkos<DeviceType>::operator()(TagFixSetTorqueAtomNonConstant, const int &i,
                                                    double_3 &toriginal_kk) const
{
  const KK_FLOAT xvalue_kk = static_cast<KK_FLOAT>(xvalue);
  const KK_FLOAT yvalue_kk = static_cast<KK_FLOAT>(yvalue);
  const KK_FLOAT zvalue_kk = static_cast<KK_FLOAT>(zvalue);
  if (d_mask[i] & groupbit) {
    if (region && !d_match[i]) return;
    toriginal_kk.d0 += static_cast<double>(d_torque(i,0));
    toriginal_kk.d1 += static_cast<double>(d_torque(i,1));
    toriginal_kk.d2 += static_cast<double>(d_torque(i,2));
    if (xstyle == ATOM) d_torque(i,0) = static_cast<KK_ACC_FLOAT>(d_storque(i,0));
    else if (xstyle) d_torque(i,0) = static_cast<KK_ACC_FLOAT>(xvalue_kk);
    if (ystyle == ATOM) d_torque(i,1) = static_cast<KK_ACC_FLOAT>(d_storque(i,1));
    else if (ystyle) d_torque(i,1) = static_cast<KK_ACC_FLOAT>(yvalue_kk);
    if (zstyle == ATOM) d_torque(i,2) = static_cast<KK_ACC_FLOAT>(d_storque(i,2));
    else if (zstyle) d_torque(i,2) = static_cast<KK_ACC_FLOAT>(zvalue_kk);
  }
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixSetTorqueAtomKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixSetTorqueAtomKokkos<LMPHostType>;
#endif
}
