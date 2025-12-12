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

#include "fix_addforce_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "domain_kokkos.h"
#include "error.h"
#include "input.h"
#include "kokkos_base.h"
#include "memory_kokkos.h"
#include "modify.h"
#include "region.h"
#include "update.h"
#include "variable.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixAddForceKokkos<DeviceType>::FixAddForceKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixAddForce(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;

  memory->destroy(sforce);
  memoryKK->create_kokkos(k_sforce,sforce,maxatom,4,"addforce:sforce");
  d_sforce = k_sforce.view<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixAddForceKokkos<DeviceType>::~FixAddForceKokkos()
{
  if (copymode) return;

  memoryKK->destroy_kokkos(k_sforce,sforce);
  sforce = nullptr;
  memoryKK->destroy_kokkos(k_vatom,vatom);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixAddForceKokkos<DeviceType>::init()
{
  FixAddForce::init();

  if (utils::strmatch(update->integrate_style,"^respa"))
    error->all(FLERR,"Cannot (yet) use respa with Kokkos");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixAddForceKokkos<DeviceType>::post_force(int vflag)
{
  atomKK->sync(execution_space, X_MASK | F_MASK | IMAGE_MASK | MASK_MASK);

  x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  image = atomKK->k_image.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();
  int nlocal = atom->nlocal;

  // virial setup

  v_init(vflag);

  // reallocate per-atom arrays if necessary

  if (vflag_atom) {
    memoryKK->destroy_kokkos(k_vatom,vatom);
    memoryKK->create_kokkos(k_vatom,vatom,maxvatom,"efield:vatom");
    d_vatom = k_vatom.template view<DeviceType>();
  }

  // update region if necessary

  if (region) {
    if (!(utils::strmatch(region->style, "^block") || utils::strmatch(region->style, "^sphere")))
      error->all(FLERR,"Cannot (yet) use {}-style region with fix addforce/kk",region->style);
    region->prematch();
    DAT::tdual_int_1d k_match = DAT::tdual_int_1d("addforce:k_match",nlocal);
    KokkosBase* regionKKBase = dynamic_cast<KokkosBase*>(region);
    regionKKBase->match_all_kokkos(groupbit,k_match);
    k_match.template sync<DeviceType>();
    d_match = k_match.template view<DeviceType>();
  }

  // reallocate sforce array if necessary

  if (varflag == ATOM && atom->nmax > maxatom) {
    maxatom = atom->nmax;
    memoryKK->destroy_kokkos(k_sforce,sforce);
    memoryKK->create_kokkos(k_sforce,sforce,maxatom,4,"addforce:sforce");
    d_sforce = k_sforce.view<DeviceType>();
  }

  foriginal[0] = foriginal[1] = foriginal[2] = foriginal[3] = 0.0;
  force_flag = 0;
  double result[10] = {0.0};
  prd = domain->prd;
  h = domain->h;
  triclinic = domain->triclinic;

  if (varflag == CONSTANT) {
    copymode = 1;
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagFixAddForceConstant>(0,nlocal),*this,result);
    copymode = 0;

  // variable force, wrap with clear/add

  } else {

    atomKK->sync(Host,ALL_MASK); // this can be removed when variable class is ported to Kokkos

    modify->clearstep_compute();

    if (xstyle == EQUAL) xvalue = input->variable->compute_equal(xvar);
    else if (xstyle == ATOM)
      input->variable->compute_atom(xvar,igroup,&sforce[0][0],4,0);
    if (ystyle == EQUAL) yvalue = input->variable->compute_equal(yvar);
    else if (ystyle == ATOM)
      input->variable->compute_atom(yvar,igroup,&sforce[0][1],4,0);
    if (zstyle == EQUAL) zvalue = input->variable->compute_equal(zvar);
    else if (zstyle == ATOM)
      input->variable->compute_atom(zvar,igroup,&sforce[0][2],4,0);
    if (estyle == ATOM) input->variable->compute_atom(evar,igroup,&sforce[0][3],4,0);

    modify->addstep_compute(update->ntimestep + 1);

    if (varflag == ATOM) {  // this can be removed when variable class is ported to Kokkos
      k_sforce.modify_host();
      k_sforce.sync<DeviceType>();
    }

    copymode = 1;
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagFixAddForceNonConstant>(0,nlocal),*this,result);
    copymode = 0;
  }

  atomKK->modified(execution_space, F_MASK);

  foriginal[0] = result[0];
  foriginal[1] = result[1];
  foriginal[2] = result[2];
  foriginal[3] = result[3];

  if (vflag_global) {
    virial[0] += result[4];
    virial[1] += result[5];
    virial[2] += result[6];
    virial[3] += result[7];
    virial[4] += result[8];
    virial[5] += result[9];
  }

  if (vflag_atom) {
    k_vatom.template modify<DeviceType>();
    k_vatom.sync_host();
  }
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixAddForceKokkos<DeviceType>::operator()(TagFixAddForceConstant, const int &i, value_type result) const {
  if (mask[i] & groupbit) {
    if (region && !d_match[i]) return;

    Few<double,3> x_i;
    x_i[0] = x(i,0);
    x_i[1] = x(i,1);
    x_i[2] = x(i,2);
    auto unwrapKK = DomainKokkos::unmap(prd,h,triclinic,x_i,image(i));

    result[0] -= xvalue * unwrapKK[0] + yvalue * unwrapKK[1] + zvalue * unwrapKK[2];
    result[1] += f(i,0);
    result[2] += f(i,1);
    result[3] += f(i,2);
    if (xstyle) f(i,0) += xvalue;
    if (ystyle) f(i,1) += yvalue;
    if (zstyle) f(i,2) += zvalue;
  }
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixAddForceKokkos<DeviceType>::operator()(TagFixAddForceNonConstant, const int &i, value_type result) const {
  if (mask[i] & groupbit) {
    if (region && !d_match[i]) return;

    Few<double,3> x_i;
    x_i[0] = x(i,0);
    x_i[1] = x(i,1);
    x_i[2] = x(i,2);
    auto unwrapKK = DomainKokkos::unmap(prd,h,triclinic,x_i,image(i));

    if (estyle == ATOM) {
      result[0] += d_sforce(i,3);
    } else {
      if (xstyle == EQUAL) result[0] -= xvalue * unwrapKK[0];
      if (ystyle == EQUAL) result[0] -= yvalue * unwrapKK[1];
      if (zstyle == EQUAL) result[0] -= zvalue * unwrapKK[2];
      if (xstyle == ATOM) result[0] -= d_sforce(i,0) * unwrapKK[0];
      if (ystyle == ATOM) result[0] -= d_sforce(i,1) * unwrapKK[1];
      if (zstyle == ATOM) result[0] -= d_sforce(i,2) * unwrapKK[2];
    }
    result[1] += f(i,0);
    result[2] += f(i,1);
    result[3] += f(i,2);
    if (xstyle == ATOM) f(i,0) += d_sforce(i,0);
    else if (xstyle) f(i,0) += xvalue;
    if (ystyle == ATOM) f(i,1) += d_sforce(i,1);
    else if (ystyle) f(i,1) += yvalue;
    if (zstyle == ATOM) f(i,2) += d_sforce(i,2);
    else if (zstyle) f(i,2) += zvalue;
  }
}

/* ----------------------------------------------------------------------
   tally virial into global and per-atom accumulators
   i = local index of atom
   v = total virial for the interaction
   increment global virial by v
   increment per-atom virial by v
   this method can be used when fix computes forces in post_force()
   and the force depends on a distance to some external object
     e.g. fix wall/lj93: compute virial only on owned atoms
------------------------------------------------------------------------- */

template <class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixAddForceKokkos<DeviceType>::v_tally(value_type result, int i, KK_FLOAT *v) const
{
  if (vflag_global) {
    result[4] += static_cast<KK_ACC_FLOAT>(v[0]);
    result[5] += static_cast<KK_ACC_FLOAT>(v[1]);
    result[6] += static_cast<KK_ACC_FLOAT>(v[2]);
    result[7] += static_cast<KK_ACC_FLOAT>(v[3]);
    result[8] += static_cast<KK_ACC_FLOAT>(v[4]);
    result[9] += static_cast<KK_ACC_FLOAT>(v[5]);
  }

  if (vflag_atom) {
    Kokkos::atomic_add(&(d_vatom(i,0)),static_cast<KK_ACC_FLOAT>(v[0]));
    Kokkos::atomic_add(&(d_vatom(i,1)),static_cast<KK_ACC_FLOAT>(v[1]));
    Kokkos::atomic_add(&(d_vatom(i,2)),static_cast<KK_ACC_FLOAT>(v[2]));
    Kokkos::atomic_add(&(d_vatom(i,3)),static_cast<KK_ACC_FLOAT>(v[3]));
    Kokkos::atomic_add(&(d_vatom(i,4)),static_cast<KK_ACC_FLOAT>(v[4]));
    Kokkos::atomic_add(&(d_vatom(i,5)),static_cast<KK_ACC_FLOAT>(v[5]));
  }
}

namespace LAMMPS_NS {
template class FixAddForceKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixAddForceKokkos<LMPHostType>;
#endif
}
