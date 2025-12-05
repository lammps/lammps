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

  // no virial support yet. remove this after adding it.
  virial_global_flag = virial_peratom_flag = 0;

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
void FixAddForceKokkos<DeviceType>::post_force(int /*vflag*/)
{
  atomKK->sync(execution_space, X_MASK | F_MASK | IMAGE_MASK | MASK_MASK);

  x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  image = atomKK->k_image.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();

  int nlocal = atom->nlocal;

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
  double_4 foriginal_kk;
  force_flag = 0;
  prd = domain->prd;
  h = domain->h;
  triclinic = domain->triclinic;

  if (varflag == CONSTANT) {
    copymode = 1;
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagFixAddForceConstant>(0,nlocal),*this,foriginal_kk);
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
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagFixAddForceNonConstant>(0,nlocal),*this,foriginal_kk);
    copymode = 0;
  }

  atomKK->modified(execution_space, F_MASK);

  foriginal[0] = foriginal_kk.d0;
  foriginal[1] = foriginal_kk.d1;
  foriginal[2] = foriginal_kk.d2;
  foriginal[3] = foriginal_kk.d3;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixAddForceKokkos<DeviceType>::operator()(TagFixAddForceConstant, const int &i, double_4& foriginal_kk) const {
  if (mask[i] & groupbit) {
    if (region && !d_match[i]) return;

    Few<double,3> x_i;
    x_i[0] = x(i,0);
    x_i[1] = x(i,1);
    x_i[2] = x(i,2);
    auto unwrapKK = DomainKokkos::unmap(prd,h,triclinic,x_i,image(i));

    foriginal_kk.d0 -= xvalue * unwrapKK[0] + yvalue * unwrapKK[1] + zvalue * unwrapKK[2];
    foriginal_kk.d1 += f(i,0);
    foriginal_kk.d2 += f(i,1);
    foriginal_kk.d3 += f(i,2);
    if (xstyle) f(i,0) += xvalue;
    if (ystyle) f(i,1) += yvalue;
    if (zstyle) f(i,2) += zvalue;
  }
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixAddForceKokkos<DeviceType>::operator()(TagFixAddForceNonConstant, const int &i, double_4& foriginal_kk) const {
  if (mask[i] & groupbit) {
    if (region && !d_match[i]) return;

    Few<double,3> x_i;
    x_i[0] = x(i,0);
    x_i[1] = x(i,1);
    x_i[2] = x(i,2);
    auto unwrapKK = DomainKokkos::unmap(prd,h,triclinic,x_i,image(i));

    if (estyle == ATOM) {
      foriginal_kk.d0 += d_sforce(i,3);
    } else {
      if (xstyle == EQUAL) foriginal_kk.d0 -= xvalue * unwrapKK[0];
      if (ystyle == EQUAL) foriginal_kk.d0 -= yvalue * unwrapKK[1];
      if (zstyle == EQUAL) foriginal_kk.d0 -= zvalue * unwrapKK[2];
      if (xstyle == ATOM) foriginal_kk.d0 -= d_sforce(i,0) * unwrapKK[0];
      if (ystyle == ATOM) foriginal_kk.d0 -= d_sforce(i,1) * unwrapKK[1];
      if (zstyle == ATOM) foriginal_kk.d0 -= d_sforce(i,2) * unwrapKK[2];
    }
    foriginal_kk.d1 += f(i,0);
    foriginal_kk.d2 += f(i,1);
    foriginal_kk.d3 += f(i,2);
    if (xstyle == ATOM) f(i,0) += d_sforce(i,0);
    else if (xstyle) f(i,0) += xvalue;
    if (ystyle == ATOM) f(i,1) += d_sforce(i,1);
    else if (ystyle) f(i,1) += yvalue;
    if (zstyle == ATOM) f(i,2) += d_sforce(i,2);
    else if (zstyle) f(i,2) += zvalue;
  }
}

namespace LAMMPS_NS {
template class FixAddForceKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixAddForceKokkos<LMPHostType>;
#endif
}

