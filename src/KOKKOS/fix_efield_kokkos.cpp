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

#include "fix_efield_kokkos.h"

#include "atom_kokkos.h"
#include "update.h"
#include "modify.h"
#include "domain_kokkos.h"
#include "region.h"
#include "input.h"
#include "variable.h"
#include "memory_kokkos.h"
#include "error.h"
#include "atom_masks.h"
#include "kokkos_base.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{NONE,CONSTANT,EQUAL,ATOM};

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixEfieldKokkos<DeviceType>::FixEfieldKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixEfield(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;

  memory->destroy(efield);
  memoryKK->create_kokkos(k_efield,efield,maxatom,4,"efield:efield");
  d_efield = k_efield.view<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixEfieldKokkos<DeviceType>::~FixEfieldKokkos()
{
  if (copymode) return;

  memoryKK->destroy_kokkos(k_efield,efield);
  efield = nullptr;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixEfieldKokkos<DeviceType>::init()
{
  FixEfield::init();

  if (utils::strmatch(update->integrate_style,"^respa"))
    error->all(FLERR,"Cannot (yet) use respa with Kokkos");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixEfieldKokkos<DeviceType>::post_force(int /*vflag*/)
{
  atomKK->sync(execution_space, X_MASK | F_MASK | Q_MASK | IMAGE_MASK | MASK_MASK);

  x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  q = atomKK->k_q.view<DeviceType>();
  image = atomKK->k_image.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();

  int nlocal = atom->nlocal;

  // update region if necessary

  if (region) {
    if (!utils::strmatch(region->style, "^block"))
      error->all(FLERR,"Cannot (yet) use {}-style region with fix efield/kk",region->style);
    region->prematch();
    DAT::tdual_int_1d k_match = DAT::tdual_int_1d("efield:k_match",nlocal);
    KokkosBase* regionKKBase = dynamic_cast<KokkosBase*>(region);
    regionKKBase->match_all_kokkos(groupbit,k_match);
    k_match.template sync<DeviceType>();
    d_match = k_match.template view<DeviceType>();
  }

  // reallocate sforce array if necessary

  if (varflag == ATOM && atom->nmax > maxatom) {
    maxatom = atom->nmax;
    memoryKK->destroy_kokkos(k_efield,efield);
    memoryKK->create_kokkos(k_efield,efield,maxatom,4,"efield:efield");
    d_efield = k_efield.view<DeviceType>();
  }

  fsum[0] = fsum[1] = fsum[2] = fsum[3] = 0.0;
  double_4 fsum_kk;
  force_flag = 0;

  if (varflag == CONSTANT) {
    copymode = 1;

    // It would be more concise to use the operators below, but there is still an issue with unwrap (TODO below)
    //Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagFixEfieldConstant>(0,nlocal),*this,fsum_kk);

    {
    // local variables for lambda capture
    auto prd = Few<double,3>(domain->prd);
    auto h = Few<double,6>(domain->h);
    auto triclinic = domain->triclinic;
    auto l_ex = ex;
    auto l_ey = ey;
    auto l_ez = ez;

    auto l_x = x;
    auto l_q = q;
    auto l_f = f;
    auto l_mask = mask;
    auto l_image = image;
    auto l_groupbit = groupbit;

    Kokkos::parallel_reduce(nlocal, LAMMPS_LAMBDA(const int& i, double_4& fsum_kk) {
      if (l_mask[i] & l_groupbit) {
        Few<double,3> x_i;
        x_i[0] = l_x(i,0);
        x_i[1] = l_x(i,1);
        x_i[2] = l_x(i,2);
        auto unwrap = DomainKokkos::unmap(prd,h,triclinic,x_i,l_image(i));
        auto qtmp = l_q(i);
        auto fx = qtmp * l_ex;
        auto fy = qtmp * l_ey;
        auto fz = qtmp * l_ez;
        l_f(i,0) += fx;
        l_f(i,1) += fy;
        l_f(i,2) += fz;
        fsum_kk.d0 -= fx * unwrap[0] + fy * unwrap[1] + fz * unwrap[2];
        fsum_kk.d1 += fx;
        fsum_kk.d2 += fy;
        fsum_kk.d3 += fz;
      }
    },fsum_kk);
    }

    copymode = 0;

  // variable force, wrap with clear/add

  } else {

    atomKK->sync(Host,ALL_MASK); // this can be removed when variable class is ported to Kokkos

    modify->clearstep_compute();

    if (xstyle == EQUAL) ex = input->variable->compute_equal(xvar);
    else if (xstyle == ATOM)
      input->variable->compute_atom(xvar,igroup,&efield[0][0],4,0);
    if (ystyle == EQUAL) ey = input->variable->compute_equal(yvar);
    else if (ystyle == ATOM)
      input->variable->compute_atom(yvar,igroup,&efield[0][1],4,0);
    if (zstyle == EQUAL) ez = input->variable->compute_equal(zvar);
    else if (zstyle == ATOM)
      input->variable->compute_atom(zvar,igroup,&efield[0][2],4,0);

    modify->addstep_compute(update->ntimestep + 1);

    if (varflag == ATOM) {  // this can be removed when variable class is ported to Kokkos
      k_efield.modify<LMPHostType>();
      k_efield.sync<DeviceType>();
    }

    copymode = 1;
    // It would be more concise to use the operators below, but there is still an issue with unwrap (TODO below)
    //Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagFixEfieldNonConstant>(0,nlocal),*this,fsum_kk);
    {
    // local variables for lambda capture
    auto prd = Few<double,3>(domain->prd);
    auto h = Few<double,6>(domain->h);
    auto triclinic = domain->triclinic;
    auto l_ex = ex;
    auto l_ey = ey;
    auto l_ez = ez;
    auto l_d_efield = d_efield;

    auto l_x = x;
    auto l_q = q;
    auto l_f = f;
    auto l_mask = mask;
    auto l_image = image;
    auto l_groupbit = groupbit;
    auto l_xstyle = xstyle;
    auto l_ystyle = ystyle;
    auto l_zstyle = zstyle;

    Kokkos::parallel_reduce(nlocal, LAMMPS_LAMBDA(const int& i, double_4& fsum_kk) {
      if (l_mask[i] & l_groupbit) {
        Few<double,3> x_i;
        x_i[0] = l_x(i,0);
        x_i[1] = l_x(i,1);
        x_i[2] = l_x(i,2);
        auto unwrap = DomainKokkos::unmap(prd,h,triclinic,x_i,l_image(i));
        auto qtmp = l_q(i);
        auto fx = qtmp * l_ex;
        auto fy = qtmp * l_ey;
        auto fz = qtmp * l_ez;
        if (l_xstyle == ATOM) l_f(i,0) += qtmp * l_d_efield(i,0);
        else if (l_xstyle) l_f(i,0) += fx;
        if (l_ystyle == ATOM) l_f(i,1) += qtmp * l_d_efield(i,1);
        else if (l_ystyle) l_f(i,1) += fy;
        if (l_zstyle == ATOM) l_f(i,2) += qtmp * l_d_efield(i,2);
        else if (l_zstyle) l_f(i,2) += fz;
        fsum_kk.d0 -= fx * unwrap[0] + fy * unwrap[1] + fz * unwrap[2];
        fsum_kk.d1 += fx;
        fsum_kk.d2 += fy;
        fsum_kk.d3 += fz;
      }
    },fsum_kk);
    }

    copymode = 0;
  }

  atomKK->modified(execution_space, F_MASK);

  fsum[0] = fsum_kk.d0;
  fsum[1] = fsum_kk.d1;
  fsum[2] = fsum_kk.d2;
  fsum[3] = fsum_kk.d3;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixEfieldKokkos<DeviceType>::operator()(TagFixEfieldConstant, const int &i, double_4& fsum_kk) const {
  if (mask[i] & groupbit) {
    if (region && !d_match[i]) return;

    auto prd = Few<double,3>(domain->prd);
    auto h = Few<double,6>(domain->h);
    auto triclinic = domain->triclinic;
    Few<double,3> x_i;
    x_i[0] = x(i,0);
    x_i[1] = x(i,1);
    x_i[2] = x(i,2);
    auto unwrap = DomainKokkos::unmap(prd,h,triclinic,x_i,image(i));
    const F_FLOAT qtmp = q(i);
    const F_FLOAT fx = qtmp * ex;
    const F_FLOAT fy = qtmp * ey;
    const F_FLOAT fz = qtmp * ez;
    f(i,0) += fx;
    f(i,1) += fy;
    f(i,2) += fz;
    // TODO: access to unwrap below crashes
    fsum_kk.d0 -= fx * unwrap[0] + fy * unwrap[1] + fz * unwrap[2];
    fsum_kk.d1 += fx;
    fsum_kk.d2 += fy;
    fsum_kk.d3 += fz;
  }
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixEfieldKokkos<DeviceType>::operator()(TagFixEfieldNonConstant, const int &i, double_4& fsum_kk) const {
  auto prd = Few<double,3>(domain->prd);
  auto h = Few<double,6>(domain->h);
  auto triclinic = domain->triclinic;
  if (mask[i] & groupbit) {
    if (region && !d_match[i]) return;
    Few<double,3> x_i;
    x_i[0] = x(i,0);
    x_i[1] = x(i,1);
    x_i[2] = x(i,2);
    auto unwrap = DomainKokkos::unmap(prd,h,triclinic,x_i,image(i));
    const F_FLOAT qtmp = q[i];
    const F_FLOAT fx = qtmp * ex;
    const F_FLOAT fy = qtmp * ey;
    const F_FLOAT fz = qtmp * ez;
    if (xstyle == ATOM) f(i,0) += d_efield(i,0);
    else if (xstyle) f(i,0) += fx;
    if (ystyle == ATOM) f(i,1) += d_efield(i,1);
    else if (ystyle) f(i,1) += fy;
    if (zstyle == ATOM) f(i,2) += d_efield(i,2);
    else if (zstyle) f(i,2) += fz;
    // TODO: access to unwrap below crashes
    fsum_kk.d0 -= fx * unwrap[0] + fy * unwrap[1] + fz * unwrap[2];
    fsum_kk.d1 += fx;
    fsum_kk.d2 += fy;
    fsum_kk.d3 += fz;
  }
}

namespace LAMMPS_NS {
template class FixEfieldKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixEfieldKokkos<LMPHostType>;
#endif
}

