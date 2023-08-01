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

#include "fix_efield_kokkos.h"

#include "atom_kokkos.h"
#include "update.h"
#include "modify.h"
#include "domain.h"
#include "region.h"
#include "input.h"
#include "variable.h"
#include "memory_kokkos.h"
#include "error.h"
#include "atom_masks.h"
#include "kokkos_base.h"

#include <cstring>

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
  atomKK->sync(execution_space, F_MASK | Q_MASK | MASK_MASK);

  f = atomKK->k_f.view<DeviceType>();
  q = atomKK->k_q.view<DeviceType>();
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
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagFixEfieldConstant>(0,nlocal),*this,fsum_kk);
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
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagFixEfieldNonConstant>(0,nlocal),*this,fsum_kk);
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
    const F_FLOAT qtmp = q[i];
    const F_FLOAT fx = qtmp * ex;
    const F_FLOAT fy = qtmp * ey;
    const F_FLOAT fz = qtmp * ez;
    f(i,0) += fx;
    f(i,1) += fy;
    f(i,2) += fz;
    //fsum_kk.d0 -= fx * unwrap[0] + fy * unwrap[1] + fz * unwrap[2];
    fsum_kk.d1 += fx;
    fsum_kk.d2 += fy;
    fsum_kk.d3 += fz;
  }
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixEfieldKokkos<DeviceType>::operator()(TagFixEfieldNonConstant, const int &i, double_4& fsum_kk) const {
  if (mask[i] & groupbit) {
    if (region && !d_match[i]) return;
    const F_FLOAT qtmp = q[i];
    const F_FLOAT fx = qtmp * ex;
    const F_FLOAT fy = qtmp * ey;
    const F_FLOAT fz = qtmp * ez;
    if (xstyle == ATOM) f(i,0) += d_efield(i,0);
    else if (xstyle) f(i,0) += fx;
    if (ystyle == ATOM) f(i,1) = d_efield(i,1);
    else if (ystyle) f(i,1) += fy;
    if (zstyle == ATOM) f(i,2) = d_efield(i,2);
    else if (zstyle) f(i,2) += fz;
    //fsum_kk.d0 -= fx * unwrap[0] + fy * unwrap[1] + fz * unwrap[2];
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

