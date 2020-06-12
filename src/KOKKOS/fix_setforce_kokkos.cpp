/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "fix_setforce_kokkos.h"
#include <cstring>
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
#include "region.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{NONE,CONSTANT,EQUAL,ATOM};

/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space>
FixSetForceKokkos<Space>::FixSetForceKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixSetForce(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = Space;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;

  memory->destroy(sforce);
  memoryKK->create_kokkos(k_sforce,sforce,maxatom,3,"setforce:sforce");
  d_sforce = DualViewHelper<Space>::view(k_sforce);
}

/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space>
FixSetForceKokkos<Space>::~FixSetForceKokkos()
{
  if (copymode) return;

  memoryKK->destroy_kokkos(k_sforce,sforce);
  sforce = NULL;
}

/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space>
void FixSetForceKokkos<Space>::init()
{
  FixSetForce::init();

  if (strstr(update->integrate_style,"respa"))
    error->all(FLERR,"Cannot (yet) use respa with Kokkos");
}

/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space>
void FixSetForceKokkos<Space>::post_force(int vflag)
{
  atomKK->sync(Space, X_MASK | F_MASK | MASK_MASK);

  x = DualViewHelper<Space>::view(atomKK->k_x);
  f = DualViewHelper<Space>::view(atomKK->k_f);
  mask = DualViewHelper<Space>::view(atomKK->k_mask);

  int nlocal = atom->nlocal;

  // update region if necessary

  region = NULL;
  if (iregion >= 0) {
    region = domain->regions[iregion];
    region->prematch();
    DAT::tdual_int_1d k_match = DAT::tdual_int_1d("setforce:k_match",nlocal);
    KokkosBase* regionKKBase = dynamic_cast<KokkosBase*>(region);
    regionKKBase->match_all_kokkos(groupbit,k_match);
    DualViewHelper<Space>::sync(k_match);
    d_match = DualViewHelper<Space>::view(k_match);
  }

  // reallocate sforce array if necessary

  if (varflag == ATOM && atom->nmax > maxatom) {
    maxatom = atom->nmax;
    memoryKK->destroy_kokkos(k_sforce,sforce);
    memoryKK->create_kokkos(k_sforce,sforce,maxatom,3,"setforce:sforce");
    d_sforce = DualViewHelper<Space>::view(k_sforce);
  }

  foriginal[0] = foriginal[1] = foriginal[2] = 0.0;
  KK_FLOAT_3 foriginal_kk;
  force_flag = 0;

  if (varflag == CONSTANT) {
    copymode = 1;
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagFixSetForceConstant>(0,nlocal),*this,foriginal_kk);
    copymode = 0;

  // variable force, wrap with clear/add

  } else {

    atomKK->sync(Host,ALL_MASK); // this can be removed when variable class is ported to Kokkos

    modify->clearstep_compute();

    if (xstyle == EQUAL) xvalue = input->variable->compute_equal(xvar);
    else if (xstyle == ATOM)
      input->variable->compute_atom(xvar,igroup,&sforce[0][0],3,0);
    if (ystyle == EQUAL) yvalue = input->variable->compute_equal(yvar);
    else if (ystyle == ATOM)
      input->variable->compute_atom(yvar,igroup,&sforce[0][1],3,0);
    if (zstyle == EQUAL) zvalue = input->variable->compute_equal(zvar);
    else if (zstyle == ATOM)
      input->variable->compute_atom(zvar,igroup,&sforce[0][2],3,0);

    modify->addstep_compute(update->ntimestep + 1);

    if (varflag == ATOM) {  // this can be removed when variable class is ported to Kokkos
      k_sforce.modify_host();
      DualViewHelper<Space>::sync(k_sforce);
    }

    copymode = 1;
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagFixSetForceNonConstant>(0,nlocal),*this,foriginal_kk);
    copymode = 0;
  }

  atomKK->modified(Space, F_MASK);

  foriginal[0] = foriginal_kk.d0;
  foriginal[1] = foriginal_kk.d1;
  foriginal[2] = foriginal_kk.d2;
}

template<ExecutionSpace Space>
KOKKOS_INLINE_FUNCTION
void FixSetForceKokkos<Space>::operator()(TagFixSetForceConstant, const int &i, KK_FLOAT_3& foriginal_kk) const {
  if (mask[i] & groupbit) {
    if (region && !d_match[i]) return;
    foriginal_kk.d0 += f(i,0);
    foriginal_kk.d1 += f(i,1);
    foriginal_kk.d2 += f(i,2);
    if (xstyle) f(i,0) = xvalue;
    if (ystyle) f(i,1) = yvalue;
    if (zstyle) f(i,2) = zvalue;
  }
}

template<ExecutionSpace Space>
KOKKOS_INLINE_FUNCTION
void FixSetForceKokkos<Space>::operator()(TagFixSetForceNonConstant, const int &i, KK_FLOAT_3& foriginal_kk) const {
  if (mask[i] & groupbit) {
    if (region && !d_match[i]) return;
    foriginal_kk.d0 += f(i,0);
    foriginal_kk.d1 += f(i,1);
    foriginal_kk.d2 += f(i,2);
    if (xstyle == ATOM) f(i,0) = d_sforce(i,0);
    else if (xstyle) f(i,0) = xvalue;
    if (ystyle == ATOM) f(i,1) = d_sforce(i,1);
    else if (ystyle) f(i,1) = yvalue;
    if (zstyle == ATOM) f(i,2) = d_sforce(i,2);
    else if (zstyle) f(i,2) = zvalue;
  }
}

namespace LAMMPS_NS {
template class FixSetForceKokkos<Device>;
template class FixSetForceKokkos<Host>;
}

