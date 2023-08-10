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

#include "fix_spring_self_kokkos.h"

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

#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixSpringSelfKokkos<DeviceType>::FixSpringSelfKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixSpringSelf(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;

  maxatom = atom->nmax;
  memory->destroy(xoriginal);
  memoryKK->create_kokkos(k_xoriginal,xoriginal,maxatom,3,"spring/self:xoriginal");
  d_xoriginal = k_xoriginal.view<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixSpringSelfKokkos<DeviceType>::~FixSpringSelfKokkos()
{
  if (copymode) return;

  memoryKK->destroy_kokkos(k_xoriginal,xoriginal);
  xoriginal = nullptr;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixSpringSelfKokkos<DeviceType>::init()
{
  FixSpringSelf::init();

  if (utils::strmatch(update->integrate_style,"^respa"))
    error->all(FLERR,"Cannot (yet) use respa with Kokkos");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixSpringSelfKokkos<DeviceType>::post_force(int /*vflag*/)
{
  atomKK->sync(execution_space, X_MASK | F_MASK | MASK_MASK);

  x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  image = atomKK->k_image.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();

  int nlocal = atom->nlocal;

  // reallocate xoriginal array if necessary

  if (atom->nmax > maxatom) {
    maxatom = atom->nmax;
    memoryKK->destroy_kokkos(k_xoriginal,xoriginal);
    memoryKK->create_kokkos(k_xoriginal,xoriginal,maxatom,3,"fix_spring/self:xoriginal");
    d_xoriginal = k_xoriginal.view<DeviceType>();
  }

  double espring_kk;


  copymode = 1;
  //Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagFixSpringSelfConstant>(0,nlocal),*this, espring_kk);
  {
  // local variables for lambda capture
  auto prd = Few<double,3>(domain->prd);
  auto h = Few<double,6>(domain->h);
  auto triclinic = domain->triclinic;
  auto l_xflag = xflag;
  auto l_yflag = yflag;
  auto l_zflag = zflag;
  auto l_k = k;
  auto l_x = x;
  auto l_xoriginal = d_xoriginal;
  auto l_f = f;
  auto l_mask = mask;
  auto l_image = image;
  auto l_groupbit = groupbit;

  Kokkos::parallel_reduce(nlocal, LAMMPS_LAMBDA(const int& i, double& espring_kk) {
    if (l_mask[i] & l_groupbit) {
      Few<double,3> x_i;
      x_i[0] = l_x(i,0);
      x_i[1] = l_x(i,1);
      x_i[2] = l_x(i,2);
      auto unwrap = DomainKokkos::unmap(prd,h,triclinic,x_i,l_image(i));
      auto dx = unwrap[0] - l_xoriginal(i, 0);
      auto dy = unwrap[1] - l_xoriginal(i, 1);
      auto dz = unwrap[2] - l_xoriginal(i, 2);
      if (!l_xflag) dx = 0.0;
      if (!l_yflag) dy = 0.0;
      if (!l_zflag) dz = 0.0;
      l_f(i,0) -= l_k*dx;
      l_f(i,1) -= l_k*dy;
      l_f(i,2) -= l_k*dz;
      espring_kk += l_k * (dx*dx + dy*dy + dz*dz);
    }
  },espring_kk);
  }

  copymode = 0;

  atomKK->modified(execution_space, F_MASK);

  espring = 0.5*espring_kk;
}

namespace LAMMPS_NS {
template class FixSpringSelfKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixSpringSelfKokkos<LMPHostType>;
#endif
}

