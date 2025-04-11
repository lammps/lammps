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
   Contributing author: Mitch Murphy (alphataubio at gmail)
------------------------------------------------------------------------- */

#include "fix_recenter_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "input.h"
#include "modify.h"
#include "update.h"
#include "domain.h"
#include "group_kokkos.h"

using namespace LAMMPS_NS;

enum{BOX,LATTICE,FRACTION};

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixRecenterKokkos<DeviceType>::FixRecenterKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixRecenter(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *)atom;
  groupKK = (GroupKokkos *)group;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  datamask_read = X_MASK | MASK_MASK;
  datamask_modify = X_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRecenterKokkos<DeviceType>::initial_integrate(int /*vflag*/)
{
  atomKK->sync(execution_space,datamask_read);
  int nlocal = atomKK->nlocal;
  if (igroup == atomKK->firstgroup) nlocal = atomKK->nfirst;

  // target COM
  // bounding box around domain works for both orthogonal and triclinic

  double xtarget,ytarget,ztarget;
  double *bboxlo,*bboxhi;

  if (scaleflag == FRACTION) {
    if (domain->triclinic == 0) {
      bboxlo = domain->boxlo;
      bboxhi = domain->boxhi;
    } else {
      bboxlo = domain->boxlo_bound;
      bboxhi = domain->boxhi_bound;
    }
  }

  if (xinitflag) xtarget = xinit;
  else if (scaleflag == FRACTION)
    xtarget = bboxlo[0] + xcom*(bboxhi[0] - bboxlo[0]);
  else xtarget = xcom;

  if (yinitflag) ytarget = yinit;
  else if (scaleflag == FRACTION)
    ytarget = bboxlo[1] + ycom*(bboxhi[1] - bboxlo[1]);
  else ytarget = ycom;

  if (zinitflag) ztarget = zinit;
  else if (scaleflag == FRACTION)
    ztarget = bboxlo[2] + zcom*(bboxhi[2] - bboxlo[2]);
  else ztarget = zcom;

  // current COM


  if (group->dynamic[igroup]) masstotal = groupKK->mass_kk<DeviceType>(igroup);
  double xcm[3];
  groupKK->xcm_kk<DeviceType>(igroup,masstotal,xcm);

  // shift coords by difference between actual COM and requested COM

  shift[0] = xflag ? (xtarget - xcm[0]) : 0.0;
  shift[1] = yflag ? (ytarget - xcm[1]) : 0.0;
  shift[2] = zflag ? (ztarget - xcm[2]) : 0.0;
  distance = sqrt(shift[0]*shift[0] + shift[1]*shift[1] + shift[2]*shift[2]);

  auto d_x = atomKK->k_x.template view<DeviceType>();
  auto d_mask = atomKK->k_mask.template view<DeviceType>();
  auto l_group2bit = group2bit;
  double l_shiftx = shift[0];
  double l_shifty = shift[1];
  double l_shiftz = shift[2];

  copymode = 1;

  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType>(0,nlocal),
    KOKKOS_LAMBDA(const int i) {
      if (d_mask[i] & l_group2bit) {
        d_x(i,0) += l_shiftx;
        d_x(i,1) += l_shifty;
        d_x(i,2) += l_shiftz;
      }
    });

  copymode = 0;
  atomKK->modified(execution_space,datamask_modify);
}

namespace LAMMPS_NS {
template class FixRecenterKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixRecenterKokkos<LMPHostType>;
#endif
}
