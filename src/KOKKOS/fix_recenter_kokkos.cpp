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
   Contributing author: Mitch Murphy (alphataubio.com)
 ------------------------------------------------------------------------- */

#include "fix_recenter_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "input.h"
#include "modify.h"
#include "update.h"
#include "domain.h"
#include "group.h"
#include "kokkos_few.h"

using namespace LAMMPS_NS;

enum{BOX,LATTICE,FRACTION};

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixRecenterKokkos<DeviceType>::FixRecenterKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixRecenter(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *)atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  //datamask_read = X_MASK | F_MASK | RMASS_MASK | MASK_MASK | TYPE_MASK;
  datamask_read = X_MASK | MASK_MASK;
  datamask_modify = X_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRecenterKokkos<DeviceType>::initial_integrate(int /*vflag*/)
{

utils::logmesg(lmp, "ok 2a\n");

  atomKK->sync(execution_space,datamask_read);
  atomKK->modified(execution_space,datamask_modify);

  x = atomKK->k_x.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();
  int nlocal = atomKK->nlocal;
  if (igroup == atomKK->firstgroup) nlocal = atomKK->nfirst;

  // FIX RECENTER
  // target COM
  // bounding box around domain works for both orthogonal and triclinic

  double xtarget = xinit;
  double ytarget = yinit;
  double ztarget = zinit;

  xflag=yflag=zflag=1;

  utils::logmesg(lmp, "ok 2b\n");

  // FIXME: only supported in KOKKOS...
  // fix ID group-ID recenter INIT INIT INIT shift all

  /*
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

  */

  // current COM

  // FIXME: make Group kokkos-aware
  //double xcm[3];
  //if (group->dynamic[igroup])
  //  masstotal = group->mass(igroup);

  //group->xcm(igroup,masstotal,xcm);

  /* this is needed because Group is not Kokkos-aware ! */
  atomKK->sync(ExecutionSpaceFromDevice<LMPHostType>::space,X_MASK);
  Few<double, 3> tmpxcm;
  group->xcm(igroup,masstotal,&tmpxcm[0]);
  const Few<double, 3> xcm(tmpxcm);


  utils::logmesg(lmp, "ok 2c, xcm={},{},{}\n", xcm[0], xcm[1], xcm[2]);

  // shift coords by difference between actual COM and requested COM

  double shiftx = xflag ? (xtarget - xcm[0]) : 0.0;
  double shifty = yflag ? (ytarget - xcm[1]) : 0.0;
  double shiftz = zflag ? (ztarget - xcm[2]) : 0.0;
  distance = sqrt(shiftx*shiftx + shifty*shifty + shiftz*shiftz);

utils::logmesg(lmp, "ok 2d, shift={},{},{}\n", shiftx, shifty, shiftz);
  // ----

  copymode = 1;

  auto group2bit_copy = group2bit;

  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType>(0,nlocal),
    LAMMPS_LAMBDA(int i) {
      if (mask[i] & group2bit_copy) {
        x(i,0) += shiftx;
        x(i,1) += shifty;
        x(i,2) += shiftz;
        utils::logmesg(lmp, "x({})={} {} {}\n", i, x(i,0), x(i,1), x(i,2));

      }
    });

      utils::logmesg(lmp, "x(1)={} {} {}\n", x(1,0), x(1,1), x(1,2));

  copymode = 0;
}


namespace LAMMPS_NS {
template class FixRecenterKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixRecenterKokkos<LMPHostType>;
#endif
}
