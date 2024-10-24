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

#include "group_kokkos.h"

#include "atom_kokkos.h"
#include "domain_kokkos.h"
#include "kokkos_few.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
GroupKokkos<DeviceType>::GroupKokkos(LAMMPS *lmp) : Group(lmp)
{
  atomKK = (AtomKokkos *)atom;
}

// ----------------------------------------------------------------------
// computations on a group of atoms
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   compute the total mass of group of atoms
   use either per-type mass or per-atom rmass
------------------------------------------------------------------------- */

template<class DeviceType>
double GroupKokkos<DeviceType>::mass(int igroup)
{
  int groupbit = bitmask[igroup];

  auto d_mass = atomKK->k_mass.template view<DeviceType>();
  auto d_rmass = atomKK->k_rmass.template view<DeviceType>();
  auto d_mask = atomKK->k_mask.template view<DeviceType>();
  auto d_type = atomKK->k_type.template view<DeviceType>();

  double one = 0.0;

  if (atomKK->rmass) {

    Kokkos::parallel_reduce(atom->nlocal, KOKKOS_LAMBDA(const int i, double &l_one) {
      if (d_mask(i) & groupbit) l_one += d_rmass(i);
    }, one);

  } else {

    Kokkos::parallel_reduce(atom->nlocal, KOKKOS_LAMBDA(const int i, double &l_one) {
      if (d_mask(i) & groupbit) l_one += d_mass(d_type(i));
    }, one);

  }

  double all;
  MPI_Allreduce(&one, &all, 1, MPI_DOUBLE, MPI_SUM, world);
  return all;
}

/* ----------------------------------------------------------------------
   compute the center-of-mass coords of group of atoms
   masstotal = total mass
   return center-of-mass coords in cm[]
   must unwrap atoms to compute center-of-mass correctly
------------------------------------------------------------------------- */

template<class DeviceType>
void GroupKokkos<DeviceType>::xcm(int igroup, double masstotal, double *cm)
{
  int groupbit = bitmask[igroup];

  auto d_x = atomKK->k_x.template view<DeviceType>();
  auto d_mask = atomKK->k_mask.template view<DeviceType>();
  auto d_type = atomKK->k_type.template view<DeviceType>();
  auto d_image = atomKK->k_image.template view<DeviceType>();
  auto l_prd = Few<double, 3>(domain->prd);
  auto l_h = Few<double, 6>(domain->h);
  auto l_triclinic = domain->triclinic;

  double cmone[3];

  if (atomKK->rmass) {

    auto d_rmass = atomKK->k_rmass.template view<DeviceType>();

    Kokkos::parallel_reduce(atom->nlocal, KOKKOS_LAMBDA(const int i, double &l_cmx, double &l_cmy, double &l_cmz) {
      if (d_mask(i) & groupbit) {
        double massone = d_rmass(i);
        Few<double,3> x_i;
        x_i[0] = d_x(i,0);
        x_i[1] = d_x(i,1);
        x_i[2] = d_x(i,2);
        auto unwrapKK = DomainKokkos::unmap(l_prd,l_h,l_triclinic,x_i,d_image(i));
        l_cmx += unwrapKK[0] * massone;
        l_cmy += unwrapKK[1] * massone;
        l_cmz += unwrapKK[2] * massone;
      }
    }, cmone[0], cmone[1], cmone[2]);

  } else {

    auto d_mass = atomKK->k_mass.template view<DeviceType>();

    Kokkos::parallel_reduce(atom->nlocal, KOKKOS_LAMBDA(const int i, double &l_cmx, double &l_cmy, double &l_cmz) {
      if (d_mask(i) & groupbit) {
        double massone = d_mass(d_type(i));
        Few<double,3> x_i;
        x_i[0] = d_x(i,0);
        x_i[1] = d_x(i,1);
        x_i[2] = d_x(i,2);
        auto unwrapKK = DomainKokkos::unmap(l_prd,l_h,l_triclinic,x_i,d_image(i));
        l_cmx += unwrapKK[0] * massone;
        l_cmy += unwrapKK[1] * massone;
        l_cmz += unwrapKK[2] * massone;
      }
    }, cmone[0], cmone[1], cmone[2]);

  }

  MPI_Allreduce(cmone, cm, 3, MPI_DOUBLE, MPI_SUM, world);
  if (masstotal > 0.0) {
    cm[0] /= masstotal;
    cm[1] /= masstotal;
    cm[2] /= masstotal;
  }
}

namespace LAMMPS_NS {
template class GroupKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class GroupKokkos<LMPHostType>;
#endif
}
