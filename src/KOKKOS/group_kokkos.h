/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_GROUP_KOKKOS_H
#define LMP_GROUP_KOKKOS_H

#include "group.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "domain_kokkos.h"
#include "kokkos_few.h"
#include "kokkos_type.h"


namespace LAMMPS_NS {

class GroupKokkos : public Group {
 public:
  GroupKokkos(LAMMPS *lmp) : Group(lmp) { atomKK = (AtomKokkos *)atom; }

// ----------------------------------------------------------------------
// computations on a group of atoms
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   compute the total mass of group of atoms
   use either per-type mass or per-atom rmass
------------------------------------------------------------------------- */

template<class DeviceType>
double mass_kk(int igroup)
{
  auto execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  int groupbit = bitmask[igroup];
  auto d_mask = atomKK->k_mask.template view<DeviceType>();
  double one = 0.0;

  if (atomKK->rmass) {

    auto d_rmass = atomKK->k_rmass.template view<DeviceType>();
    atomKK->sync(execution_space,MASK_MASK|RMASS_MASK);

    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType>(0,atom->nlocal), KOKKOS_LAMBDA(const int i, double &l_one) {
      if (d_mask(i) & groupbit) l_one += d_rmass(i);
    }, one);

  } else {

    auto d_mass = atomKK->k_mass.template view<DeviceType>();
    auto d_type = atomKK->k_type.template view<DeviceType>();
    atomKK->sync(execution_space,MASK_MASK|TYPE_MASK);
    atomKK->k_mass.template sync<DeviceType>();

    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType>(0,atom->nlocal), KOKKOS_LAMBDA(const int i, double &l_one) {
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
void xcm_kk(int igroup, double masstotal, double *xcm)
{
  auto execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  int groupbit = bitmask[igroup];
  auto d_x = atomKK->k_x.template view<DeviceType>();
  auto d_mask = atomKK->k_mask.template view<DeviceType>();
  auto d_image = atomKK->k_image.template view<DeviceType>();
  auto l_prd = Few<double, 3>(domain->prd);
  auto l_h = Few<double, 6>(domain->h);
  auto l_triclinic = domain->triclinic;
  double cmone[3] = {0.0, 0.0, 0.0};

  if (atomKK->rmass) {

    auto d_rmass = atomKK->k_rmass.template view<DeviceType>();
    atomKK->sync(execution_space,X_MASK|MASK_MASK|IMAGE_MASK|RMASS_MASK);

    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType>(0,atom->nlocal), KOKKOS_LAMBDA(const int i, double &l_cmx, double &l_cmy, double &l_cmz) {
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
    auto d_type = atomKK->k_type.template view<DeviceType>();
    atomKK->sync(execution_space,X_MASK|MASK_MASK|IMAGE_MASK|TYPE_MASK);
    atomKK->k_mass.template sync<DeviceType>();

    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType>(0,atom->nlocal), KOKKOS_LAMBDA(const int i, double &l_cmx, double &l_cmy, double &l_cmz) {
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

  MPI_Allreduce(cmone, xcm, 3, MPI_DOUBLE, MPI_SUM, world);
  if (masstotal > 0.0) {
    xcm[0] /= masstotal;
    xcm[1] /= masstotal;
    xcm[2] /= masstotal;
  }
}

/* ----------------------------------------------------------------------
   compute the center-of-mass velocity of group of atoms
   masstotal = total mass
   return center-of-mass velocity in vcm[]
------------------------------------------------------------------------- */

template<class DeviceType>
void vcm_kk(int igroup, double masstotal, double *vcm)
{
  auto execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  int groupbit = bitmask[igroup];
  auto d_v = atomKK->k_v.template view<DeviceType>();
  auto d_mask = atomKK->k_mask.template view<DeviceType>();
  auto d_image = atomKK->k_image.template view<DeviceType>();
  double p[3] = {0.0, 0.0, 0.0};

  if (atomKK->rmass) {

    auto d_rmass = atomKK->k_rmass.template view<DeviceType>();
    atomKK->sync(execution_space,V_MASK|MASK_MASK|IMAGE_MASK|RMASS_MASK);

    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType>(0,atom->nlocal), KOKKOS_LAMBDA(const int i, double &l_px, double &l_py, double &l_pz) {
      if (d_mask(i) & groupbit) {
        double massone = d_rmass(i);
        l_px += d_v(i,0) * massone;
        l_py += d_v(i,1) * massone;
        l_pz += d_v(i,2) * massone;
      }
    }, p[0], p[1], p[2]);

  } else {

    auto d_mass = atomKK->k_mass.template view<DeviceType>();
    auto d_type = atomKK->k_type.template view<DeviceType>();
    atomKK->sync(execution_space,V_MASK|MASK_MASK|IMAGE_MASK|TYPE_MASK);
    atomKK->k_mass.template sync<DeviceType>();

    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType>(0,atom->nlocal), KOKKOS_LAMBDA(const int i, double &l_px, double &l_py, double &l_pz) {
      if (d_mask(i) & groupbit) {
        double massone = d_mass(d_type(i));
        l_px += d_v(i,0) * massone;
        l_py += d_v(i,1) * massone;
        l_pz += d_v(i,2) * massone;
      }
    }, p[0], p[1], p[2]);

  }

  MPI_Allreduce(p, vcm, 3, MPI_DOUBLE, MPI_SUM, world);
  if (masstotal > 0.0) {
    vcm[0] /= masstotal;
    vcm[1] /= masstotal;
    vcm[2] /= masstotal;
  }
}

/* ----------------------------------------------------------------------
   compute the angular momentum L (lmom) of group
   around center-of-mass cm
   must unwrap atoms to compute L correctly
------------------------------------------------------------------------- */

template<class DeviceType>
void angmom_kk(int igroup, double *xcm, double *lmom)
{
  auto execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  int groupbit = bitmask[igroup];
  auto d_x = atomKK->k_x.template view<DeviceType>();
  auto d_v = atomKK->k_v.template view<DeviceType>();
  auto d_mask = atomKK->k_mask.template view<DeviceType>();
  auto d_image = atomKK->k_image.template view<DeviceType>();
  auto l_prd = Few<double, 3>(domain->prd);
  auto l_h = Few<double, 6>(domain->h);
  auto l_triclinic = domain->triclinic;
  auto l_xcm0 = xcm[0];
  auto l_xcm1 = xcm[1];
  auto l_xcm2 = xcm[2];
  double p[3] = {0.0, 0.0, 0.0};

  if (atomKK->rmass) {

    auto d_rmass = atomKK->k_rmass.template view<DeviceType>();
    atomKK->sync(execution_space,X_MASK|V_MASK|MASK_MASK|IMAGE_MASK|RMASS_MASK);

    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType>(0,atom->nlocal), KOKKOS_LAMBDA(const int i, double &l_px, double &l_py, double &l_pz) {
      if (d_mask(i) & groupbit) {
        double massone = d_rmass(i);
        Few<double,3> x_i;
        x_i[0] = d_x(i,0);
        x_i[1] = d_x(i,1);
        x_i[2] = d_x(i,2);
        auto unwrapKK = DomainKokkos::unmap(l_prd,l_h,l_triclinic,x_i,d_image(i));
        double dx = unwrapKK[0] - l_xcm0;
        double dy = unwrapKK[1] - l_xcm1;
        double dz = unwrapKK[2] - l_xcm2;
        l_px += massone * (dy * d_v(i,2) - dz * d_v(i,1));
        l_py += massone * (dz * d_v(i,0) - dx * d_v(i,2));
        l_pz += massone * (dx * d_v(i,1) - dy * d_v(i,0));
      }
    }, p[0], p[1], p[2]);

  } else {

    auto d_mass = atomKK->k_mass.template view<DeviceType>();
    auto d_type = atomKK->k_type.template view<DeviceType>();
    atomKK->sync(execution_space,X_MASK|V_MASK|MASK_MASK|IMAGE_MASK|TYPE_MASK);
    atomKK->k_mass.template sync<DeviceType>();

    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType>(0,atom->nlocal), KOKKOS_LAMBDA(const int i, double &l_px, double &l_py, double &l_pz) {
      if (d_mask(i) & groupbit) {
        double massone = d_mass(d_type(i));
        Few<double,3> x_i;
        x_i[0] = d_x(i,0);
        x_i[1] = d_x(i,1);
        x_i[2] = d_x(i,2);
        auto unwrapKK = DomainKokkos::unmap(l_prd,l_h,l_triclinic,x_i,d_image(i));
        double dx = unwrapKK[0] - l_xcm0;
        double dy = unwrapKK[1] - l_xcm1;
        double dz = unwrapKK[2] - l_xcm2;
        l_px += massone * (dy * d_v(i,2) - dz * d_v(i,1));
        l_py += massone * (dz * d_v(i,0) - dx * d_v(i,2));
        l_pz += massone * (dx * d_v(i,1) - dy * d_v(i,0));
      }
    }, p[0], p[1], p[2]);

  }
  MPI_Allreduce(p, lmom, 3, MPI_DOUBLE, MPI_SUM, world);
}

/* ----------------------------------------------------------------------
   compute moment of inertia tensor around center-of-mass xcm of group
   must unwrap atoms to compute itensor correctly
------------------------------------------------------------------------- */

template<class DeviceType>
void inertia_kk(int igroup, double *xcm, double itensor[3][3])
{
  auto execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  int groupbit = bitmask[igroup];
  auto d_x = atomKK->k_x.template view<DeviceType>();
  auto d_mask = atomKK->k_mask.template view<DeviceType>();
  auto d_image = atomKK->k_image.template view<DeviceType>();
  auto l_prd = Few<double, 3>(domain->prd);
  auto l_h = Few<double, 6>(domain->h);
  auto l_triclinic = domain->triclinic;
  auto l_xcm0 = xcm[0];
  auto l_xcm1 = xcm[1];
  auto l_xcm2 = xcm[2];

  double ione[3][3];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) ione[i][j] = 0.0;

  if (atomKK->rmass) {

    auto d_rmass = atomKK->k_rmass.template view<DeviceType>();
    atomKK->sync(execution_space,X_MASK|MASK_MASK|IMAGE_MASK|RMASS_MASK);

    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType>(0,atom->nlocal), KOKKOS_LAMBDA(const int i, double &l_i00, double &l_i11, double &l_i22, double &l_i01, double &l_i12, double &l_i02) {
      if (d_mask(i) & groupbit) {
        double massone = d_rmass(i);
        Few<double,3> x_i;
        x_i[0] = d_x(i,0);
        x_i[1] = d_x(i,1);
        x_i[2] = d_x(i,2);
        auto unwrapKK = DomainKokkos::unmap(l_prd,l_h,l_triclinic,x_i,d_image(i));
        double dx = unwrapKK[0] - l_xcm0;
        double dy = unwrapKK[1] - l_xcm1;
        double dz = unwrapKK[2] - l_xcm2;
        l_i00 += massone * (dy * dy + dz * dz);
        l_i11 += massone * (dx * dx + dz * dz);
        l_i22 += massone * (dx * dx + dy * dy);
        l_i01 -= massone * dx * dy;
        l_i12 -= massone * dy * dz;
        l_i02 -= massone * dx * dz;
      }
    }, ione[0][0], ione[1][1], ione[2][2], ione[0][1], ione[1][2], ione[0][2]);

  } else {

    auto d_mass = atomKK->k_mass.template view<DeviceType>();
    auto d_type = atomKK->k_type.template view<DeviceType>();
    atomKK->sync(execution_space,X_MASK|MASK_MASK|IMAGE_MASK|TYPE_MASK);
    atomKK->k_mass.template sync<DeviceType>();

    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType>(0,atom->nlocal), KOKKOS_LAMBDA(const int i, double &l_i00, double &l_i11, double &l_i22, double &l_i01, double &l_i12, double &l_i02) {
      if (d_mask(i) & groupbit) {
        double massone = d_mass(d_type(i));
        Few<double,3> x_i;
        x_i[0] = d_x(i,0);
        x_i[1] = d_x(i,1);
        x_i[2] = d_x(i,2);
        auto unwrapKK = DomainKokkos::unmap(l_prd,l_h,l_triclinic,x_i,d_image(i));
        double dx = unwrapKK[0] - l_xcm0;
        double dy = unwrapKK[1] - l_xcm1;
        double dz = unwrapKK[2] - l_xcm2;
        l_i00 += massone * (dy * dy + dz * dz);
        l_i11 += massone * (dx * dx + dz * dz);
        l_i22 += massone * (dx * dx + dy * dy);
        l_i01 -= massone * dx * dy;
        l_i12 -= massone * dy * dz;
        l_i02 -= massone * dx * dz;
      }
    }, ione[0][0], ione[1][1], ione[2][2], ione[0][1], ione[1][2], ione[0][2]);

  }

  ione[1][0] = ione[0][1];
  ione[2][1] = ione[1][2];
  ione[2][0] = ione[0][2];
  MPI_Allreduce(&ione[0][0], &itensor[0][0], 9, MPI_DOUBLE, MPI_SUM, world);
}

};

}    // namespace LAMMPS_NS

#endif
