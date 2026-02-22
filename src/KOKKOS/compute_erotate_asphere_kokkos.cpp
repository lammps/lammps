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

#include "compute_erotate_asphere_kokkos.h"
#include "math_extra_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "update.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
ComputeERotateAsphereKokkos<DeviceType>::ComputeERotateAsphereKokkos(LAMMPS *lmp, int narg, char **arg) :
  ComputeERotateAsphere(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  datamask_read = OMEGA_MASK | ANGMOM_MASK | MASK_MASK | RMASS_MASK | ELLIPSOID_MASK | BONUS_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ComputeERotateAsphereKokkos<DeviceType>::init()
{
  ComputeERotateAsphere::init();

  // Vanilla has if statements here checking for ellipsoid, line, tri
  // Kokkos version only supports ellipsoid this now, so this part of the if statement is not needed
  if (avec_line || avec_tri) error->all(FLERR,"No Kokkos implementation for line or tri atom styles");
  avecEllipKK = dynamic_cast<AtomVecEllipsoidKokkos *>(atom->style_match("ellipsoid"));
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
double ComputeERotateAsphereKokkos<DeviceType>::compute_scalar()
{
  atomKK->sync(execution_space,datamask_read);

  invoked_scalar = update->ntimestep;

  omega = atomKK->k_omega.view<DeviceType>();
  angmom = atomKK->k_angmom.view<DeviceType>();
  rmass = atomKK->k_rmass.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();
  int nlocal = atom->nlocal;

  ellipsoid = atomKK->k_ellipsoid.view<DeviceType>();
  bonus = avecEllipKK->k_bonus.view<DeviceType>();

  // sum rotational energy for each particle
  // point particles will not contribute, due to radius = 0.0

  double erotate = 0.0;

  {
    // local variables for lambda capture

    auto l_angmom = angmom;
    auto l_rmass = rmass;
    auto l_mask = mask;
    auto l_groupbit = groupbit;

    auto l_ellipsoid = ellipsoid;
    auto l_bonus = bonus;

    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType>(0,nlocal), LAMMPS_LAMBDA(int i, double &erotate) {
      if (l_mask[i] & l_groupbit) {
        // Vanilla has if statements here checking for ellipsoid, line, tri
        // Kokkos version only supports ellipsoid this now, so this part of the if statement is not needed
        if (l_ellipsoid(i) >= 0) {
          KK_FLOAT inertia[3], wbody[3];
          KK_FLOAT rot[3][3];
          auto shape = l_bonus(l_ellipsoid(i)).shape;
          auto quat = l_bonus(l_ellipsoid(i)).quat;

          // principal moments of inertia

          inertia[0] = l_rmass(i) * ( (shape[1]*shape[1] + shape[2]*shape[2])/5.0 );
          inertia[1] = l_rmass(i) * ( (shape[0]*shape[0] + shape[2]*shape[2])/5.0 );
          inertia[2] = l_rmass(i) * ( (shape[0]*shape[0] + shape[1]*shape[1])/5.0 );

          // wbody = angular velocity in body frame

          MathExtraKokkos::quat_to_mat(quat,rot);
          KK_FLOAT angmom_vec[3] = {l_angmom(i,0), l_angmom(i,1), l_angmom(i,2)};
          MathExtraKokkos::transpose_matvec(rot,angmom_vec,wbody);
          wbody[0] /= inertia[0];
          wbody[1] /= inertia[1];
          wbody[2] /= inertia[2];

        erotate +=
            static_cast<double>((inertia[0] * wbody[0] * wbody[0]) +
                         (inertia[1] * wbody[1] * wbody[1]) +
                         (inertia[2] * wbody[2] * wbody[2]));
        }
      }
    },erotate);
  }

  MPI_Allreduce(&erotate, &scalar, 1, MPI_DOUBLE, MPI_SUM, world);
  scalar *= pfactor;
  return scalar;
}

namespace LAMMPS_NS {
template class ComputeERotateAsphereKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class ComputeERotateAsphereKokkos<LMPHostType>;
#endif
}
