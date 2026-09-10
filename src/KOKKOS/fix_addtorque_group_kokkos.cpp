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

#include "fix_addtorque_group_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "domain_kokkos.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "input.h"
#include "modify.h"
#include "update.h"
#include "variable.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace FixConst;

enum{NONE,CONSTANT,EQUAL,ATOM};

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixAddTorqueGroupKokkos<DeviceType>::FixAddTorqueGroupKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixAddTorqueGroup(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  groupKK = (GroupKokkos *) group;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixAddTorqueGroupKokkos<DeviceType>::~FixAddTorqueGroupKokkos()
{
  if (copymode) return;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixAddTorqueGroupKokkos<DeviceType>::init()
{
  FixAddTorqueGroup::init();

  if (utils::strmatch(update->integrate_style,"^respa"))
    error->all(FLERR, Error::NOLASTLINE, "Cannot (yet) use respa with fix addtorque/group/kk");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixAddTorqueGroupKokkos<DeviceType>::post_force(int /*vflag*/)
{
  foriginal[0] = foriginal[1] = foriginal[2] = foriginal[3] = 0.0;
  force_flag = 0;

  // evaluate variables on host if needed

  if (varflag == EQUAL) {
    atomKK->sync(Host, ALL_MASK);
    modify->clearstep_compute();
    if (xstyle == EQUAL) xvalue = input->variable->compute_equal(xvar);
    if (ystyle == EQUAL) yvalue = input->variable->compute_equal(yvar);
    if (zstyle == EQUAL) zvalue = input->variable->compute_equal(zvar);
    modify->addstep_compute(update->ntimestep + 1);
  }

  // group properties (xcm, inertia, angmom) are reduced on the device by the
  // GroupKokkos helpers, which sync the arrays they read themselves.  omega()
  // is scalar linear algebra on the reduced quantities and has no per-atom part

  atom->check_mass(FLERR);
  const double masstotal = groupKK->mass_kk<DeviceType>(igroup);
  double xcm[3], inertia[3][3], angmom[3], omega[3];
  double tlocal[3], itorque[3], tcm[3], domegadt[3];

  groupKK->xcm_kk<DeviceType>(igroup, masstotal, xcm);
  groupKK->inertia_kk<DeviceType>(igroup, xcm, inertia);
  groupKK->angmom_kk<DeviceType>(igroup, xcm, angmom);
  group->omega(angmom, inertia, omega);

  const int nlocal = atom->nlocal;
  const double mvv2e = force->mvv2e;

  // store scalars read by the kernels below

  for (int d = 0; d < 3; d++) {
    l_xcm[d] = xcm[d];
    l_omega[d] = omega[d];
  }
  l_mvv2e = mvv2e;

  prd = Few<double,3>(domain->prd);
  h = Few<double,6>(domain->h);
  triclinic = domain->triclinic;

  atomKK->sync(execution_space, X_MASK | F_MASK | IMAGE_MASK | MASK_MASK | TYPE_MASK | RMASS_MASK);

  x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  image = atomKK->k_image.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();
  const int use_rmass = (atom->rmass != nullptr);
  if (use_rmass) {
    rmass = atomKK->k_rmass.view<DeviceType>();
  } else {
    atomKK->k_mass.sync<DeviceType>();
    mass = atomKK->k_mass.view<DeviceType>();
    type = atomKK->k_type.view<DeviceType>();
  }

  // torque the group's own rotation already exerts, reduced on the device

  tlocal[0] = tlocal[1] = tlocal[2] = 0.0;

  copymode = 1;
  if (use_rmass)
    Kokkos::parallel_reduce(
      Kokkos::RangePolicy<DeviceType,TagFixAddTorqueGroupItorqueRmass>(0,nlocal), *this,
      tlocal[0], tlocal[1], tlocal[2]);
  else
    Kokkos::parallel_reduce(
      Kokkos::RangePolicy<DeviceType,TagFixAddTorqueGroupItorqueMass>(0,nlocal), *this,
      tlocal[0], tlocal[1], tlocal[2]);
  copymode = 0;

  MPI_Allreduce(tlocal, itorque, 3, MPI_DOUBLE, MPI_SUM, world);

  tcm[0] = xvalue - mvv2e*itorque[0];
  tcm[1] = yvalue - mvv2e*itorque[1];
  tcm[2] = zvalue - mvv2e*itorque[2];
  group->omega(tcm, inertia, domegadt);

  for (int d = 0; d < 3; d++) l_domegadt[d] = domegadt[d];

  // apply forces on device

  double result[4] = {0.0, 0.0, 0.0, 0.0};

  copymode = 1;
  if (use_rmass)
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType,TagFixAddTorqueGroupRmass>(0,nlocal),
                            *this, result);
  else
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType,TagFixAddTorqueGroupMass>(0,nlocal),
                            *this, result);
  copymode = 0;

  atomKK->modified(execution_space, F_MASK);

  foriginal[0] = result[0];
  foriginal[1] = result[1];
  foriginal[2] = result[2];
  foriginal[3] = result[3];
}

/* ---------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   torque the group's own rotation already exerts about its centre of mass
   (mirrors the tlocal loop of FixAddTorqueGroup::post_force)
------------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixAddTorqueGroupKokkos<DeviceType>::operator()(TagFixAddTorqueGroupItorqueMass,
                                                     const int &i,
                                                     double &t0, double &t1, double &t2) const
{
  if (mask[i] & groupbit) {
    Few<double,3> x_i;
    x_i[0] = static_cast<double>(x(i,0));
    x_i[1] = static_cast<double>(x(i,1));
    x_i[2] = static_cast<double>(x(i,2));
    auto unwrap = DomainKokkos::unmap(prd,h,triclinic,x_i,image(i));
    const double dx = unwrap[0] - l_xcm[0];
    const double dy = unwrap[1] - l_xcm[1];
    const double dz = unwrap[2] - l_xcm[2];
    const double massone = static_cast<double>(mass[type[i]]);
    const double omegadotr = l_omega[0]*dx + l_omega[1]*dy + l_omega[2]*dz;
    t0 += massone * omegadotr * (dy*l_omega[2] - dz*l_omega[1]);
    t1 += massone * omegadotr * (dz*l_omega[0] - dx*l_omega[2]);
    t2 += massone * omegadotr * (dx*l_omega[1] - dy*l_omega[0]);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixAddTorqueGroupKokkos<DeviceType>::operator()(TagFixAddTorqueGroupItorqueRmass,
                                                     const int &i,
                                                     double &t0, double &t1, double &t2) const
{
  if (mask[i] & groupbit) {
    Few<double,3> x_i;
    x_i[0] = static_cast<double>(x(i,0));
    x_i[1] = static_cast<double>(x(i,1));
    x_i[2] = static_cast<double>(x(i,2));
    auto unwrap = DomainKokkos::unmap(prd,h,triclinic,x_i,image(i));
    const double dx = unwrap[0] - l_xcm[0];
    const double dy = unwrap[1] - l_xcm[1];
    const double dz = unwrap[2] - l_xcm[2];
    const double massone = static_cast<double>(rmass[i]);
    const double omegadotr = l_omega[0]*dx + l_omega[1]*dy + l_omega[2]*dz;
    t0 += massone * omegadotr * (dy*l_omega[2] - dz*l_omega[1]);
    t1 += massone * omegadotr * (dz*l_omega[0] - dx*l_omega[2]);
    t2 += massone * omegadotr * (dx*l_omega[1] - dy*l_omega[0]);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixAddTorqueGroupKokkos<DeviceType>::operator()(TagFixAddTorqueGroupMass,
                                                     const int &i,
                                                     value_type result) const
{
  if (mask[i] & groupbit) {
    Few<double,3> x_i;
    x_i[0] = static_cast<double>(x(i,0));
    x_i[1] = static_cast<double>(x(i,1));
    x_i[2] = static_cast<double>(x(i,2));
    auto unwrap = DomainKokkos::unmap(prd,h,triclinic,x_i,image(i));
    double dx = unwrap[0] - l_xcm[0];
    double dy = unwrap[1] - l_xcm[1];
    double dz = unwrap[2] - l_xcm[2];
    double vx = l_mvv2e*(dz*l_omega[1]-dy*l_omega[2]);
    double vy = l_mvv2e*(dx*l_omega[2]-dz*l_omega[0]);
    double vz = l_mvv2e*(dy*l_omega[0]-dx*l_omega[1]);
    const double massone = static_cast<double>(mass[type[i]]);
    double fx = massone * (dz*l_domegadt[1]-dy*l_domegadt[2] + vz*l_omega[1]-vy*l_omega[2]);
    double fy = massone * (dx*l_domegadt[2]-dz*l_domegadt[0] + vx*l_omega[2]-vz*l_omega[0]);
    double fz = massone * (dy*l_domegadt[0]-dx*l_domegadt[1] + vy*l_omega[0]-vx*l_omega[1]);
    result[0] -= fx*static_cast<double>(x(i,0)) + fy*static_cast<double>(x(i,1))
      + fz*static_cast<double>(x(i,2));
    result[1] += dy*static_cast<double>(f(i,2)) - dz*static_cast<double>(f(i,1));
    result[2] += dz*static_cast<double>(f(i,0)) - dx*static_cast<double>(f(i,2));
    result[3] += dx*static_cast<double>(f(i,1)) - dy*static_cast<double>(f(i,0));
    f(i,0) += static_cast<KK_ACC_FLOAT>(fx);
    f(i,1) += static_cast<KK_ACC_FLOAT>(fy);
    f(i,2) += static_cast<KK_ACC_FLOAT>(fz);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixAddTorqueGroupKokkos<DeviceType>::operator()(TagFixAddTorqueGroupRmass,
                                                     const int &i,
                                                     value_type result) const
{
  if (mask[i] & groupbit) {
    Few<double,3> x_i;
    x_i[0] = static_cast<double>(x(i,0));
    x_i[1] = static_cast<double>(x(i,1));
    x_i[2] = static_cast<double>(x(i,2));
    auto unwrap = DomainKokkos::unmap(prd,h,triclinic,x_i,image(i));
    double dx = unwrap[0] - l_xcm[0];
    double dy = unwrap[1] - l_xcm[1];
    double dz = unwrap[2] - l_xcm[2];
    double vx = l_mvv2e*(dz*l_omega[1]-dy*l_omega[2]);
    double vy = l_mvv2e*(dx*l_omega[2]-dz*l_omega[0]);
    double vz = l_mvv2e*(dy*l_omega[0]-dx*l_omega[1]);
    const double massone = static_cast<double>(rmass[i]);
    double fx = massone * (dz*l_domegadt[1]-dy*l_domegadt[2] + vz*l_omega[1]-vy*l_omega[2]);
    double fy = massone * (dx*l_domegadt[2]-dz*l_domegadt[0] + vx*l_omega[2]-vz*l_omega[0]);
    double fz = massone * (dy*l_domegadt[0]-dx*l_domegadt[1] + vy*l_omega[0]-vx*l_omega[1]);
    result[0] -= fx*static_cast<double>(x(i,0)) + fy*static_cast<double>(x(i,1))
      + fz*static_cast<double>(x(i,2));
    result[1] += dy*static_cast<double>(f(i,2)) - dz*static_cast<double>(f(i,1));
    result[2] += dz*static_cast<double>(f(i,0)) - dx*static_cast<double>(f(i,2));
    result[3] += dx*static_cast<double>(f(i,1)) - dy*static_cast<double>(f(i,0));
    f(i,0) += static_cast<KK_ACC_FLOAT>(fx);
    f(i,1) += static_cast<KK_ACC_FLOAT>(fy);
    f(i,2) += static_cast<KK_ACC_FLOAT>(fz);
  }
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixAddTorqueGroupKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixAddTorqueGroupKokkos<LMPHostType>;
#endif
}
