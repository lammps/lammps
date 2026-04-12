/* -*- c++ -*- ----------------------------------------------------------
   Shared Kokkos rigid body storage and legacy transforms for
   fix_rigid_small_kokkos and fix_rigid_nh_small_kokkos.
------------------------------------------------------------------------- */

#ifndef LMP_RIGID_BODY_KOKKOS_HPP
#define LMP_RIGID_BODY_KOKKOS_HPP

#include "fix_rigid_small.h"

using namespace LAMMPS_NS;

struct BodyKokkos {
  int natoms;
  int ilocal;
  KK_FLOAT mass;
  KK_FLOAT xcm[3];
  KK_FLOAT xgc[3];
  KK_FLOAT vcm[3];
  KK_FLOAT fcm[3];
  KK_FLOAT torque[3];
  KK_FLOAT quat[4];
  KK_FLOAT inertia[3];
  KK_FLOAT ex_space[3];
  KK_FLOAT ey_space[3];
  KK_FLOAT ez_space[3];
  KK_FLOAT xgc_body[3];
  KK_FLOAT angmom[3];
  KK_FLOAT omega[3];
  KK_FLOAT conjqm[4];
  int remapflag[4];
  imageint image;
  imageint dummy;

  // Default constructor (required since we are declaring custom constructors)
  KOKKOS_INLINE_FUNCTION
  BodyKokkos() = default;

  // Constructor from legacy LAMMPS Body (fixes static_cast<BodyKokkos>(b))
  KOKKOS_INLINE_FUNCTION
  BodyKokkos(const FixRigidSmall::Body &b) {
    natoms = b.natoms;
    ilocal = b.ilocal;
    mass   = static_cast<KK_FLOAT>(b.mass);

    xcm[0] = static_cast<KK_FLOAT>(b.xcm[0]); xcm[1] = static_cast<KK_FLOAT>(b.xcm[1]); xcm[2] = static_cast<KK_FLOAT>(b.xcm[2]);
    xgc[0] = static_cast<KK_FLOAT>(b.xgc[0]); xgc[1] = static_cast<KK_FLOAT>(b.xgc[1]); xgc[2] = static_cast<KK_FLOAT>(b.xgc[2]);
    vcm[0] = static_cast<KK_FLOAT>(b.vcm[0]); vcm[1] = static_cast<KK_FLOAT>(b.vcm[1]); vcm[2] = static_cast<KK_FLOAT>(b.vcm[2]);
    fcm[0] = static_cast<KK_FLOAT>(b.fcm[0]); fcm[1] = static_cast<KK_FLOAT>(b.fcm[1]); fcm[2] = static_cast<KK_FLOAT>(b.fcm[2]);
    torque[0] = static_cast<KK_FLOAT>(b.torque[0]); torque[1] = static_cast<KK_FLOAT>(b.torque[1]); torque[2] = static_cast<KK_FLOAT>(b.torque[2]);
    quat[0] = static_cast<KK_FLOAT>(b.quat[0]); quat[1] = static_cast<KK_FLOAT>(b.quat[1]); quat[2] = static_cast<KK_FLOAT>(b.quat[2]); quat[3] = static_cast<KK_FLOAT>(b.quat[3]);
    inertia[0] = static_cast<KK_FLOAT>(b.inertia[0]); inertia[1] = static_cast<KK_FLOAT>(b.inertia[1]); inertia[2] = static_cast<KK_FLOAT>(b.inertia[2]);
    ex_space[0] = static_cast<KK_FLOAT>(b.ex_space[0]); ex_space[1] = static_cast<KK_FLOAT>(b.ex_space[1]); ex_space[2] = static_cast<KK_FLOAT>(b.ex_space[2]);
    ey_space[0] = static_cast<KK_FLOAT>(b.ey_space[0]); ey_space[1] = static_cast<KK_FLOAT>(b.ey_space[1]); ey_space[2] = static_cast<KK_FLOAT>(b.ey_space[2]);
    ez_space[0] = static_cast<KK_FLOAT>(b.ez_space[0]); ez_space[1] = static_cast<KK_FLOAT>(b.ez_space[1]); ez_space[2] = static_cast<KK_FLOAT>(b.ez_space[2]);
    xgc_body[0] = static_cast<KK_FLOAT>(b.xgc_body[0]); xgc_body[1] = static_cast<KK_FLOAT>(b.xgc_body[1]); xgc_body[2] = static_cast<KK_FLOAT>(b.xgc_body[2]);
    angmom[0] = static_cast<KK_FLOAT>(b.angmom[0]); angmom[1] = static_cast<KK_FLOAT>(b.angmom[1]); angmom[2] = static_cast<KK_FLOAT>(b.angmom[2]);
    omega[0] = static_cast<KK_FLOAT>(b.omega[0]); omega[1] = static_cast<KK_FLOAT>(b.omega[1]); omega[2] = static_cast<KK_FLOAT>(b.omega[2]);
    conjqm[0] = static_cast<KK_FLOAT>(b.conjqm[0]); conjqm[1] = static_cast<KK_FLOAT>(b.conjqm[1]); conjqm[2] = static_cast<KK_FLOAT>(b.conjqm[2]); conjqm[3] = static_cast<KK_FLOAT>(b.conjqm[3]);

    remapflag[0] = b.remapflag[0]; remapflag[1] = b.remapflag[1]; remapflag[2] = b.remapflag[2]; remapflag[3] = b.remapflag[3];
    image = b.image;
    dummy = b.dummy;
  }

  // Conversion operator back to legacy LAMMPS Body (fixes static_cast<FixRigidSmall::Body>(bk))
  KOKKOS_INLINE_FUNCTION
  operator FixRigidSmall::Body() const {
    FixRigidSmall::Body b;
    b.natoms = natoms;
    b.ilocal = ilocal;
    b.mass   = static_cast<double>(mass);

    b.xcm[0] = static_cast<double>(xcm[0]); b.xcm[1] = static_cast<double>(xcm[1]); b.xcm[2] = static_cast<double>(xcm[2]);
    b.xgc[0] = static_cast<double>(xgc[0]); b.xgc[1] = static_cast<double>(xgc[1]); b.xgc[2] = static_cast<double>(xgc[2]);
    b.vcm[0] = static_cast<double>(vcm[0]); b.vcm[1] = static_cast<double>(vcm[1]); b.vcm[2] = static_cast<double>(vcm[2]);
    b.fcm[0] = static_cast<double>(fcm[0]); b.fcm[1] = static_cast<double>(fcm[1]); b.fcm[2] = static_cast<double>(fcm[2]);
    b.torque[0] = static_cast<double>(torque[0]); b.torque[1] = static_cast<double>(torque[1]); b.torque[2] = static_cast<double>(torque[2]);
    b.quat[0] = static_cast<double>(quat[0]); b.quat[1] = static_cast<double>(quat[1]); b.quat[2] = static_cast<double>(quat[2]); b.quat[3] = static_cast<double>(quat[3]);
    b.inertia[0] = static_cast<double>(inertia[0]); b.inertia[1] = static_cast<double>(inertia[1]); b.inertia[2] = static_cast<double>(inertia[2]);
    b.ex_space[0] = static_cast<double>(ex_space[0]); b.ex_space[1] = static_cast<double>(ex_space[1]); b.ex_space[2] = static_cast<double>(ex_space[2]);
    b.ey_space[0] = static_cast<double>(ey_space[0]); b.ey_space[1] = static_cast<double>(ey_space[1]); b.ey_space[2] = static_cast<double>(ey_space[2]);
    b.ez_space[0] = static_cast<double>(ez_space[0]); b.ez_space[1] = static_cast<double>(ez_space[1]); b.ez_space[2] = static_cast<double>(ez_space[2]);
    b.xgc_body[0] = static_cast<double>(xgc_body[0]); b.xgc_body[1] = static_cast<double>(xgc_body[1]); b.xgc_body[2] = static_cast<double>(xgc_body[2]);
    b.angmom[0] = static_cast<double>(angmom[0]); b.angmom[1] = static_cast<double>(angmom[1]); b.angmom[2] = static_cast<double>(angmom[2]);
    b.omega[0] = static_cast<double>(omega[0]); b.omega[1] = static_cast<double>(omega[1]); b.omega[2] = static_cast<double>(omega[2]);
    b.conjqm[0] = static_cast<double>(conjqm[0]); b.conjqm[1] = static_cast<double>(conjqm[1]); b.conjqm[2] = static_cast<double>(conjqm[2]); b.conjqm[3] = static_cast<double>(conjqm[3]);

    b.remapflag[0] = remapflag[0]; b.remapflag[1] = remapflag[1]; b.remapflag[2] = remapflag[2]; b.remapflag[3] = remapflag[3];
    b.image = image;
    b.dummy = dummy;

    return b;
  }
};

template<typename To, typename From>
struct Transform {
  static constexpr bool is_identity = std::is_same_v<To, From>;
  KOKKOS_INLINE_FUNCTION
  static To transform(const From &x)
  {
    return x;
  }
};

// Simplified Transform specializations leveraging the new constructors
template<>
struct Transform<BodyKokkos, FixRigidSmall::Body> {
  static constexpr bool is_identity = false;
  KOKKOS_INLINE_FUNCTION
  static BodyKokkos transform(const FixRigidSmall::Body &b)
  {
    return BodyKokkos(b);
  }
};

template<>
struct Transform<FixRigidSmall::Body, BodyKokkos> {
  static constexpr bool is_identity = false;

  KOKKOS_INLINE_FUNCTION
  static FixRigidSmall::Body transform(const BodyKokkos &bk)
  {
    return static_cast<FixRigidSmall::Body>(bk);
  }
};

#endif // !LMP_RIGID_BODY_KOKKOS_H
