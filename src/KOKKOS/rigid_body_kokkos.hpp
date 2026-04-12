/* -*- c++ -*- ----------------------------------------------------------
   Shared Kokkos rigid body storage and legacy transforms for
   fix_rigid_small_kokkos and fix_rigid_nh_small_kokkos.
------------------------------------------------------------------------- */

#ifndef LMP_RIGID_BODY_KOKKOS_HPP
#define LMP_RIGID_BODY_KOKKOS_HPP

#include "fix_rigid_small.h"
#include "kokkos_type.h"

#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>
#include <type_traits>

namespace LAMMPS_NS {

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

template<>
struct Transform<BodyKokkos, FixRigidSmall::Body> {
  static constexpr bool is_identity = false;
  KOKKOS_INLINE_FUNCTION
  static BodyKokkos transform(const FixRigidSmall::Body &b)
  {
    return BodyKokkos{b.natoms,
                      b.ilocal,
                      static_cast<KK_FLOAT>(b.mass),
                      {static_cast<KK_FLOAT>(b.xcm[0]), static_cast<KK_FLOAT>(b.xcm[1]),
                       static_cast<KK_FLOAT>(b.xcm[2])},
                      {static_cast<KK_FLOAT>(b.xgc[0]), static_cast<KK_FLOAT>(b.xgc[1]),
                       static_cast<KK_FLOAT>(b.xgc[2])},
                      {static_cast<KK_FLOAT>(b.vcm[0]), static_cast<KK_FLOAT>(b.vcm[1]),
                       static_cast<KK_FLOAT>(b.vcm[2])},
                      {static_cast<KK_FLOAT>(b.fcm[0]), static_cast<KK_FLOAT>(b.fcm[1]),
                       static_cast<KK_FLOAT>(b.fcm[2])},
                      {static_cast<KK_FLOAT>(b.torque[0]), static_cast<KK_FLOAT>(b.torque[1]),
                       static_cast<KK_FLOAT>(b.torque[2])},
                      {static_cast<KK_FLOAT>(b.quat[0]), static_cast<KK_FLOAT>(b.quat[1]),
                       static_cast<KK_FLOAT>(b.quat[2]), static_cast<KK_FLOAT>(b.quat[3])},
                      {static_cast<KK_FLOAT>(b.inertia[0]), static_cast<KK_FLOAT>(b.inertia[1]),
                       static_cast<KK_FLOAT>(b.inertia[2])},
                      {static_cast<KK_FLOAT>(b.ex_space[0]), static_cast<KK_FLOAT>(b.ex_space[1]),
                       static_cast<KK_FLOAT>(b.ex_space[2])},
                      {static_cast<KK_FLOAT>(b.ey_space[0]), static_cast<KK_FLOAT>(b.ey_space[1]),
                       static_cast<KK_FLOAT>(b.ey_space[2])},
                      {static_cast<KK_FLOAT>(b.ez_space[0]), static_cast<KK_FLOAT>(b.ez_space[1]),
                       static_cast<KK_FLOAT>(b.ez_space[2])},
                      {static_cast<KK_FLOAT>(b.xgc_body[0]), static_cast<KK_FLOAT>(b.xgc_body[1]),
                       static_cast<KK_FLOAT>(b.xgc_body[2])},
                      {static_cast<KK_FLOAT>(b.angmom[0]), static_cast<KK_FLOAT>(b.angmom[1]),
                       static_cast<KK_FLOAT>(b.angmom[2])},
                      {static_cast<KK_FLOAT>(b.omega[0]), static_cast<KK_FLOAT>(b.omega[1]),
                       static_cast<KK_FLOAT>(b.omega[2])},
                      {static_cast<KK_FLOAT>(b.conjqm[0]), static_cast<KK_FLOAT>(b.conjqm[1]),
                       static_cast<KK_FLOAT>(b.conjqm[2]), static_cast<KK_FLOAT>(b.conjqm[3])},
                      {b.remapflag[0], b.remapflag[1], b.remapflag[2], b.remapflag[3]},
                      b.image,
                      b.dummy};
  }
};

template<>
struct Transform<FixRigidSmall::Body, BodyKokkos> {
  static constexpr bool is_identity = false;

  KOKKOS_INLINE_FUNCTION
  static FixRigidSmall::Body transform(const BodyKokkos &bk)
  {
    return FixRigidSmall::Body{
        bk.natoms,
        bk.ilocal,
        static_cast<double>(bk.mass),
        {static_cast<double>(bk.xcm[0]), static_cast<double>(bk.xcm[1]),
         static_cast<double>(bk.xcm[2])},
        {static_cast<double>(bk.xgc[0]), static_cast<double>(bk.xgc[1]),
         static_cast<double>(bk.xgc[2])},
        {static_cast<double>(bk.vcm[0]), static_cast<double>(bk.vcm[1]),
         static_cast<double>(bk.vcm[2])},
        {static_cast<double>(bk.fcm[0]), static_cast<double>(bk.fcm[1]),
         static_cast<double>(bk.fcm[2])},
        {static_cast<double>(bk.torque[0]), static_cast<double>(bk.torque[1]),
         static_cast<double>(bk.torque[2])},
        {static_cast<double>(bk.quat[0]), static_cast<double>(bk.quat[1]),
         static_cast<double>(bk.quat[2]), static_cast<double>(bk.quat[3])},
        {static_cast<double>(bk.inertia[0]), static_cast<double>(bk.inertia[1]),
         static_cast<double>(bk.inertia[2])},
        {static_cast<double>(bk.ex_space[0]), static_cast<double>(bk.ex_space[1]),
         static_cast<double>(bk.ex_space[2])},
        {static_cast<double>(bk.ey_space[0]), static_cast<double>(bk.ey_space[1]),
         static_cast<double>(bk.ey_space[2])},
        {static_cast<double>(bk.ez_space[0]), static_cast<double>(bk.ez_space[1]),
         static_cast<double>(bk.ez_space[2])},
        {static_cast<double>(bk.xgc_body[0]), static_cast<double>(bk.xgc_body[1]),
         static_cast<double>(bk.xgc_body[2])},
        {static_cast<double>(bk.angmom[0]), static_cast<double>(bk.angmom[1]),
         static_cast<double>(bk.angmom[2])},
        {static_cast<double>(bk.omega[0]), static_cast<double>(bk.omega[1]),
         static_cast<double>(bk.omega[2])},
        {static_cast<double>(bk.conjqm[0]), static_cast<double>(bk.conjqm[1]),
         static_cast<double>(bk.conjqm[2]), static_cast<double>(bk.conjqm[3])},
        {bk.remapflag[0], bk.remapflag[1], bk.remapflag[2], bk.remapflag[3]},
        bk.image,
        bk.dummy};
  }
};

template<class DeviceType>
void rigid_body_copy_legacy_from_kk_host(
    Kokkos::DualView<BodyKokkos *, Kokkos::LayoutRight, DeviceType> &k_body,
    FixRigidSmall::Body *body, int nmax_body)
{
  auto hv = k_body.view_host();
  for (int i = 0; i < nmax_body; i++)
    body[i] = Transform<FixRigidSmall::Body, BodyKokkos>::transform(hv(i));
}

template<class DeviceType>
void rigid_body_copy_kk_from_legacy_host(
    Kokkos::DualView<BodyKokkos *, Kokkos::LayoutRight, DeviceType> &k_body,
    FixRigidSmall::Body *body, int nmax_body)
{
  auto hv = k_body.view_host();
  for (int i = 0; i < nmax_body; i++)
    hv(i) = Transform<BodyKokkos, FixRigidSmall::Body>::transform(body[i]);
  k_body.modify_host();
}

}    // namespace LAMMPS_NS

#endif
