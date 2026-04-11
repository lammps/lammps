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

#ifdef FIX_CLASS
// clang-format off
FixStyle(rigid/small/kk,FixRigidSmallKokkos<LMPDeviceType>);
FixStyle(rigid/small/kk/device,FixRigidSmallKokkos<LMPDeviceType>);
FixStyle(rigid/small/kk/host,FixRigidSmallKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_RIGID_SMALL_KOKKOS_H
#define LMP_FIX_RIGID_SMALL_KOKKOS_H

#include "fix_rigid_small.h"
#include "kokkos_type.h"
#include "kokkos_few.h"

#include <Kokkos_SIMD.hpp>

namespace LAMMPS_NS {

  using KK_FLOAT_SIMD = Kokkos::Experimental::simd<KK_FLOAT>;

  struct BodyKokkos {
    int natoms;            // total number of atoms in body
    int ilocal;            // index of owning atom
    KK_FLOAT mass;           // total mass of body
    KK_FLOAT_SIMD xcm[3];         // COM position
    KK_FLOAT_SIMD xgc[3];         // geometric center position
    KK_FLOAT_SIMD vcm[3];         // COM velocity
    KK_FLOAT_SIMD fcm[3];         // force on COM
    KK_FLOAT_SIMD torque[3];      // torque around COM
    KK_FLOAT_SIMD quat[4];        // quaternion for orientation of body
    KK_FLOAT_SIMD inertia[3];     // 3 principal components of inertia
    //KK_FLOAT ex_space[3];    // principal axes in space coords
    //KK_FLOAT ey_space[3];
    //KK_FLOAT ez_space[3];
    //KK_FLOAT xgc_body[3];    // geometric center relative to xcm in body coords
    KK_FLOAT_SIMD angmom[3];      // space-frame angular momentum of body
    KK_FLOAT_SIMD omega[3];       // space-frame omega of body
    KK_FLOAT_SIMD conjqm[4];      // conjugate quaternion momentum
    int remapflag[4];      // PBC remap flags
    imageint image;        // image flags of xcm
    imageint dummy;        // dummy entry for better alignment
  };


// 1. The Generic Blueprint (MUST COME FIRST)
template<typename To, typename From>
struct Transformer {
  static constexpr bool is_identity = std::is_same_v<To, From>;

  KOKKOS_INLINE_FUNCTION
  static To apply(const From& x) {
    return x; // only valid if identical
  }
};

  // 1. Legacy -> Kokkos (Downcasting & Dropping Redundant Fields)
  template<>
  struct Transformer<BodyKokkos, Body> {
  static constexpr bool is_identity = false;

  KOKKOS_INLINE_FUNCTION
  static BodyKokkos apply(const Body& b) {
    BodyKokkos bk;
    bk.natoms = b.natoms;
    bk.ilocal = b.ilocal;
    bk.mass   = static_cast<KK_FLOAT>(b.mass);

    // 3-element arrays
    for (int i = 0; i < 3; i++) {
      bk.xcm[i]     = static_cast<KK_FLOAT_SIMD>(b.xcm[i]);
      bk.xgc[i]     = static_cast<KK_FLOAT_SIMD>(b.xgc[i]);
      bk.vcm[i]     = static_cast<KK_FLOAT_SIMD>(b.vcm[i]);
      bk.fcm[i]     = static_cast<KK_FLOAT_SIMD>(b.fcm[i]);
      bk.torque[i]  = static_cast<KK_FLOAT_SIMD>(b.torque[i]);
      bk.inertia[i] = static_cast<KK_FLOAT_SIMD>(b.inertia[i]);
      bk.angmom[i]  = static_cast<KK_FLOAT_SIMD>(b.angmom[i]);
      bk.omega[i]   = static_cast<KK_FLOAT_SIMD>(b.omega[i]);
    }

    // 4-element arrays
    for (int i = 0; i < 4; i++) {
      bk.quat[i]      = static_cast<KK_FLOAT_SIMD>(b.quat[i]);
      bk.conjqm[i]    = static_cast<KK_FLOAT_SIMD>(b.conjqm[i]);
      bk.remapflag[i] = b.remapflag[i];
    }

    bk.image = b.image;
    bk.dummy = b.dummy;

    // Note: ex_space, ey_space, ez_space, and xgc_body are intentionally dropped.
    return bk;
  }
};

  // 2. Kokkos -> Legacy (Upcasting & Reconstructing Fields)
  template<>
  struct Transformer<Body, BodyKokkos> {
  static constexpr bool is_identity = false;

  KOKKOS_INLINE_FUNCTION
  static Body apply(const BodyKokkos& bk) {
    Body b;
    b.natoms = bk.natoms;
    b.ilocal = bk.ilocal;
    b.mass   = static_cast<double>(bk.mass);

    for (int i = 0; i < 3; i++) {
      b.xcm[i]     = static_cast<double>(bk.xcm[i]);
      b.xgc[i]     = static_cast<double>(bk.xgc[i]);
      b.vcm[i]     = static_cast<double>(bk.vcm[i]);
      b.fcm[i]     = static_cast<double>(bk.fcm[i]);
      b.torque[i]  = static_cast<double>(bk.torque[i]);
      b.inertia[i] = static_cast<double>(bk.inertia[i]);
      b.angmom[i]  = static_cast<double>(bk.angmom[i]);
      b.omega[i]   = static_cast<double>(bk.omega[i]);
    }

    for (int i = 0; i < 4; i++) {
      b.quat[i]      = static_cast<double>(bk.quat[i]);
      b.conjqm[i]    = static_cast<double>(bk.conjqm[i]);
      b.remapflag[i] = bk.remapflag[i];
    }

    b.image = bk.image;
    b.dummy = bk.dummy;

    // ⚠️ CRITICAL: The Kokkos struct doesn't have the space vectors.
    // For now, we zero them out to prevent uninitialized memory bugs in legacy code.
    for(int i = 0; i < 3; i++) {
        b.ex_space[i] = 0.0;
        b.ey_space[i] = 0.0;
        b.ez_space[i] = 0.0;
        b.xgc_body[i] = 0.0;
    }

    return b;
  }
};


struct TagRigidSmallInitialIntegrate {};
struct TagRigidSmallFinalIntegrate {};

template<int TRICLINIC, int NEIGHFLAG, int EVFLAG>
struct TagRigidSmallSetXV {};

template<int TRICLINIC, int NEIGHFLAG, int EVFLAG>
struct TagRigidSmallSetV {};

struct TagRigidSmallComputeForcesTorquesZero {};
struct TagRigidSmallComputeForcesTorques {};
struct TagRigidSmallImageShift {};
struct TagRigidSmallEnforce2d {};

template<class DeviceType>
class FixRigidSmallKokkos : public FixRigidSmall {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef EV_FLOAT value_type;

  FixRigidSmallKokkos(class LAMMPS *, int, char **);
  ~FixRigidSmallKokkos() override;

  void init() override;
  void setup(int) override;
  void initial_integrate(int) override;
  void post_force(int) override;
  void final_integrate() override;
  void setup_pre_neighbor() override;
  void pre_neighbor() override;

  void grow_arrays(int) override;
  void copy_arrays(int, int, int) override;
  void set_arrays(int) override;
  void set_molecule(int, tagint, int, double *, double *, double *) override;

  int pack_exchange(int, double *) override;
  int unpack_exchange(int, double *) override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;

  void reset_dt() override;
  void zero_momentum() override;
  void zero_rotation() override;

  double compute_scalar() override;

  void grow_body() override;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagRigidSmallInitialIntegrate, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagRigidSmallFinalIntegrate, const int&) const;

  template<int TRICLINIC, int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagRigidSmallSetXV<TRICLINIC,NEIGHFLAG,EVFLAG>, const int&) const;

  template<int TRICLINIC, int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagRigidSmallSetXV<TRICLINIC,NEIGHFLAG,EVFLAG>, const int&, EV_FLOAT &) const;

  template<int TRICLINIC, int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagRigidSmallSetV<TRICLINIC,NEIGHFLAG,EVFLAG>, const int&) const;

  template<int TRICLINIC, int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagRigidSmallSetV<TRICLINIC,NEIGHFLAG,EVFLAG>, const int&, EV_FLOAT &) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagRigidSmallComputeForcesTorquesZero, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagRigidSmallComputeForcesTorques, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagRigidSmallImageShift, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagRigidSmallEnforce2d, const int&) const;

 protected:
  class AtomKokkos *atomKK;
  ExecutionSpace execution_space;

  TransformView<BodyKokkos*, Body*, Kokkos::LayoutRight, DeviceType> k_body;
  typename Kokkos::DualView<BodyKokkos*, Kokkos::LayoutRight, DeviceType>::t_dev d_body;

  DAT::tdual_int_1d k_bodyown;
  DAT::tdual_tagint_1d k_bodytag;
  DAT::tdual_int_1d k_atom2body;
  DAT::tdual_imageint_1d k_xcmimage;
  DAT::tdual_kkfloat_2d k_displace;
  DAT::tdual_int_1d k_eflags;

  typename AT::t_int_1d d_bodyown;
  typename AT::t_tagint_1d d_bodytag;
  typename AT::t_int_1d d_atom2body;
  typename AT::t_imageint_1d d_xcmimage;
  typename AT::t_kkfloat_2d d_displace;
  typename AT::t_int_1d d_eflags;

  using KKDeviceType = typename KKDevice<DeviceType>::value;

  template<typename DataType, typename Layout>
  using DupScatterView =
      KKScatterView<DataType, Layout, KKDeviceType, KKScatterSum, KKScatterDuplicated>;

  template<typename DataType, typename Layout>
  using NonDupScatterView =
      KKScatterView<DataType, Layout, KKDeviceType, KKScatterSum, KKScatterNonDuplicated>;

  DupScatterView<KK_ACC_FLOAT *[6], typename AT::t_kkacc_1d_6::array_layout> dup_vatom;
  NonDupScatterView<KK_ACC_FLOAT *[6], typename AT::t_kkacc_1d_6::array_layout> ndup_vatom;

  typename AT::t_kkfloat_1d_3_lr d_x;
  typename AT::t_kkfloat_1d_3 d_v;
  typename AT::t_kkacc_1d_3 d_f;
  typename AT::t_kkfloat_1d d_rmass;
  typename AT::t_kkfloat_1d d_mass;
  typename AT::t_int_1d d_type;
  typename AT::t_int_1d d_mask;
  typename AT::t_imageint_1d d_image;

  DAT::ttransform_kkacc_1d k_eatom;
  DAT::ttransform_kkacc_1d_6 k_vatom;
  typename AT::t_kkacc_1d d_eatom;
  typename AT::t_kkacc_1d_6 d_vatom;

  Few<double,3> d_prd;
  Few<double,6> d_h;

  void sync_body_device();
  void sync_body_host();
  void sync_fix_data_device();

  void compute_forces_and_torques_kokkos();
  void enforce2d_kokkos();
  void image_shift_kokkos();
  void grow_body_kokkos();

  template<int TRICLINIC, int EVFLAG>
  void set_xv_kokkos();

  template<int TRICLINIC, int EVFLAG>
  void set_v_kokkos();
};

}    // namespace LAMMPS_NS

#endif // !LMP_FIX_RIGID_SMALL_KOKKOS_H
#endif
