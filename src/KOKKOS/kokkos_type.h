/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_LMPTYPE_KOKKOS_H
#define LMP_LMPTYPE_KOKKOS_H

#include "pointers.h"

#include <Kokkos_Core.hpp>
#include <Kokkos_TransformDualView.hpp>
#include <Kokkos_DualView.hpp>
#include <impl/Kokkos_Timer.hpp>
#include <Kokkos_Vectorization.hpp>
#include <Kokkos_ScatterView.hpp>

enum{FULL=1u,HALFTHREAD=2u,HALF=4u};

#if defined(KOKKOS_ENABLE_CXX11)
#undef ISFINITE
#define ISFINITE(x) std::isfinite(x)
#endif


#if defined(LMP_KOKKOS_SINGLE)
typedef float KK_FLOAT;
#define MPI_KK_FLOAT MPI_FLOAT
#else
typedef double KK_FLOAT;
#define MPI_KK_FLOAT MPI_DOUBLE
#endif

#define MAX_TYPES_STACKPARAMS 12
#define NeighClusterSize 8

template<class Scalar>
struct t_scalar3 {
  Scalar x,y,z;

  KOKKOS_FORCEINLINE_FUNCTION
  t_scalar3() {
    x = 0; y = 0; z = 0;
  }

  KOKKOS_FORCEINLINE_FUNCTION
  t_scalar3(const t_scalar3& rhs) {
    x = rhs.x; y = rhs.y; z = rhs.z;
  }

  KOKKOS_FORCEINLINE_FUNCTION
  t_scalar3(const Scalar& x_, const Scalar& y_, const Scalar& z_ ) {
    x = x_; y = y_; z = z_;
  }

  KOKKOS_FORCEINLINE_FUNCTION
  t_scalar3 operator= (const t_scalar3& rhs) {
    x = rhs.x; y = rhs.y; z = rhs.z;
    return *this;
  }

  KOKKOS_FORCEINLINE_FUNCTION
  t_scalar3 operator= (const volatile t_scalar3& rhs) {
    x = rhs.x; y = rhs.y; z = rhs.z;
    return *this;
  }

  KOKKOS_FORCEINLINE_FUNCTION
  t_scalar3 operator+= (const t_scalar3& rhs) {
    x += rhs.x; y += rhs.y; z += rhs.z;
    return *this;
  }

  KOKKOS_FORCEINLINE_FUNCTION
  t_scalar3 operator+= (const volatile t_scalar3& rhs) volatile {
    x += rhs.x; y += rhs.y; z += rhs.z;
    return *this;
  }
};

template<class Scalar>
KOKKOS_FORCEINLINE_FUNCTION
t_scalar3<Scalar> operator +
  (const t_scalar3<Scalar>& a, const t_scalar3<Scalar>& b) {
  return t_scalar3<Scalar>(a.x+b.x,a.y+b.y,a.z+b.z);
}

template<class Scalar>
KOKKOS_FORCEINLINE_FUNCTION
t_scalar3<Scalar> operator *
  (const t_scalar3<Scalar>& a, const Scalar& b) {
  return t_scalar3<Scalar>(a.x*b,a.y*b,a.z*b);
}

template<class Scalar>
KOKKOS_FORCEINLINE_FUNCTION
t_scalar3<Scalar> operator *
  (const Scalar& b, const t_scalar3<Scalar>& a) {
  return t_scalar3<Scalar>(a.x*b,a.y*b,a.z*b);
}

// set LMPHostype and Device from Kokkos Default Types
typedef Kokkos::DefaultExecutionSpace LMPDeviceType;
typedef Kokkos::HostSpace::execution_space LMPHostType;

#if defined(LMP_KOKKOS_LAYOUT_LEFT)
typedef LMPDeviceType::array_layout device_array_layout;
#else
typedef Kokkos::LayoutRight device_array_layout;
#endif

// Get DeviceType from ExecutionSpace

template<LAMMPS_NS::ExecutionSpace Space>
struct GetDeviceType;

template<>
struct GetDeviceType<LAMMPS_NS::Host> {
  typedef LMPHostType value;
};

template<>
struct GetDeviceType<LAMMPS_NS::Device> {
  typedef LMPDeviceType value;
};

// Get float type from ExecutionSpace

template<LAMMPS_NS::ExecutionSpace Space>
struct GetFloatType;

template<>
struct GetFloatType<LAMMPS_NS::Host> {
  typedef double type;
  static MPI_Datatype GetMPIType() {return MPI_DOUBLE;}
};

template<>
struct GetFloatType<LAMMPS_NS::Device> {
  typedef KK_FLOAT type;
  static MPI_Datatype GetMPIType() {return MPI_KK_FLOAT;}
};

// Get View from DualView depending on ExecutionSpace

template<LAMMPS_NS::ExecutionSpace Space>
class DualViewHelper {};

template<>
class DualViewHelper<LAMMPS_NS::Host> {
public:
  template <typename TYPE>
  KOKKOS_INLINE_FUNCTION
  static typename TYPE::t_host view(const TYPE &dualview) {
    return dualview.h_view;
  }

  template <typename TYPE>
  static void sync(TYPE &dualview) {
    dualview.sync_host();
  }

  template <typename TYPE>
  static void modify(TYPE &dualview) {
    dualview.modify_host();
  }
};

template<>
class DualViewHelper<LAMMPS_NS::Device> {
public:
template<typename TYPE>
  KOKKOS_INLINE_FUNCTION
  static typename TYPE::t_dev view(const TYPE dualview) {
    return dualview.d_view;
  }

  template <typename TYPE>
  static void sync(TYPE &dualview) {
    dualview.sync_device();
  }

  template <typename TYPE>
  static void modify(TYPE &dualview) {
    dualview.modify_device();
  }
};

// Need to use Cuda UVM memory space for Host execution space

template<class DeviceType>
class KKDevice {
public:
#if defined(KOKKOS_ENABLE_CUDA) && defined(KOKKOS_ENABLE_CUDA_UVM)
  typedef Kokkos::Device<DeviceType,LMPDeviceType::memory_space> value;
#else
  typedef Kokkos::Device<DeviceType, typename DeviceType::memory_space> value;
#endif
};


// Determine memory traits for force array
// Do atomic trait when running HALFTHREAD neighbor list style
template<int NEIGHFLAG>
struct AtomicF {
  enum {value = Kokkos::Unmanaged};
};

template<>
struct AtomicF<HALFTHREAD> {
  enum {value = Kokkos::Atomic|Kokkos::Unmanaged};
};


// Determine memory traits for force array
// Do atomic trait when running HALFTHREAD neighbor list style with CUDA
template<int NEIGHFLAG, class DeviceType>
struct AtomicDup {
  enum {value = Kokkos::Experimental::ScatterNonAtomic};
};

#ifdef KOKKOS_ENABLE_CUDA
template<>
struct AtomicDup<HALFTHREAD,Kokkos::Cuda> {
  enum {value = Kokkos::Experimental::ScatterAtomic};
};
#endif

#ifdef LMP_KOKKOS_USE_ATOMICS

#ifdef KOKKOS_ENABLE_OPENMP
template<>
struct AtomicDup<HALFTHREAD,Kokkos::OpenMP> {
  enum {value = Kokkos::Experimental::ScatterAtomic};
};
#endif

#ifdef KOKKOS_ENABLE_THREADS
template<>
struct AtomicDup<HALFTHREAD,Kokkos::Threads> {
  enum {value = Kokkos::Experimental::ScatterAtomic};
};
#endif

#endif


// Determine duplication traits for force array
// Use duplication when running threaded and not using atomics
template<int NEIGHFLAG, class DeviceType>
struct NeedDup {
  enum {value = Kokkos::Experimental::ScatterNonDuplicated};
};

#ifndef LMP_KOKKOS_USE_ATOMICS

#ifdef KOKKOS_ENABLE_OPENMP
template<>
struct NeedDup<HALFTHREAD,Kokkos::OpenMP> {
  enum {value = Kokkos::Experimental::ScatterDuplicated};
};
#endif

#ifdef KOKKOS_ENABLE_THREADS
template<>
struct NeedDup<HALFTHREAD,Kokkos::Threads> {
  enum {value = Kokkos::Experimental::ScatterDuplicated};
};
#endif

#endif

template<int value, typename T1, typename T2>
class ScatterViewHelper {};

template<typename T1, typename T2>
class ScatterViewHelper<Kokkos::Experimental::ScatterDuplicated,T1,T2> {
public:
  KOKKOS_INLINE_FUNCTION
  static T1 get(const T1 &dup, const T2 &nondup) {
    return dup;
  }
};

template<typename T1, typename T2>
class ScatterViewHelper<Kokkos::Experimental::ScatterNonDuplicated,T1,T2> {
public:
  KOKKOS_INLINE_FUNCTION
  static T2 get(const T1 &dup, const T2 &nondup) {
    return nondup;
  }
};

struct s_EV_FLOAT {
  KK_FLOAT evdwl;
  KK_FLOAT ecoul;
  KK_FLOAT v[6];
  KOKKOS_INLINE_FUNCTION
  s_EV_FLOAT() {
    evdwl = 0;
    ecoul = 0;
    v[0] = 0; v[1] = 0; v[2] = 0;
    v[3] = 0; v[4] = 0; v[5] = 0;
  }

  KOKKOS_INLINE_FUNCTION
  void operator+=(const s_EV_FLOAT &rhs) {
    evdwl += rhs.evdwl;
    ecoul += rhs.ecoul;
    v[0] += rhs.v[0];
    v[1] += rhs.v[1];
    v[2] += rhs.v[2];
    v[3] += rhs.v[3];
    v[4] += rhs.v[4];
    v[5] += rhs.v[5];
  }

  KOKKOS_INLINE_FUNCTION
  void operator+=(const volatile s_EV_FLOAT &rhs) volatile {
    evdwl += rhs.evdwl;
    ecoul += rhs.ecoul;
    v[0] += rhs.v[0];
    v[1] += rhs.v[1];
    v[2] += rhs.v[2];
    v[3] += rhs.v[3];
    v[4] += rhs.v[4];
    v[5] += rhs.v[5];
  }
};
typedef struct s_EV_FLOAT EV_FLOAT;

struct s_EV_FLOAT_REAX {
  KK_FLOAT evdwl;
  KK_FLOAT ecoul;
  KK_FLOAT v[6];
  KK_FLOAT ereax[10];
  KOKKOS_INLINE_FUNCTION
  s_EV_FLOAT_REAX() {
    evdwl = 0;
    ecoul = 0;
    v[0] = 0; v[1] = 0; v[2] = 0;
    v[3] = 0; v[4] = 0; v[5] = 0;
    ereax[0] = 0; ereax[1] = 0; ereax[2] = 0;
    ereax[3] = 0; ereax[4] = 0; ereax[5] = 0;
    ereax[6] = 0; ereax[7] = 0; ereax[8] = 0;
  }

  KOKKOS_INLINE_FUNCTION
  void operator+=(const s_EV_FLOAT_REAX &rhs) {
    evdwl += rhs.evdwl;
    ecoul += rhs.ecoul;
    v[0] += rhs.v[0];
    v[1] += rhs.v[1];
    v[2] += rhs.v[2];
    v[3] += rhs.v[3];
    v[4] += rhs.v[4];
    v[5] += rhs.v[5];
    ereax[0] += rhs.ereax[0];
    ereax[1] += rhs.ereax[1];
    ereax[2] += rhs.ereax[2];
    ereax[3] += rhs.ereax[3];
    ereax[4] += rhs.ereax[4];
    ereax[5] += rhs.ereax[5];
    ereax[6] += rhs.ereax[6];
    ereax[7] += rhs.ereax[7];
    ereax[8] += rhs.ereax[8];
  }

  KOKKOS_INLINE_FUNCTION
  void operator+=(const volatile s_EV_FLOAT_REAX &rhs) volatile {
    evdwl += rhs.evdwl;
    ecoul += rhs.ecoul;
    v[0] += rhs.v[0];
    v[1] += rhs.v[1];
    v[2] += rhs.v[2];
    v[3] += rhs.v[3];
    v[4] += rhs.v[4];
    v[5] += rhs.v[5];
    ereax[0] += rhs.ereax[0];
    ereax[1] += rhs.ereax[1];
    ereax[2] += rhs.ereax[2];
    ereax[3] += rhs.ereax[3];
    ereax[4] += rhs.ereax[4];
    ereax[5] += rhs.ereax[5];
    ereax[6] += rhs.ereax[6];
    ereax[7] += rhs.ereax[7];
    ereax[8] += rhs.ereax[8];
  }
};
typedef struct s_EV_FLOAT_REAX EV_FLOAT_REAX;

struct s_FEV_FLOAT {
  KK_FLOAT f[3];
  KK_FLOAT evdwl;
  KK_FLOAT ecoul;
  KK_FLOAT v[6];
  KOKKOS_INLINE_FUNCTION
  s_FEV_FLOAT() {
    f[0] = 0; f[1] = 0; f[2] = 0;
    evdwl = 0;
    ecoul = 0;
    v[0] = 0; v[1] = 0; v[2] = 0;
    v[3] = 0; v[4] = 0; v[5] = 0;
  }

  KOKKOS_INLINE_FUNCTION
  void operator+=(const s_FEV_FLOAT &rhs) {
    f[0] += rhs.f[0];
    f[1] += rhs.f[1];
    f[2] += rhs.f[2];
    evdwl += rhs.evdwl;
    ecoul += rhs.ecoul;
    v[0] += rhs.v[0];
    v[1] += rhs.v[1];
    v[2] += rhs.v[2];
    v[3] += rhs.v[3];
    v[4] += rhs.v[4];
    v[5] += rhs.v[5];
  }

  KOKKOS_INLINE_FUNCTION
  void operator+=(const volatile s_FEV_FLOAT &rhs) volatile {
    f[0] += rhs.f[0];
    f[1] += rhs.f[1];
    f[2] += rhs.f[2];
    evdwl += rhs.evdwl;
    ecoul += rhs.ecoul;
    v[0] += rhs.v[0];
    v[1] += rhs.v[1];
    v[2] += rhs.v[2];
    v[3] += rhs.v[3];
    v[4] += rhs.v[4];
    v[5] += rhs.v[5];
  }
};
typedef struct s_FEV_FLOAT FEV_FLOAT;

// ------------------------------------------------------------------------

// LAMMPS types

template <LAMMPS_NS::ExecutionSpace Space>
struct ArrayTypes;

template <>
struct ArrayTypes<LAMMPS_NS::Device> {

typedef Kokkos::
  DualView<int, LMPDeviceType::array_layout, LMPDeviceType> tdual_int_scalar;
typedef tdual_int_scalar::t_dev t_int_scalar;
typedef tdual_int_scalar::t_dev_const t_int_scalar_const;
typedef tdual_int_scalar::t_dev_um t_int_scalar_um;
typedef tdual_int_scalar::t_dev_const_um t_int_scalar_const_um;

typedef Kokkos::
  DualView<int*, LMPDeviceType::array_layout, LMPDeviceType> tdual_int_1d;
typedef tdual_int_1d::t_dev t_int_1d;
typedef tdual_int_1d::t_dev_const t_int_1d_const;
typedef tdual_int_1d::t_dev_um t_int_1d_um;
typedef tdual_int_1d::t_dev_const_um t_int_1d_const_um;
typedef tdual_int_1d::t_dev_const_randomread t_int_1d_randomread;

typedef Kokkos::
  TransformDualView<int*[3], device_array_layout, LMPDeviceType, void, int*[3], Kokkos::LayoutRight> tdual_int_1d_3;
typedef tdual_int_1d_3::t_dev t_int_1d_3;
typedef tdual_int_1d_3::t_dev_const t_int_1d_3_const;
typedef tdual_int_1d_3::t_dev_um t_int_1d_3_um;
typedef tdual_int_1d_3::t_dev_const_um t_int_1d_3_const_um;
typedef tdual_int_1d_3::t_dev_const_randomread t_int_1d_3_randomread;

typedef Kokkos::
  TransformDualView<int**, device_array_layout, LMPDeviceType, void, int**, Kokkos::LayoutRight> tdual_int_2d;
typedef tdual_int_2d::t_dev t_int_2d;
typedef tdual_int_2d::t_dev_const t_int_2d_const;
typedef tdual_int_2d::t_dev_um t_int_2d_um;
typedef tdual_int_2d::t_dev_const_um t_int_2d_const_um;
typedef tdual_int_2d::t_dev_const_randomread t_int_2d_randomread;

typedef Kokkos::
  DualView<int**, LMPDeviceType::array_layout, LMPDeviceType> tdual_int_2d_dl;
typedef tdual_int_2d_dl::t_dev t_int_2d_dl;
typedef tdual_int_2d_dl::t_dev_const t_int_2d_dl_const;
typedef tdual_int_2d_dl::t_dev_um t_int_2d_dl_um;
typedef tdual_int_2d_dl::t_dev_const_um t_int_2d_dl_const_um;
typedef tdual_int_2d_dl::t_dev_const_randomread t_int_2d_dl_randomread;

typedef Kokkos::
  DualView<int**, Kokkos::LayoutRight, LMPDeviceType> tdual_int_2d_lr;
typedef tdual_int_2d_lr::t_dev t_int_2d_lr;
typedef tdual_int_2d_lr::t_dev_const t_int_2d_lr_const;
typedef tdual_int_2d_lr::t_dev_um t_int_2d_lr_um;
typedef tdual_int_2d_lr::t_dev_const_um t_int_2d_lr_const_um;
typedef tdual_int_2d_lr::t_dev_const_randomread t_int_2d_lr_randomread;

typedef Kokkos::
  DualView<LAMMPS_NS::tagint*, LMPDeviceType::array_layout, LMPDeviceType>
  tdual_tagint_1d;
typedef tdual_tagint_1d::t_dev t_tagint_1d;
typedef tdual_tagint_1d::t_dev_const t_tagint_1d_const;
typedef tdual_tagint_1d::t_dev_um t_tagint_1d_um;
typedef tdual_tagint_1d::t_dev_const_um t_tagint_1d_const_um;
typedef tdual_tagint_1d::t_dev_const_randomread t_tagint_1d_randomread;

typedef Kokkos::
  TransformDualView<LAMMPS_NS::tagint**, device_array_layout, LMPDeviceType, void, LAMMPS_NS::tagint**, Kokkos::LayoutRight>
  tdual_tagint_2d;
typedef tdual_tagint_2d::t_dev t_tagint_2d;
typedef tdual_tagint_2d::t_dev_const t_tagint_2d_const;
typedef tdual_tagint_2d::t_dev_um t_tagint_2d_um;
typedef tdual_tagint_2d::t_dev_const_um t_tagint_2d_const_um;
typedef tdual_tagint_2d::t_dev_const_randomread t_tagint_2d_randomread;

typedef Kokkos::
  DualView<LAMMPS_NS::imageint*, LMPDeviceType::array_layout, LMPDeviceType>
  tdual_imageint_1d;
typedef tdual_imageint_1d::t_dev t_imageint_1d;
typedef tdual_imageint_1d::t_dev_const t_imageint_1d_const;
typedef tdual_imageint_1d::t_dev_um t_imageint_1d_um;
typedef tdual_imageint_1d::t_dev_const_um t_imageint_1d_const_um;
typedef tdual_imageint_1d::t_dev_const_randomread t_imageint_1d_randomread;

typedef Kokkos::
  TransformDualView<KK_FLOAT, LMPDeviceType::array_layout, LMPDeviceType, void, double, LMPDeviceType::array_layout>
  tdual_float_scalar;
typedef tdual_float_scalar::t_dev t_float_scalar;
typedef tdual_float_scalar::t_dev_const t_float_scalar_const;
typedef tdual_float_scalar::t_dev_um t_float_scalar_um;
typedef tdual_float_scalar::t_dev_const_um t_float_scalar_const_um;

typedef Kokkos::TransformDualView<KK_FLOAT*, LMPDeviceType::array_layout, LMPDeviceType, void, double*, LMPDeviceType::array_layout> tdual_float_1d;
typedef tdual_float_1d::t_dev t_float_1d;
typedef tdual_float_1d::t_dev_const t_float_1d_const;
typedef tdual_float_1d::t_dev_um t_float_1d_um;
typedef tdual_float_1d::t_dev_const_um t_float_1d_const_um;
typedef tdual_float_1d::t_dev_const_randomread t_float_1d_randomread;

typedef Kokkos::TransformDualView<KK_FLOAT**, device_array_layout, LMPDeviceType, void, double**, Kokkos::LayoutRight> tdual_float_2d;
typedef tdual_float_2d::t_dev t_float_2d;
typedef tdual_float_2d::t_dev_const t_float_2d_const;
typedef tdual_float_2d::t_dev_um t_float_2d_um;
typedef tdual_float_2d::t_dev_const_um t_float_2d_const_um;
typedef tdual_float_2d::t_dev_const_randomread t_float_2d_randomread;

typedef Kokkos::TransformDualView<KK_FLOAT**, LMPDeviceType::array_layout, LMPDeviceType, void, double**, LMPDeviceType::array_layout> tdual_float_2d_dl;
typedef tdual_float_2d_dl::t_dev t_float_2d_dl;
typedef tdual_float_2d_dl::t_dev_const t_float_2d_dl_const;
typedef tdual_float_2d_dl::t_dev_um t_float_2d_dl_um;
typedef tdual_float_2d_dl::t_dev_const_um t_float_2d_dl_const_um;
typedef tdual_float_2d_dl::t_dev_const_randomread t_float_2d_dl_randomread;

typedef Kokkos::TransformDualView<KK_FLOAT**, Kokkos::LayoutRight, LMPDeviceType, void, double**, Kokkos::LayoutRight> tdual_float_2d_lr;
typedef tdual_float_2d_lr::t_dev t_float_2d_lr;
typedef tdual_float_2d_lr::t_dev_const t_float_2d_lr_const;
typedef tdual_float_2d_lr::t_dev_um t_float_2d_lr_um;
typedef tdual_float_2d_lr::t_dev_const_um t_float_2d_lr_const_um;
typedef tdual_float_2d_lr::t_dev_const_randomread t_float_2d_lr_randomread;

typedef Kokkos::DualView<double**, Kokkos::LayoutRight, LMPDeviceType> tdual_double_2d_lr;
typedef tdual_double_2d_lr::t_dev t_double_2d_lr;
typedef tdual_double_2d_lr::t_dev_const t_double_2d_lr_const;
typedef tdual_double_2d_lr::t_dev_um t_double_2d_lr_um;
typedef tdual_double_2d_lr::t_dev_const_um t_double_2d_lr_const_um;
typedef tdual_double_2d_lr::t_dev_const_randomread t_double_2d_lr_randomread;

typedef Kokkos::TransformDualView<KK_FLOAT*[3], device_array_layout, LMPDeviceType, void, double*[3], Kokkos::LayoutRight> tdual_float_1d_3;
typedef tdual_float_1d_3::t_dev t_float_1d_3;
typedef tdual_float_1d_3::t_dev_const t_float_1d_3_const;
typedef tdual_float_1d_3::t_dev_um t_float_1d_3_um;
typedef tdual_float_1d_3::t_dev_const_um t_float_1d_3_const_um;
typedef tdual_float_1d_3::t_dev_const_randomread t_float_1d_3_randomread;

typedef Kokkos::TransformDualView<KK_FLOAT*[6], device_array_layout, LMPDeviceType, void, double*[6], Kokkos::LayoutRight> tdual_float_1d_6;
typedef tdual_float_1d_6::t_dev t_float_1d_6;
typedef tdual_float_1d_6::t_dev_const t_float_1d_6_const;
typedef tdual_float_1d_6::t_dev_um t_float_1d_6_um;
typedef tdual_float_1d_6::t_dev_const_um t_float_1d_6_const_um;
typedef tdual_float_1d_6::t_dev_const_randomread t_float_1d_6_randomread;

typedef Kokkos::DualView<KK_FLOAT**[7], Kokkos::LayoutRight, LMPDeviceType> tdual_float_2d_7;
typedef tdual_float_2d_7::t_dev_const_randomread t_float_2d_7_randomread;
typedef tdual_float_2d_7::t_dev t_float_2d_7;

//Neighbor Types

typedef Kokkos::DualView<int**, LMPDeviceType::array_layout, LMPDeviceType> tdual_neighbors_2d;
typedef tdual_neighbors_2d::t_dev t_neighbors_2d;
typedef tdual_neighbors_2d::t_dev_const t_neighbors_2d_const;
typedef tdual_neighbors_2d::t_dev_um t_neighbors_2d_um;
typedef tdual_neighbors_2d::t_dev_const_um t_neighbors_2d_const_um;
typedef tdual_neighbors_2d::t_dev_const_randomread t_neighbors_2d_randomread;

};

template <>
struct ArrayTypes<LAMMPS_NS::Host> {

//Scalar Types

typedef Kokkos::
  DualView<int, LMPDeviceType::array_layout, LMPDeviceType> tdual_int_scalar;
typedef tdual_int_scalar::t_host t_int_scalar;
typedef tdual_int_scalar::t_host_const t_int_scalar_const;
typedef tdual_int_scalar::t_host_um t_int_scalar_um;
typedef tdual_int_scalar::t_host_const_um t_int_scalar_const_um;

typedef Kokkos::
  DualView<int*, LMPDeviceType::array_layout, LMPDeviceType> tdual_int_1d;
typedef tdual_int_1d::t_host t_int_1d;
typedef tdual_int_1d::t_host_const t_int_1d_const;
typedef tdual_int_1d::t_host_um t_int_1d_um;
typedef tdual_int_1d::t_host_const_um t_int_1d_const_um;
typedef tdual_int_1d::t_host_const_randomread t_int_1d_randomread;

typedef Kokkos::
  TransformDualView<int*[3], device_array_layout, LMPDeviceType, void, int*[3], Kokkos::LayoutRight> tdual_int_1d_3;
typedef tdual_int_1d_3::t_host t_int_1d_3;
typedef tdual_int_1d_3::t_host_const t_int_1d_3_const;
typedef tdual_int_1d_3::t_host_um t_int_1d_3_um;
typedef tdual_int_1d_3::t_host_const_um t_int_1d_3_const_um;
typedef tdual_int_1d_3::t_host_const_randomread t_int_1d_3_randomread;

typedef Kokkos::
  TransformDualView<int**, device_array_layout, LMPDeviceType, void, int**, Kokkos::LayoutRight> tdual_int_2d;
typedef tdual_int_2d::t_host t_int_2d;
typedef tdual_int_2d::t_host_const t_int_2d_const;
typedef tdual_int_2d::t_host_um t_int_2d_um;
typedef tdual_int_2d::t_host_const_um t_int_2d_const_um;
typedef tdual_int_2d::t_host_const_randomread t_int_2d_randomread;

typedef Kokkos::
  DualView<int**, LMPDeviceType::array_layout, LMPDeviceType> tdual_int_2d_dl;
typedef tdual_int_2d_dl::t_host t_int_2d_dl;
typedef tdual_int_2d_dl::t_host_const t_int_2d_dl_const;
typedef tdual_int_2d_dl::t_host_um t_int_2d_dl_um;
typedef tdual_int_2d_dl::t_host_const_um t_int_2d_dl_const_um;
typedef tdual_int_2d_dl::t_host_const_randomread t_int_2d_dl_randomread;

typedef Kokkos::
  DualView<int**, Kokkos::LayoutRight, LMPDeviceType> tdual_int_2d_lr;
typedef tdual_int_2d_lr::t_host t_int_2d_lr;
typedef tdual_int_2d_lr::t_host_const t_int_2d_lr_const;
typedef tdual_int_2d_lr::t_host_um t_int_2d_lr_um;
typedef tdual_int_2d_lr::t_host_const_um t_int_2d_lr_const_um;
typedef tdual_int_2d_lr::t_host_const_randomread t_int_2d_lr_randomread;

typedef Kokkos::
  DualView<LAMMPS_NS::tagint*, LMPDeviceType::array_layout, LMPDeviceType>
  tdual_tagint_1d;
typedef tdual_tagint_1d::t_host t_tagint_1d;
typedef tdual_tagint_1d::t_host_const t_tagint_1d_const;
typedef tdual_tagint_1d::t_host_um t_tagint_1d_um;
typedef tdual_tagint_1d::t_host_const_um t_tagint_1d_const_um;
typedef tdual_tagint_1d::t_host_const_randomread t_tagint_1d_randomread;

typedef Kokkos::
  TransformDualView<LAMMPS_NS::tagint**, device_array_layout, LMPDeviceType, void, LAMMPS_NS::tagint**, Kokkos::LayoutRight>
  tdual_tagint_2d;
typedef tdual_tagint_2d::t_host t_tagint_2d;
typedef tdual_tagint_2d::t_host_const t_tagint_2d_const;
typedef tdual_tagint_2d::t_host_um t_tagint_2d_um;
typedef tdual_tagint_2d::t_host_const_um t_tagint_2d_const_um;
typedef tdual_tagint_2d::t_host_const_randomread t_tagint_2d_randomread;

typedef Kokkos::
  DualView<LAMMPS_NS::imageint*, LMPDeviceType::array_layout, LMPDeviceType>
  tdual_imageint_1d;
typedef tdual_imageint_1d::t_host t_imageint_1d;
typedef tdual_imageint_1d::t_host_const t_imageint_1d_const;
typedef tdual_imageint_1d::t_host_um t_imageint_1d_um;
typedef tdual_imageint_1d::t_host_const_um t_imageint_1d_const_um;
typedef tdual_imageint_1d::t_host_const_randomread t_imageint_1d_randomread;

typedef Kokkos::
  TransformDualView<KK_FLOAT, LMPDeviceType::array_layout, LMPDeviceType, void, double, LMPDeviceType::array_layout>
  tdual_float_scalar;
typedef tdual_float_scalar::t_host t_float_scalar;
typedef tdual_float_scalar::t_host_const t_float_scalar_const;
typedef tdual_float_scalar::t_host_um t_float_scalar_um;
typedef tdual_float_scalar::t_host_const_um t_float_scalar_const_um;

typedef Kokkos::TransformDualView<KK_FLOAT*, LMPDeviceType::array_layout, LMPDeviceType, void, double*, LMPDeviceType::array_layout> tdual_float_1d;
typedef tdual_float_1d::t_host t_float_1d;
typedef tdual_float_1d::t_host_const t_float_1d_const;
typedef tdual_float_1d::t_host_um t_float_1d_um;
typedef tdual_float_1d::t_host_const_um t_float_1d_const_um;
typedef tdual_float_1d::t_host_const_randomread t_float_1d_randomread;

typedef Kokkos::TransformDualView<KK_FLOAT**, device_array_layout, LMPDeviceType, void, double**, Kokkos::LayoutRight> tdual_float_2d;
typedef tdual_float_2d::t_host t_float_2d;
typedef tdual_float_2d::t_host_const t_float_2d_const;
typedef tdual_float_2d::t_host_um t_float_2d_um;
typedef tdual_float_2d::t_host_const_um t_float_2d_const_um;
typedef tdual_float_2d::t_host_const_randomread t_float_2d_randomread;

typedef Kokkos::TransformDualView<KK_FLOAT**, LMPDeviceType::array_layout, LMPDeviceType, void, double**, LMPDeviceType::array_layout> tdual_float_2d_dl;
typedef tdual_float_2d_dl::t_host t_float_2d_dl;
typedef tdual_float_2d_dl::t_host_const t_float_2d_dl_const;
typedef tdual_float_2d_dl::t_host_um t_float_2d_dl_um;
typedef tdual_float_2d_dl::t_host_const_um t_float_2d_dl_const_um;
typedef tdual_float_2d_dl::t_host_const_randomread t_float_2d_dl_randomread;

typedef Kokkos::TransformDualView<KK_FLOAT**, Kokkos::LayoutRight, LMPDeviceType, void, double**, Kokkos::LayoutRight> tdual_float_2d_lr;
typedef tdual_float_2d_lr::t_host t_float_2d_lr;
typedef tdual_float_2d_lr::t_host_const t_float_2d_lr_const;
typedef tdual_float_2d_lr::t_host_um t_float_2d_lr_um;
typedef tdual_float_2d_lr::t_host_const_um t_float_2d_lr_const_um;
typedef tdual_float_2d_lr::t_host_const_randomread t_float_2d_lr_randomread;

typedef Kokkos::DualView<double**, Kokkos::LayoutRight, LMPDeviceType> tdual_double_2d_lr;
typedef tdual_double_2d_lr::t_host t_double_2d_lr;
typedef tdual_double_2d_lr::t_host_const t_double_2d_lr_const;
typedef tdual_double_2d_lr::t_host_um t_double_2d_lr_um;
typedef tdual_double_2d_lr::t_host_const_um t_double_2d_lr_const_um;
typedef tdual_double_2d_lr::t_host_const_randomread t_double_2d_lr_randomread;

typedef Kokkos::TransformDualView<KK_FLOAT*[3], device_array_layout, LMPDeviceType, void, double*[3], Kokkos::LayoutRight> tdual_float_1d_3;
typedef tdual_float_1d_3::t_host t_float_1d_3;
typedef tdual_float_1d_3::t_host_const t_float_1d_3_const;
typedef tdual_float_1d_3::t_host_um t_float_1d_3_um;
typedef tdual_float_1d_3::t_host_const_um t_float_1d_3_const_um;
typedef tdual_float_1d_3::t_host_const_randomread t_float_1d_3_randomread;

typedef Kokkos::TransformDualView<KK_FLOAT*[6], device_array_layout, LMPDeviceType, void, double*[6], Kokkos::LayoutRight> tdual_float_1d_6;
typedef tdual_float_1d_6::t_host t_float_1d_6;
typedef tdual_float_1d_6::t_host_const t_float_1d_6_const;
typedef tdual_float_1d_6::t_host_um t_float_1d_6_um;
typedef tdual_float_1d_6::t_host_const_um t_float_1d_6_const_um;
typedef tdual_float_1d_6::t_host_const_randomread t_float_1d_6_randomread;

typedef Kokkos::DualView<KK_FLOAT**[7], Kokkos::LayoutRight, LMPDeviceType> tdual_float_2d_7;
typedef tdual_float_2d_7::t_host_const_randomread t_float_2d_7_randomread;
typedef tdual_float_2d_7::t_host t_float_2d_7;

//Neighbor Types

typedef Kokkos::DualView<int**, LMPDeviceType::array_layout, LMPDeviceType> tdual_neighbors_2d;
typedef tdual_neighbors_2d::t_host t_neighbors_2d;
typedef tdual_neighbors_2d::t_host_const t_neighbors_2d_const;
typedef tdual_neighbors_2d::t_host_um t_neighbors_2d_um;
typedef tdual_neighbors_2d::t_host_const_um t_neighbors_2d_const_um;
typedef tdual_neighbors_2d::t_host_const_randomread t_neighbors_2d_randomread;

};

//default LAMMPS Types
typedef struct ArrayTypes<LAMMPS_NS::Device> DAT;
typedef struct ArrayTypes<LAMMPS_NS::Host> HAT;

template<LAMMPS_NS::ExecutionSpace Space, class BufferView, class DualView>
void buffer_view(BufferView &buf, DualView &view,
                 const size_t n0,
                 const size_t n1) {
  typedef typename GetDeviceType<Space>::value DeviceType;
  buf = BufferView(DualViewHelper<Space>::view(view).data(),n0,n1);

}

template<class DeviceType>
struct MemsetZeroFunctor {
  typedef DeviceType execution_space;
  void* ptr;
  KOKKOS_INLINE_FUNCTION void operator()(const int i) const {
    ((int*)ptr)[i] = 0;
  }
};

template<class ViewType>
void memset_kokkos (ViewType &view) {
  static MemsetZeroFunctor<typename ViewType::execution_space> f;
  f.ptr = view.data();
  #ifndef KOKKOS_USING_DEPRECATED_VIEW
  Kokkos::parallel_for(view.span()*sizeof(typename ViewType::value_type)/4, f);
  #else
  Kokkos::parallel_for(view.span()*sizeof(typename ViewType::value_type)/4, f);
  #endif
  ViewType::execution_space().fence();
}

struct params_lj_coul {
  KOKKOS_INLINE_FUNCTION
  params_lj_coul(){cut_ljsq=0;cut_coulsq=0;lj1=0;lj2=0;lj3=0;lj4=0;offset=0;};
  KOKKOS_INLINE_FUNCTION
  params_lj_coul(int i){cut_ljsq=0;cut_coulsq=0;lj1=0;lj2=0;lj3=0;lj4=0;offset=0;};
  KK_FLOAT cut_ljsq,cut_coulsq,lj1,lj2,lj3,lj4,offset;
};

// Pair SNAP

typedef KK_FLOAT SNAreal;

//typedef struct { SNAreal re, im; } SNAcomplex;
template <typename real>
struct alignas(2*sizeof(real)) SNAComplex
{
  real re,im;

  KOKKOS_FORCEINLINE_FUNCTION SNAComplex() = default;

  KOKKOS_FORCEINLINE_FUNCTION SNAComplex(real re)
   : re(re), im(static_cast<real>(0.)) { ; }

  KOKKOS_FORCEINLINE_FUNCTION SNAComplex(real re, real im)
   : re(re), im(im) { ; }

  KOKKOS_FORCEINLINE_FUNCTION SNAComplex(const SNAComplex& other)
   : re(other.re), im(other.im) { ; }

  KOKKOS_FORCEINLINE_FUNCTION SNAComplex& operator=(const SNAComplex& other) {
    re = other.re; im = other.im;
    return *this;
  }

  KOKKOS_FORCEINLINE_FUNCTION SNAComplex(SNAComplex&& other)
   : re(other.re), im(other.im) { ; }

  KOKKOS_FORCEINLINE_FUNCTION SNAComplex& operator=(SNAComplex&& other) {
    re = other.re; im = other.im;
    return *this;
  }

  KOKKOS_FORCEINLINE_FUNCTION SNAComplex operator+(SNAComplex const& other) {
    return SNAComplex(re + other.re, im + other.im);
  }

  KOKKOS_FORCEINLINE_FUNCTION SNAComplex& operator+=(SNAComplex const& other) {
    re += other.re; im += other.im;
    return *this;
  }

};

template <typename real>
KOKKOS_FORCEINLINE_FUNCTION SNAComplex<real> operator*(const real& r, const SNAComplex<real>& self) {
  return SNAComplex<real>(r*self.re, r*self.im);
}

#if defined(LMP_KOKKOS_SINGLE)
template <typename real>
KOKKOS_FORCEINLINE_FUNCTION SNAComplex<real> operator*(const double& r, const SNAComplex<real>& self) {
  return SNAComplex<real>(r*self.re, r*self.im);
}
#endif

typedef SNAComplex<SNAreal> SNAcomplex;


#if defined(KOKKOS_ENABLE_CXX11)
#undef ISFINITE
#define ISFINITE(x) std::isfinite(x)
#endif

#ifdef KOKKOS_ENABLE_CUDA
#define LAMMPS_LAMBDA [=] __host__ __device__
#else
#define LAMMPS_LAMBDA [=]
#endif

#endif
