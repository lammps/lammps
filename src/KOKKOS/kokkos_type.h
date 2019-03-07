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

#include "lmptype.h"

#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>
#include <impl/Kokkos_Timer.hpp>
#include <Kokkos_Vectorization.hpp>
#include <Kokkos_ScatterView.hpp>

enum{FULL=1u,HALFTHREAD=2u,HALF=4u,N2=8u};

#if defined(KOKKOS_ENABLE_CXX11)
#undef ISFINITE
#define ISFINITE(x) std::isfinite(x)
#endif

// User-settable FFT precision

// FFT_PRECISION = 1 is single-precision complex (4-byte real, 4-byte imag)
// FFT_PRECISION = 2 is double-precision complex (8-byte real, 8-byte imag)

#ifdef FFT_SINGLE
#define FFT_PRECISION 1
#define MPI_FFT_SCALAR MPI_FLOAT
typedef float FFT_SCALAR;
#else
#define FFT_PRECISION 2
#define MPI_FFT_SCALAR MPI_DOUBLE
typedef double FFT_SCALAR;
#endif

#define MAX_TYPES_STACKPARAMS 12
#define NeighClusterSize 8

  struct lmp_float3 {
    float x,y,z;
    KOKKOS_INLINE_FUNCTION
    lmp_float3():x(0.0f),y(0.0f),z(0.0f) {}

    KOKKOS_INLINE_FUNCTION
    void operator += (const lmp_float3& tmp) {
      x+=tmp.x;
      y+=tmp.y;
      z+=tmp.z;
    }
    KOKKOS_INLINE_FUNCTION
    void operator += (const lmp_float3& tmp) volatile {
      x+=tmp.x;
      y+=tmp.y;
      z+=tmp.z;
    }
    KOKKOS_INLINE_FUNCTION
    void operator = (const lmp_float3& tmp) {
      x=tmp.x;
      y=tmp.y;
      z=tmp.z;
    }
    KOKKOS_INLINE_FUNCTION
    void operator = (const lmp_float3& tmp) volatile {
      x=tmp.x;
      y=tmp.y;
      z=tmp.z;
    }
  };

  struct lmp_double3 {
    double x,y,z;
    KOKKOS_INLINE_FUNCTION
    lmp_double3():x(0.0),y(0.0),z(0.0) {}

    KOKKOS_INLINE_FUNCTION
    void operator += (const lmp_double3& tmp) {
      x+=tmp.x;
      y+=tmp.y;
      z+=tmp.z;
    }
    KOKKOS_INLINE_FUNCTION
    void operator += (const lmp_double3& tmp) volatile {
      x+=tmp.x;
      y+=tmp.y;
      z+=tmp.z;
    }
    KOKKOS_INLINE_FUNCTION
    void operator = (const lmp_double3& tmp) {
      x=tmp.x;
      y=tmp.y;
      z=tmp.z;
    }
    KOKKOS_INLINE_FUNCTION
    void operator = (const lmp_double3& tmp) volatile {
      x=tmp.x;
      y=tmp.y;
      z=tmp.z;
    }
  };

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

#if !defined(__CUDACC__) && !defined(__VECTOR_TYPES_H__)
  struct double2 {
    double x, y;
  };
  struct float2 {
    float x, y;
  };
  struct float4 {
    float x, y, z, w;
  };
  struct double4 {
    double x, y, z, w;
  };
#endif
// set LMPHostype and LMPDeviceType from Kokkos Default Types
typedef Kokkos::DefaultExecutionSpace LMPDeviceType;
typedef Kokkos::HostSpace::execution_space LMPHostType;

// set ExecutionSpace stuct with variable "space"

template<class Device>
struct ExecutionSpaceFromDevice;

template<>
struct ExecutionSpaceFromDevice<LMPHostType> {
  static const LAMMPS_NS::ExecutionSpace space = LAMMPS_NS::Host;
};
#ifdef KOKKOS_ENABLE_CUDA
template<>
struct ExecutionSpaceFromDevice<Kokkos::Cuda> {
  static const LAMMPS_NS::ExecutionSpace space = LAMMPS_NS::Device;
};
#endif


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


// define precision
// handle global precision, force, energy, positions, kspace separately

#ifndef PRECISION
#define PRECISION 2
#endif
#if PRECISION==1
typedef float LMP_FLOAT;
typedef float2 LMP_FLOAT2;
typedef lmp_float3 LMP_FLOAT3;
typedef float4 LMP_FLOAT4;
#else
typedef double LMP_FLOAT;
typedef double2 LMP_FLOAT2;
typedef lmp_double3 LMP_FLOAT3;
typedef double4 LMP_FLOAT4;
#endif

#ifndef PREC_FORCE
#define PREC_FORCE PRECISION
#endif

#if PREC_FORCE==1
typedef float F_FLOAT;
typedef float2 F_FLOAT2;
typedef lmp_float3 F_FLOAT3;
typedef float4 F_FLOAT4;
#else
typedef double F_FLOAT;
typedef double2 F_FLOAT2;
typedef lmp_double3 F_FLOAT3;
typedef double4 F_FLOAT4;
#endif

#ifndef PREC_ENERGY
#define PREC_ENERGY PRECISION
#endif

#if PREC_ENERGY==1
typedef float E_FLOAT;
typedef float2 E_FLOAT2;
typedef float4 E_FLOAT4;
#else
typedef double E_FLOAT;
typedef double2 E_FLOAT2;
typedef double4 E_FLOAT4;
#endif

struct s_EV_FLOAT {
  E_FLOAT evdwl;
  E_FLOAT ecoul;
  E_FLOAT v[6];
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
  E_FLOAT evdwl;
  E_FLOAT ecoul;
  E_FLOAT v[6];
  E_FLOAT ereax[10];
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

#ifndef PREC_POS
#define PREC_POS PRECISION
#endif

#if PREC_POS==1
typedef float X_FLOAT;
typedef float2 X_FLOAT2;
typedef float4 X_FLOAT4;
#else
typedef double X_FLOAT;
typedef double2 X_FLOAT2;
typedef double4 X_FLOAT4;
#endif

#ifndef PREC_VELOCITIES
#define PREC_VELOCITIES PRECISION
#endif

#if PREC_VELOCITIES==1
typedef float V_FLOAT;
typedef float2 V_FLOAT2;
typedef float4 V_FLOAT4;
#else
typedef double V_FLOAT;
typedef double2 V_FLOAT2;
typedef double4 V_FLOAT4;
#endif

#if PREC_KSPACE==1
typedef float K_FLOAT;
typedef float2 K_FLOAT2;
typedef float4 K_FLOAT4;
#else
typedef double K_FLOAT;
typedef double2 K_FLOAT2;
typedef double4 K_FLOAT4;
#endif

typedef int T_INT;

// ------------------------------------------------------------------------

// LAMMPS types

template <class DeviceType>
struct ArrayTypes;

template <>
struct ArrayTypes<LMPDeviceType> {

// scalar types

typedef Kokkos::
  DualView<int, LMPDeviceType::array_layout, LMPDeviceType> tdual_int_scalar;
typedef tdual_int_scalar::t_dev t_int_scalar;
typedef tdual_int_scalar::t_dev_const t_int_scalar_const;
typedef tdual_int_scalar::t_dev_um t_int_scalar_um;
typedef tdual_int_scalar::t_dev_const_um t_int_scalar_const_um;

typedef Kokkos::
  DualView<LMP_FLOAT, LMPDeviceType::array_layout, LMPDeviceType>
  tdual_float_scalar;
typedef tdual_float_scalar::t_dev t_float_scalar;
typedef tdual_float_scalar::t_dev_const t_float_scalar_const;
typedef tdual_float_scalar::t_dev_um t_float_scalar_um;
typedef tdual_float_scalar::t_dev_const_um t_float_scalar_const_um;

// generic array types

typedef Kokkos::
  DualView<int*, LMPDeviceType::array_layout, LMPDeviceType> tdual_int_1d;
typedef tdual_int_1d::t_dev t_int_1d;
typedef tdual_int_1d::t_dev_const t_int_1d_const;
typedef tdual_int_1d::t_dev_um t_int_1d_um;
typedef tdual_int_1d::t_dev_const_um t_int_1d_const_um;
typedef tdual_int_1d::t_dev_const_randomread t_int_1d_randomread;

typedef Kokkos::
  DualView<int*[3], Kokkos::LayoutRight, LMPDeviceType> tdual_int_1d_3;
typedef tdual_int_1d_3::t_dev t_int_1d_3;
typedef tdual_int_1d_3::t_dev_const t_int_1d_3_const;
typedef tdual_int_1d_3::t_dev_um t_int_1d_3_um;
typedef tdual_int_1d_3::t_dev_const_um t_int_1d_3_const_um;
typedef tdual_int_1d_3::t_dev_const_randomread t_int_1d_3_randomread;

typedef Kokkos::
  DualView<int**, Kokkos::LayoutRight, LMPDeviceType> tdual_int_2d;
typedef tdual_int_2d::t_dev t_int_2d;
typedef tdual_int_2d::t_dev_const t_int_2d_const;
typedef tdual_int_2d::t_dev_um t_int_2d_um;
typedef tdual_int_2d::t_dev_const_um t_int_2d_const_um;
typedef tdual_int_2d::t_dev_const_randomread t_int_2d_randomread;

typedef Kokkos::
  DualView<int**, LMPDeviceType::array_layout, LMPDeviceType> tdual_int_2d_dl;
typedef tdual_int_2d_dl::t_dev t_int_2d_dl;
typedef tdual_int_2d_dl::t_dev_const t_int_2d_const_dl;
typedef tdual_int_2d_dl::t_dev_um t_int_2d_um_dl;
typedef tdual_int_2d_dl::t_dev_const_um t_int_2d_const_um_dl;
typedef tdual_int_2d_dl::t_dev_const_randomread t_int_2d_randomread_dl;

typedef Kokkos::
  DualView<LAMMPS_NS::tagint*, LMPDeviceType::array_layout, LMPDeviceType>
  tdual_tagint_1d;
typedef tdual_tagint_1d::t_dev t_tagint_1d;
typedef tdual_tagint_1d::t_dev_const t_tagint_1d_const;
typedef tdual_tagint_1d::t_dev_um t_tagint_1d_um;
typedef tdual_tagint_1d::t_dev_const_um t_tagint_1d_const_um;
typedef tdual_tagint_1d::t_dev_const_randomread t_tagint_1d_randomread;

typedef Kokkos::
  DualView<LAMMPS_NS::tagint**, Kokkos::LayoutRight, LMPDeviceType>
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
  DualView<double*, Kokkos::LayoutRight, LMPDeviceType> tdual_double_1d;
typedef tdual_double_1d::t_dev t_double_1d;
typedef tdual_double_1d::t_dev_const t_double_1d_const;
typedef tdual_double_1d::t_dev_um t_double_1d_um;
typedef tdual_double_1d::t_dev_const_um t_double_1d_const_um;
typedef tdual_double_1d::t_dev_const_randomread t_double_1d_randomread;

typedef Kokkos::
  DualView<double**, Kokkos::LayoutRight, LMPDeviceType> tdual_double_2d;
typedef tdual_double_2d::t_dev t_double_2d;
typedef tdual_double_2d::t_dev_const t_double_2d_const;
typedef tdual_double_2d::t_dev_um t_double_2d_um;
typedef tdual_double_2d::t_dev_const_um t_double_2d_const_um;
typedef tdual_double_2d::t_dev_const_randomread t_double_2d_randomread;

// 1d float array n

typedef Kokkos::DualView<LMP_FLOAT*, LMPDeviceType::array_layout, LMPDeviceType> tdual_float_1d;
typedef tdual_float_1d::t_dev t_float_1d;
typedef tdual_float_1d::t_dev_const t_float_1d_const;
typedef tdual_float_1d::t_dev_um t_float_1d_um;
typedef tdual_float_1d::t_dev_const_um t_float_1d_const_um;
typedef tdual_float_1d::t_dev_const_randomread t_float_1d_randomread;

//2d float array n
typedef Kokkos::DualView<LMP_FLOAT**, Kokkos::LayoutRight, LMPDeviceType> tdual_float_2d;
typedef tdual_float_2d::t_dev t_float_2d;
typedef tdual_float_2d::t_dev_const t_float_2d_const;
typedef tdual_float_2d::t_dev_um t_float_2d_um;
typedef tdual_float_2d::t_dev_const_um t_float_2d_const_um;
typedef tdual_float_2d::t_dev_const_randomread t_float_2d_randomread;

//Position Types
//1d X_FLOAT array n
typedef Kokkos::DualView<X_FLOAT*, LMPDeviceType::array_layout, LMPDeviceType> tdual_xfloat_1d;
typedef tdual_xfloat_1d::t_dev t_xfloat_1d;
typedef tdual_xfloat_1d::t_dev_const t_xfloat_1d_const;
typedef tdual_xfloat_1d::t_dev_um t_xfloat_1d_um;
typedef tdual_xfloat_1d::t_dev_const_um t_xfloat_1d_const_um;
typedef tdual_xfloat_1d::t_dev_const_randomread t_xfloat_1d_randomread;

//2d X_FLOAT array n*m
typedef Kokkos::DualView<X_FLOAT**, Kokkos::LayoutRight, LMPDeviceType> tdual_xfloat_2d;
typedef tdual_xfloat_2d::t_dev t_xfloat_2d;
typedef tdual_xfloat_2d::t_dev_const t_xfloat_2d_const;
typedef tdual_xfloat_2d::t_dev_um t_xfloat_2d_um;
typedef tdual_xfloat_2d::t_dev_const_um t_xfloat_2d_const_um;
typedef tdual_xfloat_2d::t_dev_const_randomread t_xfloat_2d_randomread;

//2d X_FLOAT array n*4
#ifdef LMP_KOKKOS_NO_LEGACY
typedef Kokkos::DualView<X_FLOAT*[3], Kokkos::LayoutLeft, LMPDeviceType> tdual_x_array;
#else
typedef Kokkos::DualView<X_FLOAT*[3], Kokkos::LayoutRight, LMPDeviceType> tdual_x_array;
#endif
typedef tdual_x_array::t_dev t_x_array;
typedef tdual_x_array::t_dev_const t_x_array_const;
typedef tdual_x_array::t_dev_um t_x_array_um;
typedef tdual_x_array::t_dev_const_um t_x_array_const_um;
typedef tdual_x_array::t_dev_const_randomread t_x_array_randomread;

//Velocity Types
//1d V_FLOAT array n
typedef Kokkos::DualView<V_FLOAT*, LMPDeviceType::array_layout, LMPDeviceType> tdual_vfloat_1d;
typedef tdual_vfloat_1d::t_dev t_vfloat_1d;
typedef tdual_vfloat_1d::t_dev_const t_vfloat_1d_const;
typedef tdual_vfloat_1d::t_dev_um t_vfloat_1d_um;
typedef tdual_vfloat_1d::t_dev_const_um t_vfloat_1d_const_um;
typedef tdual_vfloat_1d::t_dev_const_randomread t_vfloat_1d_randomread;

//2d V_FLOAT array n*m
typedef Kokkos::DualView<V_FLOAT**, Kokkos::LayoutRight, LMPDeviceType> tdual_vfloat_2d;
typedef tdual_vfloat_2d::t_dev t_vfloat_2d;
typedef tdual_vfloat_2d::t_dev_const t_vfloat_2d_const;
typedef tdual_vfloat_2d::t_dev_um t_vfloat_2d_um;
typedef tdual_vfloat_2d::t_dev_const_um t_vfloat_2d_const_um;
typedef tdual_vfloat_2d::t_dev_const_randomread t_vfloat_2d_randomread;

//2d V_FLOAT array n*3
typedef Kokkos::DualView<V_FLOAT*[3], Kokkos::LayoutRight, LMPDeviceType> tdual_v_array;
//typedef Kokkos::DualView<V_FLOAT*[3], LMPDeviceType::array_layout, LMPDeviceType> tdual_v_array;
typedef tdual_v_array::t_dev t_v_array;
typedef tdual_v_array::t_dev_const t_v_array_const;
typedef tdual_v_array::t_dev_um t_v_array_um;
typedef tdual_v_array::t_dev_const_um t_v_array_const_um;
typedef tdual_v_array::t_dev_const_randomread t_v_array_randomread;

//Force Types
//1d F_FLOAT array n

typedef Kokkos::DualView<F_FLOAT*, LMPDeviceType::array_layout, LMPDeviceType> tdual_ffloat_1d;
typedef tdual_ffloat_1d::t_dev t_ffloat_1d;
typedef tdual_ffloat_1d::t_dev_const t_ffloat_1d_const;
typedef tdual_ffloat_1d::t_dev_um t_ffloat_1d_um;
typedef tdual_ffloat_1d::t_dev_const_um t_ffloat_1d_const_um;
typedef tdual_ffloat_1d::t_dev_const_randomread t_ffloat_1d_randomread;

//2d F_FLOAT array n*m

typedef Kokkos::DualView<F_FLOAT**, Kokkos::LayoutRight, LMPDeviceType> tdual_ffloat_2d;
typedef tdual_ffloat_2d::t_dev t_ffloat_2d;
typedef tdual_ffloat_2d::t_dev_const t_ffloat_2d_const;
typedef tdual_ffloat_2d::t_dev_um t_ffloat_2d_um;
typedef tdual_ffloat_2d::t_dev_const_um t_ffloat_2d_const_um;
typedef tdual_ffloat_2d::t_dev_const_randomread t_ffloat_2d_randomread;

//2d F_FLOAT array n*m, device layout

typedef Kokkos::DualView<F_FLOAT**, LMPDeviceType::array_layout, LMPDeviceType> tdual_ffloat_2d_dl;
typedef tdual_ffloat_2d_dl::t_dev t_ffloat_2d_dl;
typedef tdual_ffloat_2d_dl::t_dev_const t_ffloat_2d_const_dl;
typedef tdual_ffloat_2d_dl::t_dev_um t_ffloat_2d_um_dl;
typedef tdual_ffloat_2d_dl::t_dev_const_um t_ffloat_2d_const_um_dl;
typedef tdual_ffloat_2d_dl::t_dev_const_randomread t_ffloat_2d_randomread_dl;

//2d F_FLOAT array n*3

typedef Kokkos::DualView<F_FLOAT*[3], Kokkos::LayoutRight, LMPDeviceType> tdual_f_array;
//typedef Kokkos::DualView<F_FLOAT*[3], LMPDeviceType::array_layout, LMPDeviceType> tdual_f_array;
typedef tdual_f_array::t_dev t_f_array;
typedef tdual_f_array::t_dev_const t_f_array_const;
typedef tdual_f_array::t_dev_um t_f_array_um;
typedef tdual_f_array::t_dev_const_um t_f_array_const_um;
typedef tdual_f_array::t_dev_const_randomread t_f_array_randomread;

//2d F_FLOAT array n*6 (for virial)

typedef Kokkos::DualView<F_FLOAT*[6], Kokkos::LayoutRight, LMPDeviceType> tdual_virial_array;
typedef tdual_virial_array::t_dev t_virial_array;
typedef tdual_virial_array::t_dev_const t_virial_array_const;
typedef tdual_virial_array::t_dev_um t_virial_array_um;
typedef tdual_virial_array::t_dev_const_um t_virial_array_const_um;
typedef tdual_virial_array::t_dev_const_randomread t_virial_array_randomread;

//Energy Types
//1d E_FLOAT array n

typedef Kokkos::DualView<E_FLOAT*, LMPDeviceType::array_layout, LMPDeviceType> tdual_efloat_1d;
typedef tdual_efloat_1d::t_dev t_efloat_1d;
typedef tdual_efloat_1d::t_dev_const t_efloat_1d_const;
typedef tdual_efloat_1d::t_dev_um t_efloat_1d_um;
typedef tdual_efloat_1d::t_dev_const_um t_efloat_1d_const_um;
typedef tdual_efloat_1d::t_dev_const_randomread t_efloat_1d_randomread;

//2d E_FLOAT array n*m

typedef Kokkos::DualView<E_FLOAT**, Kokkos::LayoutRight, LMPDeviceType> tdual_efloat_2d;
typedef tdual_efloat_2d::t_dev t_efloat_2d;
typedef tdual_efloat_2d::t_dev_const t_efloat_2d_const;
typedef tdual_efloat_2d::t_dev_um t_efloat_2d_um;
typedef tdual_efloat_2d::t_dev_const_um t_efloat_2d_const_um;
typedef tdual_efloat_2d::t_dev_const_randomread t_efloat_2d_randomread;

//2d E_FLOAT array n*3

typedef Kokkos::DualView<E_FLOAT*[3], Kokkos::LayoutRight, LMPDeviceType> tdual_e_array;
typedef tdual_e_array::t_dev t_e_array;
typedef tdual_e_array::t_dev_const t_e_array_const;
typedef tdual_e_array::t_dev_um t_e_array_um;
typedef tdual_e_array::t_dev_const_um t_e_array_const_um;
typedef tdual_e_array::t_dev_const_randomread t_e_array_randomread;

//Neighbor Types

typedef Kokkos::DualView<int**, LMPDeviceType::array_layout, LMPDeviceType> tdual_neighbors_2d;
typedef tdual_neighbors_2d::t_dev t_neighbors_2d;
typedef tdual_neighbors_2d::t_dev_const t_neighbors_2d_const;
typedef tdual_neighbors_2d::t_dev_um t_neighbors_2d_um;
typedef tdual_neighbors_2d::t_dev_const_um t_neighbors_2d_const_um;
typedef tdual_neighbors_2d::t_dev_const_randomread t_neighbors_2d_randomread;

//Kspace

typedef Kokkos::
  DualView<FFT_SCALAR*, Kokkos::LayoutRight, LMPDeviceType> tdual_FFT_SCALAR_1d;
typedef tdual_FFT_SCALAR_1d::t_dev t_FFT_SCALAR_1d;
typedef tdual_FFT_SCALAR_1d::t_dev_um t_FFT_SCALAR_1d_um;

typedef Kokkos::DualView<FFT_SCALAR**,Kokkos::LayoutRight,LMPDeviceType> tdual_FFT_SCALAR_2d;
typedef tdual_FFT_SCALAR_2d::t_dev t_FFT_SCALAR_2d;

typedef Kokkos::DualView<FFT_SCALAR**[3],Kokkos::LayoutRight,LMPDeviceType> tdual_FFT_SCALAR_2d_3;
typedef tdual_FFT_SCALAR_2d_3::t_dev t_FFT_SCALAR_2d_3;

typedef Kokkos::DualView<FFT_SCALAR***,Kokkos::LayoutRight,LMPDeviceType> tdual_FFT_SCALAR_3d;
typedef tdual_FFT_SCALAR_3d::t_dev t_FFT_SCALAR_3d;

typedef Kokkos::
  DualView<FFT_SCALAR*[2], Kokkos::LayoutRight, LMPDeviceType> tdual_FFT_DATA_1d;
typedef tdual_FFT_DATA_1d::t_dev t_FFT_DATA_1d;
typedef tdual_FFT_DATA_1d::t_dev_um t_FFT_DATA_1d_um;

typedef Kokkos::
  DualView<int*, LMPDeviceType::array_layout, LMPDeviceType> tdual_int_64;
typedef tdual_int_64::t_dev t_int_64;
typedef tdual_int_64::t_dev_um t_int_64_um;

};

#ifdef KOKKOS_ENABLE_CUDA
template <>
struct ArrayTypes<LMPHostType> {

//Scalar Types

typedef Kokkos::DualView<int, LMPDeviceType::array_layout, LMPDeviceType> tdual_int_scalar;
typedef tdual_int_scalar::t_host t_int_scalar;
typedef tdual_int_scalar::t_host_const t_int_scalar_const;
typedef tdual_int_scalar::t_host_um t_int_scalar_um;
typedef tdual_int_scalar::t_host_const_um t_int_scalar_const_um;

typedef Kokkos::DualView<LMP_FLOAT, LMPDeviceType::array_layout, LMPDeviceType> tdual_float_scalar;
typedef tdual_float_scalar::t_host t_float_scalar;
typedef tdual_float_scalar::t_host_const t_float_scalar_const;
typedef tdual_float_scalar::t_host_um t_float_scalar_um;
typedef tdual_float_scalar::t_host_const_um t_float_scalar_const_um;

//Generic ArrayTypes
typedef Kokkos::DualView<int*, LMPDeviceType::array_layout, LMPDeviceType> tdual_int_1d;
typedef tdual_int_1d::t_host t_int_1d;
typedef tdual_int_1d::t_host_const t_int_1d_const;
typedef tdual_int_1d::t_host_um t_int_1d_um;
typedef tdual_int_1d::t_host_const_um t_int_1d_const_um;
typedef tdual_int_1d::t_host_const_randomread t_int_1d_randomread;

typedef Kokkos::DualView<int*[3], Kokkos::LayoutRight, LMPDeviceType> tdual_int_1d_3;
typedef tdual_int_1d_3::t_host t_int_1d_3;
typedef tdual_int_1d_3::t_host_const t_int_1d_3_const;
typedef tdual_int_1d_3::t_host_um t_int_1d_3_um;
typedef tdual_int_1d_3::t_host_const_um t_int_1d_3_const_um;
typedef tdual_int_1d_3::t_host_const_randomread t_int_1d_3_randomread;

typedef Kokkos::DualView<int**, Kokkos::LayoutRight, LMPDeviceType> tdual_int_2d;
typedef tdual_int_2d::t_host t_int_2d;
typedef tdual_int_2d::t_host_const t_int_2d_const;
typedef tdual_int_2d::t_host_um t_int_2d_um;
typedef tdual_int_2d::t_host_const_um t_int_2d_const_um;
typedef tdual_int_2d::t_host_const_randomread t_int_2d_randomread;

typedef Kokkos::DualView<int**, LMPDeviceType::array_layout, LMPDeviceType> tdual_int_2d_dl;
typedef tdual_int_2d_dl::t_host t_int_2d_dl;
typedef tdual_int_2d_dl::t_host_const t_int_2d_const_dl;
typedef tdual_int_2d_dl::t_host_um t_int_2d_um_dl;
typedef tdual_int_2d_dl::t_host_const_um t_int_2d_const_um_dl;
typedef tdual_int_2d_dl::t_host_const_randomread t_int_2d_randomread_dl;

typedef Kokkos::DualView<LAMMPS_NS::tagint*, LMPDeviceType::array_layout, LMPDeviceType> tdual_tagint_1d;
typedef tdual_tagint_1d::t_host t_tagint_1d;
typedef tdual_tagint_1d::t_host_const t_tagint_1d_const;
typedef tdual_tagint_1d::t_host_um t_tagint_1d_um;
typedef tdual_tagint_1d::t_host_const_um t_tagint_1d_const_um;
typedef tdual_tagint_1d::t_host_const_randomread t_tagint_1d_randomread;

typedef Kokkos::
  DualView<LAMMPS_NS::tagint**, Kokkos::LayoutRight, LMPDeviceType>
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
  DualView<double*, Kokkos::LayoutRight, LMPDeviceType> tdual_double_1d;
typedef tdual_double_1d::t_host t_double_1d;
typedef tdual_double_1d::t_host_const t_double_1d_const;
typedef tdual_double_1d::t_host_um t_double_1d_um;
typedef tdual_double_1d::t_host_const_um t_double_1d_const_um;
typedef tdual_double_1d::t_host_const_randomread t_double_1d_randomread;

typedef Kokkos::
  DualView<double**, Kokkos::LayoutRight, LMPDeviceType> tdual_double_2d;
typedef tdual_double_2d::t_host t_double_2d;
typedef tdual_double_2d::t_host_const t_double_2d_const;
typedef tdual_double_2d::t_host_um t_double_2d_um;
typedef tdual_double_2d::t_host_const_um t_double_2d_const_um;
typedef tdual_double_2d::t_host_const_randomread t_double_2d_randomread;

//1d float array n
typedef Kokkos::DualView<LMP_FLOAT*, LMPDeviceType::array_layout, LMPDeviceType> tdual_float_1d;
typedef tdual_float_1d::t_host t_float_1d;
typedef tdual_float_1d::t_host_const t_float_1d_const;
typedef tdual_float_1d::t_host_um t_float_1d_um;
typedef tdual_float_1d::t_host_const_um t_float_1d_const_um;
typedef tdual_float_1d::t_host_const_randomread t_float_1d_randomread;

//2d float array n
typedef Kokkos::DualView<LMP_FLOAT**, Kokkos::LayoutRight, LMPDeviceType> tdual_float_2d;
typedef tdual_float_2d::t_host t_float_2d;
typedef tdual_float_2d::t_host_const t_float_2d_const;
typedef tdual_float_2d::t_host_um t_float_2d_um;
typedef tdual_float_2d::t_host_const_um t_float_2d_const_um;
typedef tdual_float_2d::t_host_const_randomread t_float_2d_randomread;

//Position Types
//1d X_FLOAT array n
typedef Kokkos::DualView<X_FLOAT*, LMPDeviceType::array_layout, LMPDeviceType> tdual_xfloat_1d;
typedef tdual_xfloat_1d::t_host t_xfloat_1d;
typedef tdual_xfloat_1d::t_host_const t_xfloat_1d_const;
typedef tdual_xfloat_1d::t_host_um t_xfloat_1d_um;
typedef tdual_xfloat_1d::t_host_const_um t_xfloat_1d_const_um;
typedef tdual_xfloat_1d::t_host_const_randomread t_xfloat_1d_randomread;

//2d X_FLOAT array n*m
typedef Kokkos::DualView<X_FLOAT**, Kokkos::LayoutRight, LMPDeviceType> tdual_xfloat_2d;
typedef tdual_xfloat_2d::t_host t_xfloat_2d;
typedef tdual_xfloat_2d::t_host_const t_xfloat_2d_const;
typedef tdual_xfloat_2d::t_host_um t_xfloat_2d_um;
typedef tdual_xfloat_2d::t_host_const_um t_xfloat_2d_const_um;
typedef tdual_xfloat_2d::t_host_const_randomread t_xfloat_2d_randomread;

//2d X_FLOAT array n*3
typedef Kokkos::DualView<X_FLOAT*[3], Kokkos::LayoutRight, LMPDeviceType> tdual_x_array;
typedef tdual_x_array::t_host t_x_array;
typedef tdual_x_array::t_host_const t_x_array_const;
typedef tdual_x_array::t_host_um t_x_array_um;
typedef tdual_x_array::t_host_const_um t_x_array_const_um;
typedef tdual_x_array::t_host_const_randomread t_x_array_randomread;

//Velocity Types
//1d V_FLOAT array n
typedef Kokkos::DualView<V_FLOAT*, LMPDeviceType::array_layout, LMPDeviceType> tdual_vfloat_1d;
typedef tdual_vfloat_1d::t_host t_vfloat_1d;
typedef tdual_vfloat_1d::t_host_const t_vfloat_1d_const;
typedef tdual_vfloat_1d::t_host_um t_vfloat_1d_um;
typedef tdual_vfloat_1d::t_host_const_um t_vfloat_1d_const_um;
typedef tdual_vfloat_1d::t_host_const_randomread t_vfloat_1d_randomread;

//2d V_FLOAT array n*m
typedef Kokkos::DualView<V_FLOAT**, Kokkos::LayoutRight, LMPDeviceType> tdual_vfloat_2d;
typedef tdual_vfloat_2d::t_host t_vfloat_2d;
typedef tdual_vfloat_2d::t_host_const t_vfloat_2d_const;
typedef tdual_vfloat_2d::t_host_um t_vfloat_2d_um;
typedef tdual_vfloat_2d::t_host_const_um t_vfloat_2d_const_um;
typedef tdual_vfloat_2d::t_host_const_randomread t_vfloat_2d_randomread;

//2d V_FLOAT array n*3
typedef Kokkos::DualView<V_FLOAT*[3], Kokkos::LayoutRight, LMPDeviceType> tdual_v_array;
//typedef Kokkos::DualView<V_FLOAT*[3], LMPDeviceType::array_layout, LMPDeviceType> tdual_v_array;
typedef tdual_v_array::t_host t_v_array;
typedef tdual_v_array::t_host_const t_v_array_const;
typedef tdual_v_array::t_host_um t_v_array_um;
typedef tdual_v_array::t_host_const_um t_v_array_const_um;
typedef tdual_v_array::t_host_const_randomread t_v_array_randomread;

//Force Types
//1d F_FLOAT array n
typedef Kokkos::DualView<F_FLOAT*, LMPDeviceType::array_layout, LMPDeviceType> tdual_ffloat_1d;
typedef tdual_ffloat_1d::t_host t_ffloat_1d;
typedef tdual_ffloat_1d::t_host_const t_ffloat_1d_const;
typedef tdual_ffloat_1d::t_host_um t_ffloat_1d_um;
typedef tdual_ffloat_1d::t_host_const_um t_ffloat_1d_const_um;
typedef tdual_ffloat_1d::t_host_const_randomread t_ffloat_1d_randomread;

//2d F_FLOAT array n*m
typedef Kokkos::DualView<F_FLOAT**, Kokkos::LayoutRight, LMPDeviceType> tdual_ffloat_2d;
typedef tdual_ffloat_2d::t_host t_ffloat_2d;
typedef tdual_ffloat_2d::t_host_const t_ffloat_2d_const;
typedef tdual_ffloat_2d::t_host_um t_ffloat_2d_um;
typedef tdual_ffloat_2d::t_host_const_um t_ffloat_2d_const_um;
typedef tdual_ffloat_2d::t_host_const_randomread t_ffloat_2d_randomread;

//2d F_FLOAT array n*m, device layout
typedef Kokkos::DualView<F_FLOAT**, LMPDeviceType::array_layout, LMPDeviceType> tdual_ffloat_2d_dl;
typedef tdual_ffloat_2d_dl::t_host t_ffloat_2d_dl;
typedef tdual_ffloat_2d_dl::t_host_const t_ffloat_2d_const_dl;
typedef tdual_ffloat_2d_dl::t_host_um t_ffloat_2d_um_dl;
typedef tdual_ffloat_2d_dl::t_host_const_um t_ffloat_2d_const_um_dl;
typedef tdual_ffloat_2d_dl::t_host_const_randomread t_ffloat_2d_randomread_dl;

//2d F_FLOAT array n*3
typedef Kokkos::DualView<F_FLOAT*[3], Kokkos::LayoutRight, LMPDeviceType> tdual_f_array;
//typedef Kokkos::DualView<F_FLOAT*[3], LMPDeviceType::array_layout, LMPDeviceType> tdual_f_array;
typedef tdual_f_array::t_host t_f_array;
typedef tdual_f_array::t_host_const t_f_array_const;
typedef tdual_f_array::t_host_um t_f_array_um;
typedef tdual_f_array::t_host_const_um t_f_array_const_um;
typedef tdual_f_array::t_host_const_randomread t_f_array_randomread;

//2d F_FLOAT array n*6 (for virial)
typedef Kokkos::DualView<F_FLOAT*[6], Kokkos::LayoutRight, LMPDeviceType> tdual_virial_array;
typedef tdual_virial_array::t_host t_virial_array;
typedef tdual_virial_array::t_host_const t_virial_array_const;
typedef tdual_virial_array::t_host_um t_virial_array_um;
typedef tdual_virial_array::t_host_const_um t_virial_array_const_um;
typedef tdual_virial_array::t_host_const_randomread t_virial_array_randomread;



//Energy Types
//1d E_FLOAT array n
typedef Kokkos::DualView<E_FLOAT*, LMPDeviceType::array_layout, LMPDeviceType> tdual_efloat_1d;
typedef tdual_efloat_1d::t_host t_efloat_1d;
typedef tdual_efloat_1d::t_host_const t_efloat_1d_const;
typedef tdual_efloat_1d::t_host_um t_efloat_1d_um;
typedef tdual_efloat_1d::t_host_const_um t_efloat_1d_const_um;
typedef tdual_efloat_1d::t_host_const_randomread t_efloat_1d_randomread;

//2d E_FLOAT array n*m
typedef Kokkos::DualView<E_FLOAT**, Kokkos::LayoutRight, LMPDeviceType> tdual_efloat_2d;
typedef tdual_efloat_2d::t_host t_efloat_2d;
typedef tdual_efloat_2d::t_host_const t_efloat_2d_const;
typedef tdual_efloat_2d::t_host_um t_efloat_2d_um;
typedef tdual_efloat_2d::t_host_const_um t_efloat_2d_const_um;
typedef tdual_efloat_2d::t_host_const_randomread t_efloat_2d_randomread;

//2d E_FLOAT array n*3
typedef Kokkos::DualView<E_FLOAT*[3], Kokkos::LayoutRight, LMPDeviceType> tdual_e_array;
typedef tdual_e_array::t_host t_e_array;
typedef tdual_e_array::t_host_const t_e_array_const;
typedef tdual_e_array::t_host_um t_e_array_um;
typedef tdual_e_array::t_host_const_um t_e_array_const_um;
typedef tdual_e_array::t_host_const_randomread t_e_array_randomread;

//Neighbor Types
typedef Kokkos::DualView<int**, LMPDeviceType::array_layout, LMPDeviceType> tdual_neighbors_2d;
typedef tdual_neighbors_2d::t_host t_neighbors_2d;
typedef tdual_neighbors_2d::t_host_const t_neighbors_2d_const;
typedef tdual_neighbors_2d::t_host_um t_neighbors_2d_um;
typedef tdual_neighbors_2d::t_host_const_um t_neighbors_2d_const_um;
typedef tdual_neighbors_2d::t_host_const_randomread t_neighbors_2d_randomread;


//Kspace

typedef Kokkos::
  DualView<FFT_SCALAR*, Kokkos::LayoutRight, LMPDeviceType> tdual_FFT_SCALAR_1d;
typedef tdual_FFT_SCALAR_1d::t_host t_FFT_SCALAR_1d;
typedef tdual_FFT_SCALAR_1d::t_host_um t_FFT_SCALAR_1d_um;

typedef Kokkos::DualView<FFT_SCALAR**,Kokkos::LayoutRight,LMPDeviceType> tdual_FFT_SCALAR_2d;
typedef tdual_FFT_SCALAR_2d::t_host t_FFT_SCALAR_2d;

typedef Kokkos::DualView<FFT_SCALAR**[3],Kokkos::LayoutRight,LMPDeviceType> tdual_FFT_SCALAR_2d_3;
typedef tdual_FFT_SCALAR_2d_3::t_host t_FFT_SCALAR_2d_3;

typedef Kokkos::DualView<FFT_SCALAR***,Kokkos::LayoutRight,LMPDeviceType> tdual_FFT_SCALAR_3d;
typedef tdual_FFT_SCALAR_3d::t_host t_FFT_SCALAR_3d;

typedef Kokkos::
  DualView<FFT_SCALAR*[2], Kokkos::LayoutRight, LMPDeviceType> tdual_FFT_DATA_1d;
typedef tdual_FFT_DATA_1d::t_host t_FFT_DATA_1d;
typedef tdual_FFT_DATA_1d::t_host_um t_FFT_DATA_1d_um;

typedef Kokkos::
  DualView<int*, LMPDeviceType::array_layout, LMPDeviceType> tdual_int_64;
typedef tdual_int_64::t_host t_int_64;
typedef tdual_int_64::t_host_um t_int_64_um;

};
#endif
//default LAMMPS Types
typedef struct ArrayTypes<LMPDeviceType> DAT;
typedef struct ArrayTypes<LMPHostType> HAT;

template<class DeviceType, class BufferView, class DualView>
void buffer_view(BufferView &buf, DualView &view,
                 const size_t n0,
                 const size_t n1 = 0,
                 const size_t n2 = 0,
                 const size_t n3 = 0,
                 const size_t n4 = 0,
                 const size_t n5 = 0,
                 const size_t n6 = 0,
                 const size_t n7 = 0) {

  buf = BufferView(
          view.template view<DeviceType>().data(),
          n0,n1,n2,n3,n4,n5,n6,n7);

}

template<class DeviceType>
struct MemsetZeroFunctor {
  typedef DeviceType  execution_space ;
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
  Kokkos::parallel_for(view.capacity()*sizeof(typename ViewType::value_type)/4, f);
  #endif
  ViewType::execution_space::fence();
}

struct params_lj_coul {
  KOKKOS_INLINE_FUNCTION
  params_lj_coul(){cut_ljsq=0;cut_coulsq=0;lj1=0;lj2=0;lj3=0;lj4=0;offset=0;};
  KOKKOS_INLINE_FUNCTION
  params_lj_coul(int i){cut_ljsq=0;cut_coulsq=0;lj1=0;lj2=0;lj3=0;lj4=0;offset=0;};
  F_FLOAT cut_ljsq,cut_coulsq,lj1,lj2,lj3,lj4,offset;
};

#if defined(KOKKOS_ENABLE_CXX11)
#undef ISFINITE
#define ISFINITE(x) std::isfinite(x)
#endif

#ifdef KOKKOS_ENABLE_CUDA
#define LAMMPS_LAMBDA [=] __device__
#else
#define LAMMPS_LAMBDA [=]
#endif

#endif
