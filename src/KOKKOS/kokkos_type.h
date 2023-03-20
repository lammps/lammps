// clang-format off
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

#ifndef LMP_LMPTYPE_KOKKOS_H
#define LMP_LMPTYPE_KOKKOS_H

#include "pointers.h"
#include "lmptype.h"

#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>
#include <Kokkos_Timer.hpp>
#include <Kokkos_ScatterView.hpp>
#include <Kokkos_UnorderedMap.hpp>

constexpr int FULL = 1;
constexpr int HALFTHREAD = 2;
constexpr int HALF = 4;

#if defined(KOKKOS_ENABLE_CXX11)
#undef ISFINITE
#define ISFINITE(x) std::isfinite(x)
#endif

#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP) || defined(KOKKOS_ENABLE_SYCL) || defined(KOKKOS_ENABLE_OPENMPTARGET)
#define LMP_KOKKOS_GPU
#endif

#if defined(LMP_KOKKOS_GPU)
#define KOKKOS_GPU_ARG(x) x
#else
#define KOKKOS_GPU_ARG(x)
#endif

#define MAX_TYPES_STACKPARAMS 12
static constexpr LAMMPS_NS::bigint LMP_KOKKOS_AV_DELTA = 10;

namespace Kokkos {
  static auto NoInit = [](std::string const& label) {
    return Kokkos::view_alloc(Kokkos::WithoutInitializing, label);
  };
}

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

// set LMPHostype and LMPDeviceType from Kokkos Default Types
typedef Kokkos::DefaultExecutionSpace LMPDeviceType;
typedef Kokkos::HostSpace::execution_space LMPHostType;


// Need to use Cuda UVM memory space for Host execution space

template<class DeviceType>
class KKDevice {
public:
#if defined(KOKKOS_ENABLE_CUDA) && defined(KOKKOS_ENABLE_CUDA_UVM)
  typedef Kokkos::Device<DeviceType,LMPDeviceType::memory_space> value;
#else
  typedef Kokkos::Device<DeviceType,typename DeviceType::memory_space> value;
#endif
};

// Helpers for readability

using KKScatterSum = Kokkos::Experimental::ScatterSum;
using KKScatterDuplicated = Kokkos::Experimental::ScatterDuplicated;
using KKScatterNonDuplicated = Kokkos::Experimental::ScatterNonDuplicated;

template<typename DataType, typename Layout, typename Device, typename... Args>
using KKScatterView = Kokkos::Experimental::ScatterView<DataType, Layout, Device, Args...>;


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
#elif defined(KOKKOS_ENABLE_HIP)
template<>
struct ExecutionSpaceFromDevice<Kokkos::Experimental::HIP> {
  static const LAMMPS_NS::ExecutionSpace space = LAMMPS_NS::Device;
};
#elif defined(KOKKOS_ENABLE_SYCL)
template<>
struct ExecutionSpaceFromDevice<Kokkos::Experimental::SYCL> {
  static const LAMMPS_NS::ExecutionSpace space = LAMMPS_NS::Device;
};
#elif defined(KOKKOS_ENABLE_OPENMPTARGET)
template<>
struct ExecutionSpaceFromDevice<Kokkos::Experimental::OpenMPTarget> {
  static const LAMMPS_NS::ExecutionSpace space = LAMMPS_NS::Device;
};
#endif

// set host pinned space
#if defined(KOKKOS_ENABLE_CUDA)
typedef Kokkos::CudaHostPinnedSpace LMPPinnedHostType;
#elif defined(KOKKOS_ENABLE_HIP)
typedef Kokkos::Experimental::HIPHostPinnedSpace LMPPinnedHostType;
#elif defined(KOKKOS_ENABLE_SYCL)
typedef Kokkos::Experimental::SYCLHostUSMSpace LMPPinnedHostType;
#elif defined(KOKKOS_ENABLE_OPENMPTARGET)
typedef Kokkos::Serial LMPPinnedHostType;
#else
typedef LMPHostType LMPPinnedHostType;
#endif

// create simple LMPDeviceSpace typedef for non CUDA-, HIP-, or SYCL-specific
// behaviour
#if defined(KOKKOS_ENABLE_CUDA)
typedef Kokkos::Cuda LMPDeviceSpace;
#elif defined(KOKKOS_ENABLE_HIP)
typedef Kokkos::Experimental::HIP LMPDeviceSpace;
#elif defined(KOKKOS_ENABLE_SYCL)
typedef Kokkos::Experimental::SYCL LMPDeviceSpace;
#elif defined(KOKKOS_ENABLE_OPENMPTARGET)
typedef Kokkos::Experimental::OpenMPTarget LMPDeviceSpace;
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
  using value = Kokkos::Experimental::ScatterNonAtomic;
};

template<int NEIGHFLAG, class DeviceType>
using AtomicDup_v = typename AtomicDup<NEIGHFLAG, DeviceType>::value;

#ifdef KOKKOS_ENABLE_CUDA
template<>
struct AtomicDup<HALFTHREAD,Kokkos::Cuda> {
  using value = Kokkos::Experimental::ScatterAtomic;
};
#elif defined(KOKKOS_ENABLE_HIP)
template<>
struct AtomicDup<HALFTHREAD,Kokkos::Experimental::HIP> {
  using value = Kokkos::Experimental::ScatterAtomic;
};
#elif defined(KOKKOS_ENABLE_SYCL)
template<>
struct AtomicDup<HALFTHREAD,Kokkos::Experimental::SYCL> {
  using value = Kokkos::Experimental::ScatterAtomic;
};
#elif defined(KOKKOS_ENABLE_OPENMPTARGET)
template<>
struct AtomicDup<HALFTHREAD,Kokkos::Experimental::OpenMPTarget> {
  using value = Kokkos::Experimental::ScatterAtomic;
};
#endif

#ifdef LMP_KOKKOS_USE_ATOMICS

#ifdef KOKKOS_ENABLE_OPENMP
template<>
struct AtomicDup<HALFTHREAD,Kokkos::OpenMP> {
  using value = Kokkos::Experimental::ScatterAtomic;
};
#endif

#ifdef KOKKOS_ENABLE_THREADS
template<>
struct AtomicDup<HALFTHREAD,Kokkos::Threads> {
  using value = Kokkos::Experimental::ScatterAtomic;
};
#endif

#endif


// Determine duplication traits for force array
// Use duplication when running threaded and not using atomics
template<int NEIGHFLAG, class DeviceType>
struct NeedDup {
  using value = Kokkos::Experimental::ScatterNonDuplicated;
};

template<int NEIGHFLAG, class DeviceType>
using NeedDup_v = typename NeedDup<NEIGHFLAG,DeviceType>::value;

#ifndef LMP_KOKKOS_USE_ATOMICS

#ifdef KOKKOS_ENABLE_OPENMP
template<>
struct NeedDup<HALFTHREAD,Kokkos::OpenMP> {
  using value = Kokkos::Experimental::ScatterDuplicated;
};
#endif

#ifdef KOKKOS_ENABLE_THREADS
template<>
struct NeedDup<HALFTHREAD,Kokkos::Threads> {
  using value = Kokkos::Experimental::ScatterDuplicated;
};
#endif

#endif

template<typename value, typename T1, typename T2>
class ScatterViewHelper {};

template<typename T1, typename T2>
class ScatterViewHelper<Kokkos::Experimental::ScatterDuplicated,T1,T2> {
public:
  KOKKOS_INLINE_FUNCTION
  static T1 get(const T1 &dup, const T2 & /*nondup*/) {
    return dup;
  }
};

template<typename T1, typename T2>
class ScatterViewHelper<Kokkos::Experimental::ScatterNonDuplicated,T1,T2> {
public:
  KOKKOS_INLINE_FUNCTION
  static T2 get(const T1 & /*dup*/, const T2 &nondup) {
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
#else
typedef double LMP_FLOAT;
#endif

#ifndef PREC_FORCE
#define PREC_FORCE PRECISION
#endif

#if PREC_FORCE==1
typedef float F_FLOAT;
#else
typedef double F_FLOAT;
#endif

#ifndef PREC_ENERGY
#define PREC_ENERGY PRECISION
#endif

#if PREC_ENERGY==1
typedef float E_FLOAT;
#else
typedef double E_FLOAT;
#endif

struct s_EV_FLOAT {
  E_FLOAT evdwl;
  E_FLOAT ecoul;
  E_FLOAT v[6];
  KOKKOS_INLINE_FUNCTION
  s_EV_FLOAT() {
    evdwl = 0;
    ecoul = 0;
    for (int i = 0; i < 6; ++i)
      v[i] = 0;
  }

  KOKKOS_INLINE_FUNCTION
  void operator+=(const s_EV_FLOAT &rhs) {
    evdwl += rhs.evdwl;
    ecoul += rhs.ecoul;
    for (int i = 0; i < 6; ++i)
      v[i] += rhs.v[i];
  }

  KOKKOS_INLINE_FUNCTION
  void operator+=(const volatile s_EV_FLOAT &rhs) volatile {
    evdwl += rhs.evdwl;
    ecoul += rhs.ecoul;
    for (int i = 0; i < 6; ++i)
      v[i] += rhs.v[i];
  }
};
typedef struct s_EV_FLOAT EV_FLOAT;

struct s_EV_FLOAT_REAX {
  E_FLOAT evdwl;
  E_FLOAT ecoul;
  E_FLOAT v[6];
  E_FLOAT ereax[9];
  KOKKOS_INLINE_FUNCTION
  s_EV_FLOAT_REAX() {
    evdwl = 0;
    ecoul = 0;
    for (int i = 0; i < 6; ++i)
      v[i] = 0;
    for (int i = 0; i < 9; ++i)
      ereax[i] = 0;
  }

  KOKKOS_INLINE_FUNCTION
  void operator+=(const s_EV_FLOAT_REAX &rhs) {
    evdwl += rhs.evdwl;
    ecoul += rhs.ecoul;
    for (int i = 0; i < 6; ++i)
      v[i] += rhs.v[i];
    for (int i = 0; i < 9; ++i)
      ereax[i] += rhs.ereax[i];
  }

  KOKKOS_INLINE_FUNCTION
  void operator+=(const volatile s_EV_FLOAT_REAX &rhs) volatile {
    evdwl += rhs.evdwl;
    ecoul += rhs.ecoul;
    for (int i = 0; i < 6; ++i)
      v[i] += rhs.v[i];
    for (int i = 0; i < 9; ++i)
      ereax[i] += rhs.ereax[i];
  }
};
typedef struct s_EV_FLOAT_REAX EV_FLOAT_REAX;

struct s_FEV_FLOAT {
  F_FLOAT f[3];
  E_FLOAT evdwl;
  E_FLOAT ecoul;
  E_FLOAT v[6];
  KOKKOS_INLINE_FUNCTION
  s_FEV_FLOAT() {
    evdwl = 0;
    ecoul = 0;
    for (int i = 0; i < 6; ++i)
      v[i] = 0;
    for (int i = 0; i < 3; ++i)
      f[i] = 0;
  }

  KOKKOS_INLINE_FUNCTION
  void operator+=(const s_FEV_FLOAT &rhs) {
    evdwl += rhs.evdwl;
    ecoul += rhs.ecoul;
    for (int i = 0; i < 6; ++i)
      v[i] += rhs.v[i];
    for (int i = 0; i < 3; ++i)
      f[i] += rhs.f[i];
  }

  KOKKOS_INLINE_FUNCTION
  void operator+=(const volatile s_FEV_FLOAT &rhs) volatile {
    evdwl += rhs.evdwl;
    ecoul += rhs.ecoul;
    for (int i = 0; i < 6; ++i)
      v[i] += rhs.v[i];
    for (int i = 0; i < 3; ++i)
      f[i] += rhs.f[i];
  }
};
typedef struct s_FEV_FLOAT FEV_FLOAT;

struct alignas(2*sizeof(F_FLOAT)) s_FLOAT2 {
  F_FLOAT v[2];

  KOKKOS_INLINE_FUNCTION
  s_FLOAT2() {
    v[0] = v[1] = 0.0;
  }

  KOKKOS_INLINE_FUNCTION
  s_FLOAT2(const s_FLOAT2 & rhs) {
    for (int i = 0; i < 2; i++){
      v[i] = rhs.v[i];
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator+=(const s_FLOAT2 &rhs) {
    v[0] += rhs.v[0];
    v[1] += rhs.v[1];
  }

  KOKKOS_INLINE_FUNCTION
  void operator+=(const volatile s_FLOAT2 &rhs) volatile {
    v[0] += rhs.v[0];
    v[1] += rhs.v[1];
  }
};
typedef struct s_FLOAT2 F_FLOAT2;

#ifndef PREC_POS
#define PREC_POS PRECISION
#endif

#if PREC_POS==1
typedef float X_FLOAT;
#else
typedef double X_FLOAT;
#endif

#ifndef PREC_VELOCITIES
#define PREC_VELOCITIES PRECISION
#endif

#if PREC_VELOCITIES==1
typedef float V_FLOAT;
#else
typedef double V_FLOAT;
#endif

#if PREC_KSPACE==1
typedef float K_FLOAT;
#else
typedef double K_FLOAT;
#endif

typedef int T_INT;

// ------------------------------------------------------------------------

// LAMMPS types

typedef Kokkos::UnorderedMap<LAMMPS_NS::tagint,int,LMPDeviceType> hash_type;
typedef hash_type::HostMirror host_hash_type;

struct dual_hash_type {
  hash_type d_view;
  host_hash_type h_view;

  template<class DeviceType>
  KOKKOS_INLINE_FUNCTION
  std::enable_if_t<(std::is_same<DeviceType,LMPDeviceType>::value || Kokkos::SpaceAccessibility<LMPDeviceType::memory_space,LMPHostType::memory_space>::accessible),hash_type&> view() {return d_view;}

  template<class DeviceType>
  KOKKOS_INLINE_FUNCTION
  std::enable_if_t<!(std::is_same<DeviceType,LMPDeviceType>::value || Kokkos::SpaceAccessibility<LMPDeviceType::memory_space,LMPHostType::memory_space>::accessible),host_hash_type&> view() {return h_view;}

};

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

//3d float array n
typedef Kokkos::DualView<LMP_FLOAT***, Kokkos::LayoutRight, LMPDeviceType> tdual_float_3d;
typedef tdual_float_3d::t_dev t_float_3d;
typedef tdual_float_3d::t_dev_const t_float_3d_const;
typedef tdual_float_3d::t_dev_um t_float_3d_um;
typedef tdual_float_3d::t_dev_const_um t_float_3d_const_um;
typedef tdual_float_3d::t_dev_const_randomread t_float_3d_randomread;

#ifdef LMP_KOKKOS_NO_LEGACY
typedef Kokkos::DualView<X_FLOAT*[4], Kokkos::LayoutLeft, LMPDeviceType> tdual_float_1d_4;
#else
typedef Kokkos::DualView<X_FLOAT*[4], Kokkos::LayoutRight, LMPDeviceType> tdual_float_1d_4;
#endif

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

// 1d F_FLOAT2 array n
typedef Kokkos::DualView<F_FLOAT*[2], Kokkos::LayoutRight, LMPDeviceType> tdual_ffloat2_1d;
typedef tdual_ffloat2_1d::t_dev t_ffloat2_1d;
typedef tdual_ffloat2_1d::t_dev_const t_ffloat2_1d_const;
typedef tdual_ffloat2_1d::t_dev_um t_ffloat2_1d_um;
typedef tdual_ffloat2_1d::t_dev_const_um t_ffloat2_1d_const_um;
typedef tdual_ffloat2_1d::t_dev_const_randomread t_ffloat2_1d_randomread;

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

//2d F_FLOAT array n*4 (for dipoles and quaterions)

typedef tdual_float_1d_4::t_dev t_mu_array;
typedef tdual_float_1d_4::t_dev_const t_mu_array_const;
typedef tdual_float_1d_4::t_dev_um t_mu_array_um;
typedef tdual_float_1d_4::t_dev_const_um t_mu_array_const_um;
typedef tdual_float_1d_4::t_dev_const_randomread t_mu_array_randomread;

//2d F_FLOAT array n*6 (for virial)

typedef Kokkos::DualView<F_FLOAT*[6], Kokkos::LayoutRight, LMPDeviceType> tdual_virial_array;
typedef tdual_virial_array::t_dev t_virial_array;
typedef tdual_virial_array::t_dev_const t_virial_array_const;
typedef tdual_virial_array::t_dev_um t_virial_array_um;
typedef tdual_virial_array::t_dev_const_um t_virial_array_const_um;
typedef tdual_virial_array::t_dev_const_randomread t_virial_array_randomread;

// Spin Types

//3d SP_FLOAT array n*4

typedef tdual_float_1d_4::t_dev t_sp_array;
typedef tdual_float_1d_4::t_dev_const t_sp_array_const;
typedef tdual_float_1d_4::t_dev_um t_sp_array_um;
typedef tdual_float_1d_4::t_dev_const_um t_sp_array_const_um;
typedef tdual_float_1d_4::t_dev_const_randomread t_sp_array_randomread;

//3d FM_FLOAT array n*3

typedef tdual_f_array::t_dev t_fm_array;
typedef tdual_f_array::t_dev_const t_fm_array_const;
typedef tdual_f_array::t_dev_um t_fm_array_um;
typedef tdual_f_array::t_dev_const_um t_fm_array_const_um;
typedef tdual_f_array::t_dev_const_randomread t_fm_array_randomread;

//3d FML_FLOAT array n*3

typedef tdual_f_array::t_dev t_fm_long_array;
typedef tdual_f_array::t_dev_const t_fm_long_array_const;
typedef tdual_f_array::t_dev_um t_fm_long_array_um;
typedef tdual_f_array::t_dev_const_um t_fm_long_array_const_um;
typedef tdual_f_array::t_dev_const_randomread t_fm_long_array_randomread;

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

typedef Kokkos::DualView<int**, Kokkos::LayoutRight, LMPDeviceType> tdual_neighbors_2d_lr;
typedef tdual_neighbors_2d_lr::t_dev t_neighbors_2d_lr;
typedef tdual_neighbors_2d_lr::t_dev_const t_neighbors_2d_const_lr;
typedef tdual_neighbors_2d_lr::t_dev_um t_neighbors_2d_um_lr;
typedef tdual_neighbors_2d_lr::t_dev_const_um t_neighbors_2d_const_um_lr;
typedef tdual_neighbors_2d_lr::t_dev_const_randomread t_neighbors_2d_randomread_lr;

};

#ifdef LMP_KOKKOS_GPU
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

#ifdef LMP_KOKKOS_NO_LEGACY
typedef Kokkos::DualView<X_FLOAT*[4], Kokkos::LayoutLeft, LMPDeviceType> tdual_float_1d_4;
#else
typedef Kokkos::DualView<X_FLOAT*[4], Kokkos::LayoutRight, LMPDeviceType> tdual_float_1d_4;
#endif

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

// 1d F_FLOAT2 array n
typedef Kokkos::DualView<F_FLOAT*[2], Kokkos::LayoutRight, LMPDeviceType> tdual_ffloat2_1d;
typedef tdual_ffloat2_1d::t_host t_ffloat2_1d;
typedef tdual_ffloat2_1d::t_host_const t_ffloat2_1d_const;
typedef tdual_ffloat2_1d::t_host_um t_ffloat2_1d_um;
typedef tdual_ffloat2_1d::t_host_const_um t_ffloat2_1d_const_um;
typedef tdual_ffloat2_1d::t_host_const_randomread t_ffloat2_1d_randomread;

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

//2d F_FLOAT array n*4 (for dipoles and quaterions)

typedef tdual_float_1d_4::t_host t_mu_array;
typedef tdual_float_1d_4::t_host_const t_mu_array_const;
typedef tdual_float_1d_4::t_host_um t_mu_array_um;
typedef tdual_float_1d_4::t_host_const_um t_mu_array_const_um;
typedef tdual_float_1d_4::t_host_const_randomread t_mu_array_randomread;

//2d F_FLOAT array n*6 (for virial)
typedef Kokkos::DualView<F_FLOAT*[6], Kokkos::LayoutRight, LMPDeviceType> tdual_virial_array;
typedef tdual_virial_array::t_host t_virial_array;
typedef tdual_virial_array::t_host_const t_virial_array_const;
typedef tdual_virial_array::t_host_um t_virial_array_um;
typedef tdual_virial_array::t_host_const_um t_virial_array_const_um;
typedef tdual_virial_array::t_host_const_randomread t_virial_array_randomread;

// Spin types

//2d X_FLOAT array n*4
typedef tdual_float_1d_4::t_host t_sp_array;
typedef tdual_float_1d_4::t_host_const t_sp_array_const;
typedef tdual_float_1d_4::t_host_um t_sp_array_um;
typedef tdual_float_1d_4::t_host_const_um t_sp_array_const_um;
typedef tdual_float_1d_4::t_host_const_randomread t_sp_array_randomread;

//2d F_FLOAT array n*3
typedef tdual_f_array::t_host t_fm_array;
typedef tdual_f_array::t_host_const t_fm_array_const;
typedef tdual_f_array::t_host_um t_fm_array_um;
typedef tdual_f_array::t_host_const_um t_fm_array_const_um;
typedef tdual_f_array::t_host_const_randomread t_fm_array_randomread;

//2d F_FLOAT array n*3
typedef tdual_f_array::t_host t_fm_long_array;
typedef tdual_f_array::t_host_const t_fm_long_array_const;
typedef tdual_f_array::t_host_um t_fm_long_array_um;
typedef tdual_f_array::t_host_const_um t_fm_long_array_const_um;
typedef tdual_f_array::t_host_const_randomread t_fm_long_array_randomread;


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

typedef Kokkos::DualView<int**, Kokkos::LayoutRight, LMPDeviceType> tdual_neighbors_2d_lr;
typedef tdual_neighbors_2d_lr::t_host t_neighbors_2d_lr;
typedef tdual_neighbors_2d_lr::t_host_const t_neighbors_2d_const_lr;
typedef tdual_neighbors_2d_lr::t_host_um t_neighbors_2d_um_lr;
typedef tdual_neighbors_2d_lr::t_host_const_um t_neighbors_2d_const_um_lr;
typedef tdual_neighbors_2d_lr::t_host_const_randomread t_neighbors_2d_randomread_lr;

};
#endif
//default LAMMPS Types
typedef struct ArrayTypes<LMPDeviceType> DAT;
typedef struct ArrayTypes<LMPHostType> HAT;

template<class DeviceType, class BufferView, class DualView>
void buffer_view(BufferView &buf, DualView &view,
                 const size_t n0,
                 const size_t n1) {

  buf = BufferView(view.template view<DeviceType>().data(),n0,n1);

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
  Kokkos::parallel_for(view.span()*sizeof(typename ViewType::value_type)/4, f);
  #endif
  ViewType::execution_space().fence();
}

struct params_lj_coul {
  KOKKOS_INLINE_FUNCTION
  params_lj_coul() {cut_ljsq=0;cut_coulsq=0;lj1=0;lj2=0;lj3=0;lj4=0;offset=0;};
  KOKKOS_INLINE_FUNCTION
  params_lj_coul(int /*i*/) {cut_ljsq=0;cut_coulsq=0;lj1=0;lj2=0;lj3=0;lj4=0;offset=0;};
  F_FLOAT cut_ljsq,cut_coulsq,lj1,lj2,lj3,lj4,offset;
};

// ReaxFF

struct alignas(4 * sizeof(int)) reax_int4 {
  int i0, i1, i2, i3;
};

// Pair SNAP

#define SNAP_KOKKOS_REAL double
#define SNAP_KOKKOS_HOST_VECLEN 1

#ifdef LMP_KOKKOS_GPU
#define SNAP_KOKKOS_DEVICE_VECLEN 32
#else
#define SNAP_KOKKOS_DEVICE_VECLEN 1
#endif


// intentional: SNAreal/complex gets reused beyond SNAP
typedef double SNAreal;

//typedef struct { SNAreal re, im; } SNAcomplex;
template <typename real_type_>
struct alignas(2*sizeof(real_type_)) SNAComplex
{
  using real_type = real_type_;
  using complex = SNAComplex<real_type>;
  real_type re,im;

  KOKKOS_FORCEINLINE_FUNCTION SNAComplex()
   : re(static_cast<real_type>(0.)), im(static_cast<real_type>(0.)) { ; }

  KOKKOS_FORCEINLINE_FUNCTION SNAComplex(real_type re)
   : re(re), im(static_cast<real_type>(0.)) { ; }

  KOKKOS_FORCEINLINE_FUNCTION SNAComplex(real_type re, real_type im)
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

  KOKKOS_INLINE_FUNCTION
  static constexpr complex zero() { return complex(static_cast<real_type>(0.), static_cast<real_type>(0.)); }

  KOKKOS_INLINE_FUNCTION
  static constexpr complex one() { return complex(static_cast<real_type>(1.), static_cast<real_type>(0.)); }

  KOKKOS_INLINE_FUNCTION
  const complex conj() const { return complex(re, -im); }

  KOKKOS_INLINE_FUNCTION
  const real_type real_part_product(const complex &cm2) { return re * cm2.re - im * cm2.im; }

  KOKKOS_INLINE_FUNCTION
  const real_type real_part_product(const real_type &r) const { return re * r; }
};

template <typename real_type>
KOKKOS_FORCEINLINE_FUNCTION SNAComplex<real_type> operator*(const real_type& r, const SNAComplex<real_type>& self) {
  return SNAComplex<real_type>(r*self.re, r*self.im);
}

template <typename real_type>
KOKKOS_FORCEINLINE_FUNCTION SNAComplex<real_type> operator*(const SNAComplex<real_type>& self, const real_type& r) {
  return SNAComplex<real_type>(r*self.re, r*self.im);
}

template <typename real_type>
KOKKOS_FORCEINLINE_FUNCTION SNAComplex<real_type> operator*(const SNAComplex<real_type>& self, const SNAComplex<real_type>& cm2) {
  return SNAComplex<real_type>(self.re*cm2.re - self.im*cm2.im, self.re*cm2.im + self.im*cm2.re);
}

typedef SNAComplex<SNAreal> SNAcomplex;

#if defined(KOKKOS_ENABLE_CXX11)
#undef ISFINITE
#define ISFINITE(x) std::isfinite(x)
#endif

#define LAMMPS_LAMBDA KOKKOS_LAMBDA

#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
#define LAMMPS_DEVICE_FUNCTION __device__
#else
#define LAMMPS_DEVICE_FUNCTION
#endif

#ifdef LMP_KOKKOS_GPU
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__) || defined(__SYCL_DEVICE_ONLY__)
#define LMP_KK_DEVICE_COMPILE
#endif
#endif

#endif
