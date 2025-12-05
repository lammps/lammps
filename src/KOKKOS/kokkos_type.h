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

namespace LAMMPS_NS {
  union d_ubuf {
    double d;
    int64_t i;
    KOKKOS_INLINE_FUNCTION d_ubuf(double arg) : d(arg) {}
    KOKKOS_INLINE_FUNCTION d_ubuf(int64_t arg) : i(arg) {}
    KOKKOS_INLINE_FUNCTION d_ubuf(int arg) : i(arg) {}
  };
}

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
    void operator = (const lmp_float3& tmp) {
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
    void operator = (const lmp_double3& tmp) {
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
  t_scalar3 operator+= (const t_scalar3& rhs) {
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

// set default device layout

#if !defined (LMP_KOKKOS_LAYOUT_LEGACY) && !defined (LMP_KOKKOS_LAYOUT_DEFAULT)
#define LMP_KOKKOS_LAYOUT_LEGACY
#endif

#if defined(LMP_KOKKOS_LAYOUT_LEGACY)
typedef Kokkos::LayoutRight LMPDeviceLayout;
#else
typedef LMPDeviceType::array_layout LMPDeviceLayout;
//typedef Kokkos::LayoutLeft LMPDeviceLayout;
#endif

// If unified memory, need to use device memory space for host execution space

template<class DeviceType>
class KKDevice {
 public:
#if ((defined(KOKKOS_ENABLE_CUDA) && defined(KOKKOS_ENABLE_CUDA_UVM)) || \
     (defined(KOKKOS_ENABLE_HIP) && defined(KOKKOS_ARCH_AMD_GFX942_APU)))
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


// set ExecutionSpace struct with variable "space"

// Warning: the templated version here returns the Kokkos host
//  "HostKK" space. If the legacy Host space is required instead, do not use
//  these templates (use "Host" space directly)

template<class Device>
struct ExecutionSpaceFromDevice;

template<>
struct ExecutionSpaceFromDevice<LMPHostType> {
  static const LAMMPS_NS::ExecutionSpace space = LAMMPS_NS::HostKK;
};

#ifdef KOKKOS_ENABLE_CUDA
template<>
struct ExecutionSpaceFromDevice<Kokkos::Cuda> {
  static const LAMMPS_NS::ExecutionSpace space = LAMMPS_NS::Device;
};
#elif defined(KOKKOS_ENABLE_HIP)
template<>
struct ExecutionSpaceFromDevice<Kokkos::HIP> {
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
typedef Kokkos::HIPHostPinnedSpace LMPPinnedHostType;
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
typedef Kokkos::HIP LMPDeviceSpace;
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
struct AtomicDup<HALFTHREAD,Kokkos::HIP> {
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

// LMP_KOKKOS_SINGLE_SINGLE: Single precision for all calculations
// LMP_KOKKOS_DOUBLE_DOUBLE: Double precision for all calculations
// LMP_KOKKOS_SINGLE_DOUBLE: Accumulation of forces, etc. in double

#if !(defined LMP_KOKKOS_SINGLE_SINGLE) && !(defined LMP_KOKKOS_DOUBLE_DOUBLE) \
  && !(defined LMP_KOKKOS_SINGLE_DOUBLE)
#define LMP_KOKKOS_DOUBLE_DOUBLE
#endif

#if defined (LMP_KOKKOS_SINGLE_SINGLE) // single
typedef float KK_FLOAT;
typedef float KK_ACC_FLOAT;
#elif defined (LMP_KOKKOS_DOUBLE_DOUBLE)  // double
typedef double KK_FLOAT;
typedef double KK_ACC_FLOAT;
#elif defined (LMP_KOKKOS_SINGLE_DOUBLE) // mixed
typedef float KK_FLOAT;
typedef double KK_ACC_FLOAT;
#endif

struct s_EV_FLOAT {
  KK_ACC_FLOAT evdwl;
  KK_ACC_FLOAT ecoul;
  KK_ACC_FLOAT v[6];
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
};
typedef struct s_EV_FLOAT EV_FLOAT;

struct s_EV_FLOAT_REAX {
  KK_ACC_FLOAT evdwl;
  KK_ACC_FLOAT ecoul;
  KK_ACC_FLOAT v[6];
  KK_ACC_FLOAT ereax[9];
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
};
typedef struct s_EV_FLOAT_REAX EV_FLOAT_REAX;

struct s_FEV_FLOAT {
  KK_ACC_FLOAT f[3];
  KK_ACC_FLOAT evdwl;
  KK_ACC_FLOAT ecoul;
  KK_ACC_FLOAT v[6];
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
};
typedef struct s_FEV_FLOAT FEV_FLOAT;

struct alignas(2*sizeof(double)) s_KK_double2 {
  double v[2];

  KOKKOS_INLINE_FUNCTION
  s_KK_double2() {
    v[0] = v[1] = 0.0;
  }

  KOKKOS_INLINE_FUNCTION
  void operator+=(const s_KK_double2 &rhs) {
    v[0] += rhs.v[0];
    v[1] += rhs.v[1];
  }
};
typedef struct s_KK_double2 KK_double2;

struct alignas(2*sizeof(KK_FLOAT)) s_KK_FLOAT2 {
  KK_FLOAT v[2];

  KOKKOS_INLINE_FUNCTION
  s_KK_FLOAT2() {
    v[0] = v[1] = 0;
  }

  KOKKOS_INLINE_FUNCTION
  void operator+=(const s_KK_FLOAT2 &rhs) {
    v[0] += rhs.v[0];
    v[1] += rhs.v[1];
  }
};
typedef struct s_KK_FLOAT2 KK_FLOAT2;

template <class KeyViewType>
struct BinOp3DLAMMPS {
  int max_bins_[3] = {};
  double mul_[3]   = {};
  double min_[3]   = {};

  BinOp3DLAMMPS() = default;

  BinOp3DLAMMPS(int max_bins__[], double min[],
          double max[]) {
    max_bins_[0] = max_bins__[0];
    max_bins_[1] = max_bins__[1];
    max_bins_[2] = max_bins__[2];
    mul_[0]      = static_cast<double>(max_bins__[0]) /
              (static_cast<double>(max[0]) - static_cast<double>(min[0]));
    mul_[1] = static_cast<double>(max_bins__[1]) /
              (static_cast<double>(max[1]) - static_cast<double>(min[1]));
    mul_[2] = static_cast<double>(max_bins__[2]) /
              (static_cast<double>(max[2]) - static_cast<double>(min[2]));
    min_[0] = static_cast<double>(min[0]);
    min_[1] = static_cast<double>(min[1]);
    min_[2] = static_cast<double>(min[2]);
  }

  template <class ViewType>
  KOKKOS_INLINE_FUNCTION int bin(ViewType& keys, const int& i) const {
    int ix = static_cast<int> ((keys(i, 0) - min_[0]) * mul_[0]);
    int iy = static_cast<int> ((keys(i, 1) - min_[1]) * mul_[1]);
    int iz = static_cast<int> ((keys(i, 2) - min_[2]) * mul_[2]);
    ix = MAX(ix,0);
    iy = MAX(iy,0);
    iz = MAX(iz,0);
    ix = MIN(ix,max_bins_[0]-1);
    iy = MIN(iy,max_bins_[1]-1);
    iz = MIN(iz,max_bins_[2]-1);
    const int ibin = iz*max_bins_[1]*max_bins_[0] + iy*max_bins_[0] + ix;
    return ibin;
  }

  KOKKOS_INLINE_FUNCTION
  int max_bins() const { return max_bins_[0] * max_bins_[1] * max_bins_[2]; }

  template <class ViewType, typename iType1, typename iType2>
  KOKKOS_INLINE_FUNCTION bool operator()(ViewType& keys, iType1& i1,
                                         iType2& i2) const {
    if (keys(i1, 2) > keys(i2, 2))
      return true;
    else if (keys(i1, 2) == keys(i2, 2)) {
      if (keys(i1, 1) > keys(i2, 1))
        return true;
      else if (keys(i1, 1) == keys(i2, 1)) {
        if (keys(i1, 0) > keys(i2, 0)) return true;
      }
    }
    return false;
  }
};

typedef int T_INT;

// ------------------------------------------------------------------------

// LAMMPS types

typedef Kokkos::UnorderedMap<LAMMPS_NS::tagint,int,LMPDeviceType> hash_type;
typedef hash_type::HostMirror host_hash_type;

struct dual_hash_type {

 private:
  hash_type d_view;
  host_hash_type h_view;

 public:
  bool modified_device;
  bool modified_host;

  dual_hash_type() {
    modified_device = modified_host = false;
    d_view = hash_type();
    h_view = host_hash_type();
 }

  dual_hash_type(int capacity) {
    modified_device = modified_host = false;
    d_view = hash_type(capacity);
    h_view = host_hash_type(capacity);
 }

  template<class DeviceType>
  std::enable_if_t<(std::is_same_v<DeviceType,LMPDeviceType> || Kokkos::SpaceAccessibility<LMPDeviceType::memory_space,LMPHostType::memory_space>::accessible),hash_type&> view() {return d_view;}

  template<class DeviceType>
  std::enable_if_t<!(std::is_same_v<DeviceType,LMPDeviceType> || Kokkos::SpaceAccessibility<LMPDeviceType::memory_space,LMPHostType::memory_space>::accessible),host_hash_type&> view() {return h_view;}

  template<class DeviceType>
  KOKKOS_INLINE_FUNCTION
  std::enable_if_t<(std::is_same_v<DeviceType,LMPDeviceType> || Kokkos::SpaceAccessibility<LMPDeviceType::memory_space,LMPHostType::memory_space>::accessible),const hash_type&> const_view() const {return d_view;}

  template<class DeviceType>
  KOKKOS_INLINE_FUNCTION
  std::enable_if_t<!(std::is_same_v<DeviceType,LMPDeviceType> || Kokkos::SpaceAccessibility<LMPDeviceType::memory_space,LMPHostType::memory_space>::accessible),const host_hash_type&> const_view() const {return h_view;}

  void modify_device()
  {
    modified_device = true;
    if (modified_device && modified_host)
      Kokkos::abort("Concurrent modification of host and device hashes");
  }

  void modify_host()
  {
    modified_host = true;
    if (modified_device && modified_host)
      Kokkos::abort("Concurrent modification of host and device hashes");
  }

  void sync_device()
  {
    if (modified_host) {
      Kokkos::deep_copy(d_view,h_view);
      modified_host = false;
    }
  }

  void sync_host()
  {
    if (modified_device) {
      Kokkos::deep_copy(h_view,d_view);
      modified_device = false;
    }
  }

  template<class DeviceType>
  std::enable_if_t<(std::is_same<DeviceType,LMPDeviceType>::value || Kokkos::SpaceAccessibility<LMPDeviceType::memory_space,LMPHostType::memory_space>::accessible),void> sync() {sync_device();}

  template<class DeviceType>
  std::enable_if_t<!(std::is_same<DeviceType,LMPDeviceType>::value || Kokkos::SpaceAccessibility<LMPDeviceType::memory_space,LMPHostType::memory_space>::accessible),void> sync() {sync_host();}

  KOKKOS_INLINE_FUNCTION
  host_hash_type& view_host() { return h_view; }

  KOKKOS_INLINE_FUNCTION
  hash_type& view_device() { return d_view; }

};


template<class KKType, class LegacyType, class KKLayout, class KKSpace = LMPDeviceType>
struct TransformView {

  static constexpr int NEED_TRANSFORM = !(std::is_same<KKType,LegacyType>::value && std::is_same<KKLayout,Kokkos::LayoutRight>::value);

  typedef Kokkos::DualView<KKType, KKLayout, KKSpace> kk_view;
  typedef typename Kokkos::DualView<LegacyType, Kokkos::LayoutRight, KKSpace>::t_host legacy_view;

 private:
  kk_view k_view;
  typename kk_view::t_dev d_view;
  typename kk_view::t_host h_viewkk;
  legacy_view h_view;

 public:
  // does the DualView have only one device

  static constexpr int SINGLE_DEVICE =
    std::is_same_v<typename kk_view::t_dev::device_type, typename kk_view::t_host::device_type>;

  typedef typename legacy_view::value_type value_type;
  typedef typename legacy_view::array_layout array_layout;

  typedef Kokkos::View<typename kk_view::data_type,
               typename kk_view::array_layout,
               typename std::conditional<
                 std::is_same_v<typename kk_view::execution_space,LMPDeviceType>,
                 LMPPinnedHostType,typename kk_view::memory_space>::type,
               Kokkos::MemoryTraits<Kokkos::Unmanaged> > pinned_mirror_type;

  int modified_hostkk_legacy; // Kokkos host was modified wrt legacy host
  int modified_device_legacy; // device was modified wrt legacy host
  int modified_legacy_hostkk; // legacy host was modified wrt Kokkos host
  int modified_legacy_device; // legacy host was modified wrt device

  TransformView() {
    modified_hostkk_legacy = 0;
    modified_device_legacy = 0;
    modified_legacy_hostkk = 0;
    modified_legacy_device = 0;
    k_view = {};
    d_view = {};
    h_viewkk = {};
    h_view = {};
  }

  template <typename... Indices>
  TransformView(std::string name, Indices... ns) {
    modified_hostkk_legacy = 0;
    modified_device_legacy = 0;
    modified_legacy_hostkk = 0;
    modified_legacy_device = 0;
    k_view = kk_view(name, ns...);
    d_view = k_view.view_device();
    h_viewkk = k_view.view_host();
    if constexpr (NEED_TRANSFORM)
      h_view = legacy_view(name, ns...);
    else
      h_view = h_viewkk;
  }

  template <typename... Indices>
  void resize(Indices... ns) {
    if (modified_legacy_device) {
      sync_hostkk();
      k_view.modify_host(); // force resize on host
    }

    k_view.resize(ns...);
    d_view = k_view.view_device();
    h_viewkk = k_view.view_host();
    if constexpr (NEED_TRANSFORM) {
      Kokkos::resize(h_view,ns...);
      if (k_view.need_sync_host()) {
        modify_legacy_hostkk();
      } else if (k_view.need_sync_device()) {
        modify_legacy_device();
      }
    } else
      h_view = h_viewkk;
  }

  // mark device as modified wrt legacy

  void modify_device_legacy()
  {
    if (SINGLE_DEVICE) return modify_hostkk_legacy();

    if constexpr (NEED_TRANSFORM) {
      if (!d_view.data()) return;

      modified_device_legacy = 1;

      if (modified_legacy_device) {
        std::string msg = "TransformView::modify ERROR: ";
        msg += "Concurrent modification of legacy host and device views ";
        msg += "in TransformView \"";
        msg += d_view.label();
        msg += "\"\n";
        Kokkos::abort(msg.c_str());
      }
    }
  }

  // mark device as modified wrt all

  void modify_device()
  {
    if (SINGLE_DEVICE) return modify_hostkk();

    k_view.modify_device();
    modify_device_legacy();
  }

  void modify_hostkk_legacy() {
    if constexpr (NEED_TRANSFORM) {
      if (!h_viewkk.data()) return;

      modified_hostkk_legacy = 1;

      if (modified_legacy_hostkk) {
        std::string msg = "TransformView::modify ERROR: ";
        msg += "Concurrent modification of legacy host and Kokkos host views ";
        msg += "in TransformView \"";
        msg += d_view.label();
        msg += "\"\n";
        Kokkos::abort(msg.c_str());
      }
    }
  }

  // mark Kokkos host as modified wrt all

  void modify_hostkk()
  {
    k_view.modify_host();
    modify_hostkk_legacy();
  }

  void modify_legacy_device()
  {
    if (SINGLE_DEVICE) return modify_legacy_hostkk();

    if constexpr (NEED_TRANSFORM) {
      if (!h_view.data()) return;

      modified_legacy_device = 1;

      if (modified_device_legacy) {
        std::string msg = "TransformView::modify ERROR: ";
        msg += "Concurrent modification of device and legacy host views ";
        msg += "in TransformView \"";
        msg += d_view.label();
        msg += "\"\n";
        Kokkos::abort(msg.c_str());
      }
    }
  }

  // mark legacy host modified wrt Kokkos host

  void modify_legacy_hostkk()
  {
    if constexpr (NEED_TRANSFORM) {
      if (!h_view.data()) return;

      modified_legacy_hostkk = 1;

      if (modified_hostkk_legacy) {
        std::string msg = "TransformView::modify ERROR: ";
        msg += "Concurrent modification of Kokkos host and legacy host views ";
        msg += "in TransformView \"";
        msg += d_view.label();
        msg += "\"\n";
        Kokkos::abort(msg.c_str());
      }
    }
  }

  // mark legacy host as modified wrt all

  void modify_host() {
    if constexpr (NEED_TRANSFORM) {
      if (!SINGLE_DEVICE)
        modify_legacy_device();
      modify_legacy_hostkk();
    } else {
     modify_hostkk();
    }
  }

  void sync_legacy_to_device(void* buffer = nullptr, int async_flag = 0) {

    if (SINGLE_DEVICE) return sync_legacy_to_hostkk();

    if constexpr (NEED_TRANSFORM) {
      if (!d_view.data()) return;

      if (modified_legacy_device) {
        if (buffer) {
          pinned_mirror_type tmp_view((typename kk_view::value_type*)buffer, d_view.layout());
          Kokkos::deep_copy(LMPHostType(),tmp_view,h_view);
          Kokkos::deep_copy(LMPHostType(),d_view,tmp_view);

          if (modified_legacy_hostkk || modified_hostkk_legacy)
            k_view.modify_device();
          else
            k_view.clear_sync_state();

          if (!async_flag) Kokkos::fence();
        } else {
          Kokkos::deep_copy(h_viewkk,h_view);
          k_view.clear_sync_state();
          k_view.modify_host();
          k_view.sync_device();
          modified_legacy_hostkk = 0;
          modified_hostkk_legacy = 0;
        }
        modified_legacy_device = 0;
      }
    }
  }

  // sync all to device

  void sync_device(void* buffer = nullptr, int async_flag = 0)
  {
    if (SINGLE_DEVICE) return sync_hostkk(buffer,async_flag);

    if (!d_view.data()) return;

    // prevent double copy

    if constexpr (NEED_TRANSFORM) {
      if (modified_legacy_device && k_view.need_sync_device()) {
        if (modified_hostkk_legacy) {
          modified_legacy_device = 0;
        } else if (modified_legacy_hostkk) {
          k_view.clear_sync_state();
        }
      }
    }

    if (buffer && k_view.need_sync_device()) {
      pinned_mirror_type tmp_view((typename kk_view::value_type*)buffer, d_view.layout());
      Kokkos::deep_copy(LMPHostType(),tmp_view,h_viewkk);
      Kokkos::deep_copy(LMPHostType(),d_view,tmp_view);
      k_view.clear_sync_state();
      if (!async_flag) Kokkos::fence();
    } else {
      k_view.sync_device();
    }

    sync_legacy_to_device(buffer,async_flag);
  }

  // sync legacy host to Kokkos host

  void sync_legacy_to_hostkk() {
    if constexpr (NEED_TRANSFORM) {
      if (!h_viewkk.data()) return;

      if (modified_legacy_hostkk) {
        Kokkos::deep_copy(h_viewkk,h_view);
        if (modified_legacy_device || modified_device_legacy)
          k_view.modify_host();
        else
          k_view.clear_sync_state();
        modified_legacy_hostkk = 0;
      }
    }
  }

  // sync all to Kokkos host

  void sync_hostkk(void* buffer = nullptr, int async_flag = 0)
  {
    if (!h_viewkk.data()) return;

    // prevent double copy

    if constexpr (NEED_TRANSFORM) {
      if (modified_legacy_hostkk && k_view.need_sync_host()) {
        if (modified_legacy_device) {
          k_view.clear_sync_state();
        } else if (modified_device_legacy) {
          modified_legacy_hostkk = 0;
        }
      }
    }

    if (buffer && k_view.need_sync_host()) {
      pinned_mirror_type tmp_view((typename kk_view::value_type*)buffer, d_view.layout());
      Kokkos::deep_copy(LMPHostType(),tmp_view,d_view);
      Kokkos::deep_copy(LMPHostType(),h_viewkk,tmp_view);
      k_view.clear_sync_state();
      if (!async_flag) Kokkos::fence();
    } else {
      k_view.sync_host();
    }
    sync_legacy_to_hostkk();
  }

  // sync device to legacy host

  void sync_device_to_legacy(void* buffer = nullptr, int async_flag = 0)
  {
    if (SINGLE_DEVICE) return sync_hostkk_to_legacy();

    if constexpr (NEED_TRANSFORM) {
      if (!h_view.data()) return;

      if (modified_device_legacy) {
        if (buffer) {
          pinned_mirror_type tmp_view((typename kk_view::value_type*)buffer, d_view.layout());
          Kokkos::deep_copy(LMPHostType(),tmp_view,d_view);
          Kokkos::deep_copy(LMPHostType(),h_view,tmp_view);

          if (k_view.need_sync_host() || k_view.need_sync_device())
            modify_legacy_hostkk();
          else
            modified_legacy_hostkk = 0;

          if (!async_flag) Kokkos::fence();
        } else {
          k_view.modify_device();
          k_view.sync_host();
          Kokkos::deep_copy(h_view,h_viewkk);
          modified_hostkk_legacy = 0;
          modified_legacy_hostkk = 0;
        }
        modified_device_legacy = 0;
      }
    }
  }

  // sync Kokkos host to legacy host

  void sync_hostkk_to_legacy()
  {
    if constexpr (NEED_TRANSFORM) {

      if (!h_view.data()) return;

      if (modified_hostkk_legacy) {
        Kokkos::deep_copy(h_view,h_viewkk);
        modified_hostkk_legacy = 0;
        if (k_view.need_sync_host() || k_view.need_sync_device()) {
          if (SINGLE_DEVICE)
            modify_legacy_device();
        } else
          modified_legacy_device = 0;
      }
    }
  }

  // sync all to legacy host

  void sync_host(void* buffer = nullptr, int async_flag = 0) {
    if constexpr (NEED_TRANSFORM) {
      if (!h_view.data()) return;

      // prevent double copy

      if (modified_device_legacy && modified_hostkk_legacy) {
        if (k_view.need_sync_device()) {
          modified_device_legacy = 0;
        } else if (k_view.need_sync_host()) {
          modified_hostkk_legacy = 0;
        }
      }

      if (!SINGLE_DEVICE)
        sync_device_to_legacy(buffer,async_flag);
      sync_hostkk_to_legacy();
    } else {
      sync_hostkk(buffer,async_flag);
    }
  }

  // Warning: these templated calls (e.g. t_view.view<DeviceType>()
  //  and t_view.sync<DeviceType>()) specify Kokkos host "HostKK" views
  //  (i.e. h_viewkk). Use non-templated calls (e.g. t_view.view_host() and
  //  t_view.sync_host()) to specify legacy "Host" views (h_view)
  //  instead

  template <typename iType>
  KOKKOS_INLINE_FUNCTION constexpr std::enable_if_t<std::is_integral_v<iType>,
                                                    size_t>
  extent(const iType& r) const {
    return k_view.extent(r);
  }

  template <typename iType>
  KOKKOS_INLINE_FUNCTION constexpr std::enable_if_t<std::is_integral_v<iType>,
                                                    int>
  extent_int(const iType& r) const {
    return static_cast<int>(k_view.extent(r));
  }

  template<class DeviceType>
  std::enable_if_t<(std::is_same_v<DeviceType,LMPDeviceType> || Kokkos::SpaceAccessibility<LMPDeviceType::memory_space,LMPHostType::memory_space>::accessible),typename kk_view::t_dev&> view() {return d_view;}

  template<class DeviceType>
  std::enable_if_t<!(std::is_same_v<DeviceType,LMPDeviceType> || Kokkos::SpaceAccessibility<LMPDeviceType::memory_space,LMPHostType::memory_space>::accessible),typename kk_view::t_host&> view() {return h_viewkk;}

  template<class DeviceType>
  KOKKOS_INLINE_FUNCTION std::enable_if_t<(std::is_same_v<DeviceType,LMPDeviceType> || Kokkos::SpaceAccessibility<LMPDeviceType::memory_space,LMPHostType::memory_space>::accessible),const typename kk_view::t_dev&> view() const {return d_view;}

  template<class DeviceType>
  KOKKOS_INLINE_FUNCTION std::enable_if_t<!(std::is_same_v<DeviceType,LMPDeviceType> || Kokkos::SpaceAccessibility<LMPDeviceType::memory_space,LMPHostType::memory_space>::accessible),const typename kk_view::t_host&> view() const {return h_viewkk;}

  template<class DeviceType>
  std::enable_if_t<(std::is_same<DeviceType,LMPDeviceType>::value || Kokkos::SpaceAccessibility<LMPDeviceType::memory_space,LMPHostType::memory_space>::accessible),void> modify() {modify_device();}

  template<class DeviceType>
  std::enable_if_t<!(std::is_same<DeviceType,LMPDeviceType>::value || Kokkos::SpaceAccessibility<LMPDeviceType::memory_space,LMPHostType::memory_space>::accessible),void> modify() {modify_hostkk();}

  template<class DeviceType>
  std::enable_if_t<(std::is_same<DeviceType,LMPDeviceType>::value || Kokkos::SpaceAccessibility<LMPDeviceType::memory_space,LMPHostType::memory_space>::accessible),void> sync() {sync_device();}

  template<class DeviceType>
  std::enable_if_t<!(std::is_same<DeviceType,LMPDeviceType>::value || Kokkos::SpaceAccessibility<LMPDeviceType::memory_space,LMPHostType::memory_space>::accessible),void> sync() {sync_hostkk();}

  void clear_sync_state()
  {
    k_view.clear_sync_state();

    modified_hostkk_legacy = 0;
    modified_device_legacy = 0;
    modified_legacy_hostkk = 0;
    modified_legacy_device = 0;
  }

  bool need_sync_device()
  {
    return (k_view.need_sync_device() || modified_legacy_device);
  }

  bool need_sync_host()
  {
    return (k_view.need_sync_host() || modified_device_legacy || modified_legacy_hostkk || modified_hostkk_legacy);
  }

  bool need_sync_device_kk()
  {
    return k_view.need_sync_device();
  }

  bool need_sync_host_kk()
  {
    return k_view.need_sync_host();
  }

  KOKKOS_INLINE_FUNCTION
  const legacy_view& view_host() const { return h_view; }

  KOKKOS_INLINE_FUNCTION
  const typename kk_view::t_host& view_hostkk() const { return h_viewkk; }

  KOKKOS_INLINE_FUNCTION
  const typename kk_view::t_dev& view_device() const { return d_view; }

};

// --------------------------------------------------------------------------------

// For device views with fully qualified types
#define KOKKOS_DEVICE_DUALVIEW(TYPE, LAYOUT, SUFFIX) \
typedef Kokkos::DualView<TYPE, LAYOUT, LMPDeviceType> tdual_##SUFFIX; \
typedef typename tdual_##SUFFIX::t_dev t_##SUFFIX; \
typedef typename tdual_##SUFFIX::t_dev_const t_##SUFFIX##_const; \
typedef typename tdual_##SUFFIX::t_dev_um t_##SUFFIX##_um; \
typedef typename tdual_##SUFFIX::t_dev_const_um t_##SUFFIX##_const_um; \
typedef typename tdual_##SUFFIX::t_dev_const_randomread t_##SUFFIX##_randomread;

// For host views with fully qualified types
#define KOKKOS_HOST_DUALVIEW(TYPE, LAYOUT, SUFFIX) \
typedef Kokkos::DualView<TYPE, LAYOUT, LMPDeviceType> tdual_##SUFFIX; \
typedef typename tdual_##SUFFIX::t_host t_##SUFFIX; \
typedef typename tdual_##SUFFIX::t_host_const t_##SUFFIX##_const; \
typedef typename tdual_##SUFFIX::t_host_um t_##SUFFIX##_um; \
typedef typename tdual_##SUFFIX::t_host_const_um t_##SUFFIX##_const_um; \
typedef typename tdual_##SUFFIX::t_host_const_randomread t_##SUFFIX##_randomread;

using LAMMPS_NS::bigint;
using LAMMPS_NS::tagint;
using LAMMPS_NS::imageint;

template <class DeviceType>
struct ArrayTypes;

template <>
struct ArrayTypes<LMPDeviceType> {

// scalar types

KOKKOS_DEVICE_DUALVIEW(int, Kokkos::LayoutRight, int_scalar)
KOKKOS_DEVICE_DUALVIEW(bigint, Kokkos::LayoutRight, bigint_scalar)
KOKKOS_DEVICE_DUALVIEW(tagint, Kokkos::LayoutRight, tagint_scalar)
KOKKOS_DEVICE_DUALVIEW(imageint, Kokkos::LayoutRight, imageint_scalar)
KOKKOS_DEVICE_DUALVIEW(double, Kokkos::LayoutRight, double_scalar)
KOKKOS_DEVICE_DUALVIEW(KK_FLOAT, Kokkos::LayoutRight, kkfloat_scalar)

// 1D view types

KOKKOS_DEVICE_DUALVIEW(int*, Kokkos::LayoutRight, int_1d)
KOKKOS_DEVICE_DUALVIEW(bigint*, Kokkos::LayoutRight, bigint_1d)
KOKKOS_DEVICE_DUALVIEW(tagint*, Kokkos::LayoutRight, tagint_1d)
KOKKOS_DEVICE_DUALVIEW(imageint*, Kokkos::LayoutRight, imageint_1d)
KOKKOS_DEVICE_DUALVIEW(double*, Kokkos::LayoutRight, double_1d)
KOKKOS_DEVICE_DUALVIEW(KK_FLOAT*, Kokkos::LayoutRight, kkfloat_1d)
KOKKOS_DEVICE_DUALVIEW(KK_ACC_FLOAT*, Kokkos::LayoutRight, kkacc_1d)

typedef TransformView<KK_FLOAT*, double*, Kokkos::LayoutRight> ttransform_kkfloat_1d;

// 2D view types

KOKKOS_DEVICE_DUALVIEW(int**, LMPDeviceLayout, int_2d)
KOKKOS_DEVICE_DUALVIEW(int**, Kokkos::LayoutRight, int_2d_lr)
KOKKOS_DEVICE_DUALVIEW(int**, LMPDeviceType::array_layout, int_2d_dl)
KOKKOS_DEVICE_DUALVIEW(int*[3], LMPDeviceLayout, int_1d_3)
KOKKOS_DEVICE_DUALVIEW(tagint**, LMPDeviceLayout, tagint_2d)
KOKKOS_DEVICE_DUALVIEW(double**, Kokkos::LayoutRight, double_2d_lr)
KOKKOS_DEVICE_DUALVIEW(double**, LMPDeviceLayout, double_2d)
KOKKOS_DEVICE_DUALVIEW(double*[2], Kokkos::LayoutRight, double_1d_2_lr)
KOKKOS_DEVICE_DUALVIEW(double*[3], Kokkos::LayoutRight, double_1d_3_lr)
KOKKOS_DEVICE_DUALVIEW(double*[4], Kokkos::LayoutRight, double_1d_4_lr)
KOKKOS_DEVICE_DUALVIEW(double*[6], Kokkos::LayoutRight, double_1d_6_lr)
KOKKOS_DEVICE_DUALVIEW(KK_FLOAT**, LMPDeviceLayout, kkfloat_2d)
KOKKOS_DEVICE_DUALVIEW(KK_FLOAT**, Kokkos::LayoutRight, kkfloat_2d_lr)
KOKKOS_DEVICE_DUALVIEW(KK_FLOAT**, LMPDeviceType::array_layout, kkfloat_2d_dl)
KOKKOS_DEVICE_DUALVIEW(KK_FLOAT*[2], LMPDeviceLayout, kkfloat_1d_2)
KOKKOS_DEVICE_DUALVIEW(KK_FLOAT*[3], LMPDeviceLayout, kkfloat_1d_3)
KOKKOS_DEVICE_DUALVIEW(KK_FLOAT*[3], Kokkos::LayoutRight, kkfloat_1d_3_lr)
KOKKOS_DEVICE_DUALVIEW(KK_ACC_FLOAT*[3], LMPDeviceLayout, kkacc_1d_3)
KOKKOS_DEVICE_DUALVIEW(KK_FLOAT*[4], LMPDeviceLayout, kkfloat_1d_4)
KOKKOS_DEVICE_DUALVIEW(KK_FLOAT*[6], LMPDeviceLayout, kkfloat_1d_6)
KOKKOS_DEVICE_DUALVIEW(KK_ACC_FLOAT*[6], LMPDeviceLayout, kkacc_1d_6)
KOKKOS_DEVICE_DUALVIEW(KK_ACC_FLOAT*[9], LMPDeviceLayout, kkacc_1d_9)

typedef TransformView<KK_ACC_FLOAT*, double*, LMPDeviceLayout> ttransform_kkacc_1d;
typedef TransformView<int**, int**, LMPDeviceLayout> ttransform_int_2d;
typedef TransformView<LAMMPS_NS::tagint**, LAMMPS_NS::tagint**, LMPDeviceLayout> ttransform_tagint_2d;
typedef TransformView<KK_FLOAT**, double**, LMPDeviceLayout> ttransform_kkfloat_2d;
typedef TransformView<KK_FLOAT**, double**, Kokkos::LayoutRight> ttransform_kkfloat_2d_lr;
typedef TransformView<KK_FLOAT*[2], double*[2], LMPDeviceLayout> ttransform_kkfloat_1d_2;
typedef TransformView<KK_FLOAT*[3], double*[3], LMPDeviceLayout> ttransform_kkfloat_1d_3;
typedef TransformView<KK_FLOAT*[3], double*[3], Kokkos::LayoutRight> ttransform_kkfloat_1d_3_lr;
typedef TransformView<KK_ACC_FLOAT*[3], double*[3], LMPDeviceLayout> ttransform_kkacc_1d_3;
typedef TransformView<KK_FLOAT*[4], double*[4], LMPDeviceLayout> ttransform_kkfloat_1d_4;
typedef TransformView<KK_FLOAT*[6], double*[6], LMPDeviceLayout> ttransform_kkfloat_1d_6;
typedef TransformView<KK_ACC_FLOAT*[6], double*[6], LMPDeviceLayout> ttransform_kkacc_1d_6;
typedef TransformView<KK_ACC_FLOAT*[9], double*[9], LMPDeviceLayout> ttransform_kkacc_1d_9;

// 3D view types

KOKKOS_DEVICE_DUALVIEW(int***, LMPDeviceLayout, int_3d)
KOKKOS_DEVICE_DUALVIEW(int***, Kokkos::LayoutRight, int_3d_lr)
KOKKOS_DEVICE_DUALVIEW(double***, Kokkos::LayoutRight, double_3d_lr)
KOKKOS_DEVICE_DUALVIEW(KK_FLOAT***, LMPDeviceLayout, kkfloat_3d)

typedef TransformView<KK_FLOAT***, double***, LMPDeviceLayout> ttransform_kkfloat_3d;

// 4D view types

KOKKOS_DEVICE_DUALVIEW(double****, Kokkos::LayoutRight, double_4d_lr)
KOKKOS_DEVICE_DUALVIEW(KK_FLOAT****, LMPDeviceLayout, kkfloat_4d)

typedef TransformView<KK_FLOAT****, double****, LMPDeviceLayout> ttransform_kkfloat_4d;

// Neighbor Types

typedef tdual_int_2d_dl tdual_neighbors_2d;
typedef tdual_neighbors_2d::t_dev t_neighbors_2d;
typedef tdual_neighbors_2d::t_dev_const t_neighbors_2d_const;
typedef tdual_neighbors_2d::t_dev_um t_neighbors_2d_um;
typedef tdual_neighbors_2d::t_dev_const_um t_neighbors_2d_const_um;
typedef tdual_neighbors_2d::t_dev_const_randomread t_neighbors_2d_randomread;

typedef tdual_int_2d_lr tdual_neighbors_2d_lr;
typedef tdual_neighbors_2d_lr::t_dev t_neighbors_2d_lr;
typedef tdual_neighbors_2d_lr::t_dev_const t_neighbors_2d_const_lr;
typedef tdual_neighbors_2d_lr::t_dev_um t_neighbors_2d_um_lr;
typedef tdual_neighbors_2d_lr::t_dev_const_um t_neighbors_2d_const_um_lr;
typedef tdual_neighbors_2d_lr::t_dev_const_randomread t_neighbors_2d_randomread_lr;

};

#ifdef LMP_KOKKOS_GPU
template <>
struct ArrayTypes<LMPHostType> {

// scalar types

KOKKOS_HOST_DUALVIEW(int, Kokkos::LayoutRight, int_scalar)
KOKKOS_HOST_DUALVIEW(bigint, Kokkos::LayoutRight, bigint_scalar)
KOKKOS_HOST_DUALVIEW(tagint, Kokkos::LayoutRight, tagint_scalar)
KOKKOS_HOST_DUALVIEW(imageint, Kokkos::LayoutRight, imageint_scalar)
KOKKOS_HOST_DUALVIEW(double, Kokkos::LayoutRight, double_scalar)
KOKKOS_HOST_DUALVIEW(KK_FLOAT, Kokkos::LayoutRight, kkfloat_scalar)

// 1D view types

KOKKOS_HOST_DUALVIEW(int*, Kokkos::LayoutRight, int_1d)
KOKKOS_HOST_DUALVIEW(bigint*, Kokkos::LayoutRight, bigint_1d)
KOKKOS_HOST_DUALVIEW(tagint*, Kokkos::LayoutRight, tagint_1d)
KOKKOS_HOST_DUALVIEW(imageint*, Kokkos::LayoutRight, imageint_1d)
KOKKOS_HOST_DUALVIEW(double*, Kokkos::LayoutRight, double_1d)
KOKKOS_HOST_DUALVIEW(KK_FLOAT*, Kokkos::LayoutRight, kkfloat_1d)
KOKKOS_HOST_DUALVIEW(KK_ACC_FLOAT*, Kokkos::LayoutRight, kkacc_1d)

// 2D view types

KOKKOS_HOST_DUALVIEW(int**, LMPDeviceLayout, int_2d)
KOKKOS_HOST_DUALVIEW(int**, Kokkos::LayoutRight, int_2d_lr)
KOKKOS_HOST_DUALVIEW(int**, LMPDeviceType::array_layout, int_2d_dl)
KOKKOS_HOST_DUALVIEW(int*[3], LMPDeviceLayout, int_1d_3)
KOKKOS_HOST_DUALVIEW(tagint**, LMPDeviceLayout, tagint_2d)
KOKKOS_HOST_DUALVIEW(double**, Kokkos::LayoutRight, double_2d_lr)
KOKKOS_HOST_DUALVIEW(double**, LMPDeviceLayout, double_2d)
KOKKOS_HOST_DUALVIEW(double*[2], Kokkos::LayoutRight, double_1d_2_lr)
KOKKOS_HOST_DUALVIEW(double*[3], Kokkos::LayoutRight, double_1d_3_lr)
KOKKOS_HOST_DUALVIEW(double*[4], Kokkos::LayoutRight, double_1d_4_lr)
KOKKOS_HOST_DUALVIEW(double*[6], Kokkos::LayoutRight, double_1d_6_lr)
KOKKOS_HOST_DUALVIEW(KK_FLOAT**, LMPDeviceLayout, kkfloat_2d)
KOKKOS_HOST_DUALVIEW(KK_FLOAT**, Kokkos::LayoutRight, kkfloat_2d_lr)
KOKKOS_HOST_DUALVIEW(KK_FLOAT**, LMPDeviceType::array_layout, kkfloat_2d_dl)
KOKKOS_HOST_DUALVIEW(KK_FLOAT*[2], LMPDeviceLayout, kkfloat_1d_2)
KOKKOS_HOST_DUALVIEW(KK_FLOAT*[3], LMPDeviceLayout, kkfloat_1d_3)
KOKKOS_HOST_DUALVIEW(KK_FLOAT*[3], Kokkos::LayoutRight, kkfloat_1d_3_lr)
KOKKOS_HOST_DUALVIEW(KK_ACC_FLOAT*[3], LMPDeviceLayout, kkacc_1d_3)
KOKKOS_HOST_DUALVIEW(KK_FLOAT*[4], LMPDeviceLayout, kkfloat_1d_4)
KOKKOS_HOST_DUALVIEW(KK_FLOAT*[6], LMPDeviceLayout, kkfloat_1d_6)
KOKKOS_HOST_DUALVIEW(KK_ACC_FLOAT*[6], LMPDeviceLayout, kkacc_1d_6)
KOKKOS_HOST_DUALVIEW(KK_ACC_FLOAT*[9], LMPDeviceLayout, kkacc_1d_9)

// 3D view types

KOKKOS_HOST_DUALVIEW(int***, LMPDeviceLayout, int_3d)
KOKKOS_HOST_DUALVIEW(int***, Kokkos::LayoutRight, int_3d_lr)
KOKKOS_HOST_DUALVIEW(double***, Kokkos::LayoutRight, double_3d_lr)
KOKKOS_HOST_DUALVIEW(KK_FLOAT***, LMPDeviceLayout, kkfloat_3d)

// 4D view types

KOKKOS_HOST_DUALVIEW(double****, Kokkos::LayoutRight, double_4d_lr)
KOKKOS_HOST_DUALVIEW(KK_FLOAT****, LMPDeviceLayout, kkfloat_4d)

// Neighbor Types

typedef tdual_int_2d_dl tdual_neighbors_2d;
typedef tdual_neighbors_2d::t_host t_neighbors_2d;
typedef tdual_neighbors_2d::t_host_const t_neighbors_2d_const;
typedef tdual_neighbors_2d::t_host_um t_neighbors_2d_um;
typedef tdual_neighbors_2d::t_host_const_um t_neighbors_2d_const_um;
typedef tdual_neighbors_2d::t_host_const_randomread t_neighbors_2d_randomread;

typedef tdual_int_2d_lr tdual_neighbors_2d_lr;
typedef tdual_neighbors_2d_lr::t_host t_neighbors_2d_lr;
typedef tdual_neighbors_2d_lr::t_host_const t_neighbors_2d_lr_const;
typedef tdual_neighbors_2d_lr::t_host_um t_neighbors_2d_lr_um;
typedef tdual_neighbors_2d_lr::t_host_const_um t_neighbors_2d_lr_const_um;
typedef tdual_neighbors_2d_lr::t_host_const_randomread t_neighbors_2d_lr_randomread;

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
  Kokkos::parallel_for(view.span()*sizeof(typename ViewType::value_type)/4, f);
  ViewType::execution_space().fence();
}

struct params_lj_coul {
  KOKKOS_INLINE_FUNCTION
  params_lj_coul() {cut_ljsq=0;cut_coulsq=0;lj1=0;lj2=0;lj3=0;lj4=0;offset=0;};
  KOKKOS_INLINE_FUNCTION
  params_lj_coul(int /*i*/) {cut_ljsq=0;cut_coulsq=0;lj1=0;lj2=0;lj3=0;lj4=0;offset=0;};
  KK_FLOAT cut_ljsq,cut_coulsq,lj1,lj2,lj3,lj4,offset;
};

// ReaxFF

struct alignas(4 * sizeof(int)) reax_int4 {
  int i0, i1, i2, i3;
};

// Pair SNAP

#define SNAP_KOKKOS_REAL KK_FLOAT
#define SNAP_KOKKOS_ACCUM KK_ACC_FLOAT
#define SNAP_KOKKOS_HOST_VECLEN 1

#ifdef LMP_KOKKOS_GPU
  #if defined(KOKKOS_ENABLE_SYCL)
    #define SNAP_KOKKOS_DEVICE_VECLEN 16
  #else
    #define SNAP_KOKKOS_DEVICE_VECLEN 32
  #endif
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
#define LAMMPS_CLASS_LAMBDA KOKKOS_CLASS_LAMBDA

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
