// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

module;

#include <Kokkos_Core.hpp>
#include <KokkosExp_InterOp.hpp>

export module kokkos.core;

export {
  namespace Kokkos {
  // hwloc
  namespace hwloc {
  using ::Kokkos::hwloc::available;
  using ::Kokkos::hwloc::bind_this_thread;
  using ::Kokkos::hwloc::can_bind_threads;
  using ::Kokkos::hwloc::get_available_cores_per_numa;
  using ::Kokkos::hwloc::get_available_numa_count;
  using ::Kokkos::hwloc::get_available_threads_per_core;
  using ::Kokkos::hwloc::get_this_thread_coordinate;
  using ::Kokkos::hwloc::thread_mapping;
  using ::Kokkos::hwloc::unbind_this_thread;
  }  // namespace hwloc

  // execution/memory spaces
#ifdef KOKKOS_ENABLE_SERIAL
  using ::Kokkos::NewInstance;
  using ::Kokkos::Serial;
#endif
#ifdef KOKKOS_ENABLE_OPENMP
  using ::Kokkos::OpenMP;
#endif
#ifdef KOKKOS_ENABLE_THREADS
  using ::Kokkos::Threads;
#endif
#ifdef KOKKOS_ENABLE_CUDA
  using ::Kokkos::Cuda;
#endif
#ifdef KOKKOS_ENABLE_HIP
  using ::Kokkos::HIP;
#endif
#ifdef KOKKOS_ENABLE_SYCL
  using ::Kokkos::SYCL;
#endif
  namespace Experimental {
#ifdef KOKKOS_ENABLE_HPX
  using ::Kokkos::Experimental::HPX;
#endif
#ifdef KOKKOS_ENABLE_OPENMPTARGET
  using ::Kokkos::Experimental::OpenMPTarget;
#endif
#ifdef KOKKOS_ENABLE_OPENACC
  using ::Kokkos::Experimental::OpenACC;
#endif
  }  // namespace Experimental
  namespace Experimental {
  using ::Kokkos::Experimental::partition_space;
  }  // namespace Experimental
  using ::Kokkos::AnonymousSpace;
  using ::Kokkos::DefaultExecutionSpace;
  using ::Kokkos::DefaultHostExecutionSpace;
  using ::Kokkos::Device;
  using ::Kokkos::device_id;
  using ::Kokkos::has_shared_host_pinned_space;
  using ::Kokkos::has_shared_space;
  using ::Kokkos::HostSpace;
  using ::Kokkos::is_device;
  using ::Kokkos::is_device_v;
  using ::Kokkos::is_execution_space;
  using ::Kokkos::is_execution_space_v;
  using ::Kokkos::is_memory_space;
  using ::Kokkos::is_memory_space_v;
  using ::Kokkos::is_space;
  using ::Kokkos::ScratchMemorySpace;
  using ::Kokkos::ScratchRequest;
#ifdef KOKKOS_HAS_SHARED_SPACE
  using ::Kokkos::SharedSpace;
#endif
  using ::Kokkos::fence;
  using ::Kokkos::kokkos_free;
  using ::Kokkos::kokkos_malloc;
  using ::Kokkos::kokkos_realloc;
  using ::Kokkos::MemoryPool;  // FIXME
  using ::Kokkos::SpaceAccessibility;
#ifdef KOKKOS_HAS_SHARED_HOST_PINNED_SPACE
  using ::Kokkos::SharedHostPinnedSpace;
#endif

  // View-related
  using ::Kokkos::ALL;
  using ::Kokkos::ALL_t;
  using ::Kokkos::AllowPadding;
  using ::Kokkos::common_view_alloc_prop;
  using ::Kokkos::create_mirror;
  using ::Kokkos::create_mirror_view;
  using ::Kokkos::create_mirror_view_and_copy;
  using ::Kokkos::DeducedCommonPropsType;
  using ::Kokkos::deep_copy;
  using ::Kokkos::InvalidType;
  using ::Kokkos::is_always_assignable;
  using ::Kokkos::is_always_assignable_v;
  using ::Kokkos::is_array_layout;
  using ::Kokkos::is_array_layout_v;
  using ::Kokkos::is_assignable;
  using ::Kokkos::is_memory_traits;
  using ::Kokkos::is_memory_traits_v;
  using ::Kokkos::is_view;
  using ::Kokkos::is_view_v;
  using ::Kokkos::LayoutLeft;
  using ::Kokkos::LayoutRight;
  using ::Kokkos::LayoutStride;
  using ::Kokkos::MemoryRandomAccess;
  using ::Kokkos::MemoryTraits;
  using ::Kokkos::MemoryTraitsFlags;
  using ::Kokkos::MemoryUnmanaged;
  using ::Kokkos::rank;
  using ::Kokkos::Rank;
  using ::Kokkos::realloc;
  using ::Kokkos::resize;
  using ::Kokkos::SequentialHostInit;
  using ::Kokkos::subview;
  using ::Kokkos::Subview;
  using ::Kokkos::Unmanaged;
  using ::Kokkos::View;
  using ::Kokkos::view_alloc;
  using ::Kokkos::view_wrap;
  using ::Kokkos::ViewAllocateWithoutInitializing;
  using ::Kokkos::ViewTraits;
  using ::Kokkos::WithoutInitializing;
  namespace Experimental {
  using ::Kokkos::Experimental::AppendExtent;
  using ::Kokkos::Experimental::Extents;
  using ::Kokkos::Experimental::is_hooks_policy;
  using ::Kokkos::Experimental::is_hooks_policy_v;
  using ::Kokkos::Experimental::local_deep_copy;
  using ::Kokkos::Experimental::local_deep_copy_contiguous;
  using ::Kokkos::Experimental::PrependExtent;
  using ::Kokkos::Experimental::SubscribableViewHooks;
  }  // namespace Experimental

  // execution policies
  using ::Kokkos::AUTO;
  using ::Kokkos::AUTO_t;
  using ::Kokkos::ChunkSize;
  using ::Kokkos::default_inner_direction;
  using ::Kokkos::default_outer_direction;
  using ::Kokkos::Dynamic;
  using ::Kokkos::IndexType;
  using ::Kokkos::is_execution_policy;
  using ::Kokkos::is_execution_policy_v;
  using ::Kokkos::is_team_handle;
  using ::Kokkos::is_team_handle_v;
  using ::Kokkos::Iterate;
  using ::Kokkos::LaunchBounds;
  using ::Kokkos::MDRangePolicy;
  using ::Kokkos::parallel_for;
  using ::Kokkos::parallel_reduce;
  using ::Kokkos::parallel_scan;
  using ::Kokkos::ParallelForTag;
  using ::Kokkos::ParallelReduceTag;
  using ::Kokkos::ParallelScanTag;
  using ::Kokkos::PerTeam;
  using ::Kokkos::PerThread;
  using ::Kokkos::RangePolicy;
  using ::Kokkos::Schedule;
  using ::Kokkos::single;
  using ::Kokkos::Static;
  using ::Kokkos::team_policy_check_valid_storage_level_argument;
  using ::Kokkos::TeamPolicy;
  using ::Kokkos::TeamThreadMDRange;
  using ::Kokkos::TeamThreadRange;
  using ::Kokkos::TeamVectorMDRange;
  using ::Kokkos::TeamVectorRange;
  using ::Kokkos::ThreadVectorMDRange;
  using ::Kokkos::ThreadVectorRange;
  using ::Kokkos::WorkGraphPolicy;
  namespace Experimental {
  using ::Kokkos::Experimental::DesiredOccupancy;
  using ::Kokkos::Experimental::is_work_item_property;
  using ::Kokkos::Experimental::is_work_item_property_v;
  using ::Kokkos::Experimental::MaximizeOccupancy;
  using ::Kokkos::Experimental::partition_space;
  using ::Kokkos::Experimental::prefer;
  using ::Kokkos::Experimental::require;
  using ::Kokkos::Experimental::StaticBatchSize;
  using ::Kokkos::Experimental::WorkItemProperty;
  }  // namespace Experimental

  // miscellaneous
  using ::Kokkos::detected_or_t;
  using ::Kokkos::detected_t;
  using ::Kokkos::is_detected;
  using ::Kokkos::is_detected_convertible;
  using ::Kokkos::is_detected_convertible_v;
  using ::Kokkos::is_detected_exact;
  using ::Kokkos::is_detected_exact_v;
  using ::Kokkos::is_detected_v;
  using ::Kokkos::nonesuch;
  using ::Kokkos::num_devices;
  using ::Kokkos::num_threads;
  using ::Kokkos::print_configuration;
  using ::Kokkos::Timer;
  namespace Experimental {
  using ::Kokkos::Experimental::python_view_type;
  }

  // initialization/finalization
  using ::Kokkos::finalize;
  using ::Kokkos::InitializationSettings;
  using ::Kokkos::initialize;
  using ::Kokkos::is_finalized;
  using ::Kokkos::is_initialized;
  using ::Kokkos::push_finalize_hook;
  using ::Kokkos::ScopeGuard;
  using ::Kokkos::show_warnings;
  using ::Kokkos::tune_internals;

  // std replacements (other than math)
  using ::Kokkos::abort;
  using ::Kokkos::Array;
  using ::Kokkos::begin;
  using ::Kokkos::clamp;
  using ::Kokkos::complex;
  using ::Kokkos::conj;
  using ::Kokkos::end;
  using ::Kokkos::get;
  using ::Kokkos::imag;
  using ::Kokkos::kokkos_swap;
  using ::Kokkos::make_pair;
  using ::Kokkos::max;
  using ::Kokkos::min;
  using ::Kokkos::minmax;
  using ::Kokkos::pair;
  using ::Kokkos::polar;
  using ::Kokkos::printf;
  using ::Kokkos::real;
  using ::Kokkos::tie;
  using ::Kokkos::to_array;

  // reducers
  using ::Kokkos::BAnd;
  using ::Kokkos::BOr;
  using ::Kokkos::FirstLoc;
  using ::Kokkos::FirstLocScalar;
  using ::Kokkos::is_reducer;
  using ::Kokkos::is_reducer_v;
  using ::Kokkos::LAnd;
  using ::Kokkos::LastLoc;
  using ::Kokkos::LastLocScalar;
  using ::Kokkos::LOr;
  using ::Kokkos::Max;
  using ::Kokkos::MaxFirstLoc;
  using ::Kokkos::MaxFirstLocCustomComparator;
  using ::Kokkos::MaxLoc;
  using ::Kokkos::Min;
  using ::Kokkos::MinFirstLoc;
  using ::Kokkos::MinFirstLocCustomComparator;
  using ::Kokkos::MinLoc;
  using ::Kokkos::MinMax;
  using ::Kokkos::MinMaxFirstLastLoc;
  using ::Kokkos::MinMaxFirstLastLocCustomComparator;
  using ::Kokkos::MinMaxLoc;
  using ::Kokkos::MinMaxLocScalar;
  using ::Kokkos::MinMaxScalar;
  using ::Kokkos::Prod;
  using ::Kokkos::reduction_identity;
  using ::Kokkos::StdIsPartitioned;    // FIXME Move to algorithms
  using ::Kokkos::StdIsPartScalar;     // FIXME Move to algorithms
  using ::Kokkos::StdPartitionPoint;   // FIXME Move to algorithms
  using ::Kokkos::StdPartPointScalar;  // FIXME Move to algorithms
  using ::Kokkos::Sum;
  using ::Kokkos::ValLocScalar;

  // half types
  using ::Kokkos::test_fallback_bhalf;  // FIXME
  using ::Kokkos::test_fallback_half;   // FIXME
  namespace Experimental {
  using ::Kokkos::Experimental::bhalf_t;
  using ::Kokkos::Experimental::cast_from_bhalf;
  using ::Kokkos::Experimental::cast_from_half;
  using ::Kokkos::Experimental::cast_to_bhalf;
  using ::Kokkos::Experimental::cast_to_half;
  using ::Kokkos::Experimental::half_t;
  }  // namespace Experimental

  // bit
  namespace Experimental {
  using ::Kokkos::Experimental::bit_cast_builtin;
  using ::Kokkos::Experimental::bit_ceil_builtin;
  using ::Kokkos::Experimental::bit_floor_builtin;
  using ::Kokkos::Experimental::bit_width_builtin;
  using ::Kokkos::Experimental::byteswap_builtin;
  using ::Kokkos::Experimental::countl_one_builtin;
  using ::Kokkos::Experimental::countl_zero_builtin;
  using ::Kokkos::Experimental::countr_one_builtin;
  using ::Kokkos::Experimental::countr_zero_builtin;
  using ::Kokkos::Experimental::has_single_bit_builtin;
  using ::Kokkos::Experimental::popcount_builtin;
  using ::Kokkos::Experimental::rotl_builtin;
  using ::Kokkos::Experimental::rotr_builtin;
  }  // namespace Experimental
  using ::Kokkos::bit_cast;
  using ::Kokkos::bit_ceil;
  using ::Kokkos::bit_floor;
  using ::Kokkos::bit_width;
  using ::Kokkos::byteswap;
  using ::Kokkos::countl_one;
  using ::Kokkos::countl_zero;
  using ::Kokkos::countr_one;
  using ::Kokkos::countr_zero;
  using ::Kokkos::has_single_bit;
  using ::Kokkos::popcount;
  using ::Kokkos::rotl;
  using ::Kokkos::rotr;

  // numeric limits
  namespace Experimental {
  using ::Kokkos::Experimental::denorm_min;
  using ::Kokkos::Experimental::denorm_min_v;
  using ::Kokkos::Experimental::digits;
  using ::Kokkos::Experimental::digits10;
  using ::Kokkos::Experimental::digits10_v;
  using ::Kokkos::Experimental::digits_v;
  using ::Kokkos::Experimental::epsilon;
  using ::Kokkos::Experimental::epsilon_v;
  using ::Kokkos::Experimental::finite_max;
  using ::Kokkos::Experimental::finite_max_v;
  using ::Kokkos::Experimental::finite_min;
  using ::Kokkos::Experimental::finite_min_v;
  using ::Kokkos::Experimental::infinity;
  using ::Kokkos::Experimental::infinity_v;
  using ::Kokkos::Experimental::max_digits10;
  using ::Kokkos::Experimental::max_digits10_v;
  using ::Kokkos::Experimental::max_exponent;
  using ::Kokkos::Experimental::max_exponent10;
  using ::Kokkos::Experimental::max_exponent10_v;
  using ::Kokkos::Experimental::max_exponent_v;
  using ::Kokkos::Experimental::min_exponent;
  using ::Kokkos::Experimental::min_exponent10;
  using ::Kokkos::Experimental::min_exponent10_v;
  using ::Kokkos::Experimental::min_exponent_v;
  using ::Kokkos::Experimental::norm_min;
  using ::Kokkos::Experimental::norm_min_v;
  using ::Kokkos::Experimental::quiet_NaN;
  using ::Kokkos::Experimental::quiet_NaN_v;
  using ::Kokkos::Experimental::radix;
  using ::Kokkos::Experimental::radix_v;
  using ::Kokkos::Experimental::round_error;
  using ::Kokkos::Experimental::round_error_v;
  using ::Kokkos::Experimental::signaling_NaN;
  using ::Kokkos::Experimental::signaling_NaN_v;
  }  // namespace Experimental

  // atomics
  using ::Kokkos::atomic_add;
  using ::Kokkos::atomic_add_fetch;
  using ::Kokkos::atomic_and;
  using ::Kokkos::atomic_and_fetch;
  using ::Kokkos::atomic_compare_exchange;
  using ::Kokkos::atomic_dec;
  using ::Kokkos::atomic_dec_fetch;
  using ::Kokkos::atomic_div;
  using ::Kokkos::atomic_div_fetch;
  using ::Kokkos::atomic_exchange;
  using ::Kokkos::atomic_fetch_add;
  using ::Kokkos::atomic_fetch_and;
  using ::Kokkos::atomic_fetch_dec;
  using ::Kokkos::atomic_fetch_div;
  using ::Kokkos::atomic_fetch_inc;
  using ::Kokkos::atomic_fetch_lshift;
  using ::Kokkos::atomic_fetch_max;
  using ::Kokkos::atomic_fetch_min;
  using ::Kokkos::atomic_fetch_mod;
  using ::Kokkos::atomic_fetch_mul;
  using ::Kokkos::atomic_fetch_nand;
  using ::Kokkos::atomic_fetch_or;
  using ::Kokkos::atomic_fetch_rshift;
  using ::Kokkos::atomic_fetch_sub;
  using ::Kokkos::atomic_fetch_xor;
  using ::Kokkos::atomic_inc;
  using ::Kokkos::atomic_inc_fetch;
  using ::Kokkos::atomic_load;
  using ::Kokkos::atomic_lshift;
  using ::Kokkos::atomic_lshift_fetch;
  using ::Kokkos::atomic_max;
  using ::Kokkos::atomic_max_fetch;
  using ::Kokkos::atomic_min;
  using ::Kokkos::atomic_min_fetch;
  using ::Kokkos::atomic_mod;
  using ::Kokkos::atomic_mod_fetch;
  using ::Kokkos::atomic_mul;
  using ::Kokkos::atomic_mul_fetch;
  using ::Kokkos::atomic_nand;
  using ::Kokkos::atomic_nand_fetch;
  using ::Kokkos::atomic_or;
  using ::Kokkos::atomic_or_fetch;
  using ::Kokkos::atomic_rshift;
  using ::Kokkos::atomic_rshift_fetch;
  using ::Kokkos::atomic_store;
  using ::Kokkos::atomic_sub;
  using ::Kokkos::atomic_sub_fetch;
  using ::Kokkos::atomic_xor;
  using ::Kokkos::atomic_xor_fetch;
  using ::Kokkos::load_fence;
  using ::Kokkos::memory_fence;
  using ::Kokkos::safe_load;
  using ::Kokkos::store_fence;
  using ::Kokkos::volatile_load;
  using ::Kokkos::volatile_store;

  // math functions
  using ::Kokkos::abs;
  using ::Kokkos::acos;
  using ::Kokkos::acosf;
  using ::Kokkos::acosh;
  using ::Kokkos::acoshf;
  using ::Kokkos::acoshl;
  using ::Kokkos::acosl;
  using ::Kokkos::asin;
  using ::Kokkos::asinf;
  using ::Kokkos::asinh;
  using ::Kokkos::asinhf;
  using ::Kokkos::asinhl;
  using ::Kokkos::asinl;
  using ::Kokkos::atan;
  using ::Kokkos::atan2;
  using ::Kokkos::atan2f;
  using ::Kokkos::atan2l;
  using ::Kokkos::atanf;
  using ::Kokkos::atanh;
  using ::Kokkos::atanhf;
  using ::Kokkos::atanhl;
  using ::Kokkos::atanl;
  using ::Kokkos::cbrt;
  using ::Kokkos::cbrtf;
  using ::Kokkos::cbrtl;
  using ::Kokkos::ceil;
  using ::Kokkos::ceilf;
  using ::Kokkos::ceill;
  using ::Kokkos::copysign;
  using ::Kokkos::copysignf;
  using ::Kokkos::copysignl;
  using ::Kokkos::cos;
  using ::Kokkos::cosf;
  using ::Kokkos::cosh;
  using ::Kokkos::coshf;
  using ::Kokkos::coshl;
  using ::Kokkos::cosl;
  using ::Kokkos::erf;
  using ::Kokkos::erfc;
  using ::Kokkos::erfcf;
  using ::Kokkos::erfcl;
  using ::Kokkos::erff;
  using ::Kokkos::erfl;
  using ::Kokkos::exp;
  using ::Kokkos::exp2;
  using ::Kokkos::exp2f;
  using ::Kokkos::exp2l;
  using ::Kokkos::expf;
  using ::Kokkos::expl;
  using ::Kokkos::expm1;
  using ::Kokkos::expm1f;
  using ::Kokkos::expm1l;
  using ::Kokkos::fabs;
  using ::Kokkos::fabsf;
  using ::Kokkos::fabsl;
  using ::Kokkos::fdim;
  using ::Kokkos::fdimf;
  using ::Kokkos::fdiml;
  using ::Kokkos::floor;
  using ::Kokkos::floorf;
  using ::Kokkos::floorl;
  using ::Kokkos::fma;
  using ::Kokkos::fmaf;
  using ::Kokkos::fmal;
  using ::Kokkos::fmax;
  using ::Kokkos::fmaxf;
  using ::Kokkos::fmaxl;
  using ::Kokkos::fmin;
  using ::Kokkos::fminf;
  using ::Kokkos::fminl;
  using ::Kokkos::fmod;
  using ::Kokkos::fmodf;
  using ::Kokkos::fmodl;
  using ::Kokkos::hypot;
  using ::Kokkos::hypotf;
  using ::Kokkos::hypotl;
  using ::Kokkos::isfinite;
  using ::Kokkos::isinf;
  using ::Kokkos::isnan;
  using ::Kokkos::lgamma;
  using ::Kokkos::lgammaf;
  using ::Kokkos::lgammal;
  using ::Kokkos::log;
  using ::Kokkos::log10;
  using ::Kokkos::log10f;
  using ::Kokkos::log10l;
  using ::Kokkos::log1p;
  using ::Kokkos::log1pf;
  using ::Kokkos::log1pl;
  using ::Kokkos::log2;
  using ::Kokkos::log2f;
  using ::Kokkos::log2l;
  using ::Kokkos::logb;
  using ::Kokkos::logbf;
  using ::Kokkos::logbl;
  using ::Kokkos::logf;
  using ::Kokkos::logl;
  using ::Kokkos::nan;
  using ::Kokkos::nanf;
  using ::Kokkos::nanl;
  using ::Kokkos::nearbyint;
  using ::Kokkos::nearbyintf;
  using ::Kokkos::nearbyintl;
  using ::Kokkos::nextafter;
  using ::Kokkos::nextafterf;
  using ::Kokkos::nextafterl;
  using ::Kokkos::pow;
  using ::Kokkos::powf;
  using ::Kokkos::powl;
  using ::Kokkos::remainder;
  using ::Kokkos::remainderf;
  using ::Kokkos::remainderl;
  using ::Kokkos::round;
  using ::Kokkos::roundf;
  using ::Kokkos::roundl;
  using ::Kokkos::rsqrt;
  using ::Kokkos::rsqrtf;
  using ::Kokkos::rsqrtl;
  using ::Kokkos::signbit;
  using ::Kokkos::sin;
  using ::Kokkos::sinf;
  using ::Kokkos::sinh;
  using ::Kokkos::sinhf;
  using ::Kokkos::sinhl;
  using ::Kokkos::sinl;
  using ::Kokkos::sqrt;
  using ::Kokkos::sqrtf;
  using ::Kokkos::sqrtl;
  using ::Kokkos::tan;
  using ::Kokkos::tanf;
  using ::Kokkos::tanh;
  using ::Kokkos::tanhf;
  using ::Kokkos::tanhl;
  using ::Kokkos::tanl;
  using ::Kokkos::tgamma;
  using ::Kokkos::tgammaf;
  using ::Kokkos::tgammal;
  using ::Kokkos::trunc;
  using ::Kokkos::truncf;
  using ::Kokkos::truncl;
  namespace Experimental {
  using ::Kokkos::Experimental::cyl_bessel_h10;
  using ::Kokkos::Experimental::cyl_bessel_h11;
  using ::Kokkos::Experimental::cyl_bessel_h20;
  using ::Kokkos::Experimental::cyl_bessel_h21;
  using ::Kokkos::Experimental::cyl_bessel_i0;
  using ::Kokkos::Experimental::cyl_bessel_i1;
  using ::Kokkos::Experimental::cyl_bessel_j0;
  using ::Kokkos::Experimental::cyl_bessel_j1;
  using ::Kokkos::Experimental::cyl_bessel_k0;
  using ::Kokkos::Experimental::cyl_bessel_k1;
  using ::Kokkos::Experimental::cyl_bessel_y0;
  using ::Kokkos::Experimental::cyl_bessel_y1;
  using ::Kokkos::Experimental::erf;
  using ::Kokkos::Experimental::erfcx;
  using ::Kokkos::Experimental::expint1;
  }  // namespace Experimental

  // numbers
  namespace numbers {
  using ::Kokkos::numbers::e;
  using ::Kokkos::numbers::e_v;
  using ::Kokkos::numbers::egamma;
  using ::Kokkos::numbers::egamma_v;
  using ::Kokkos::numbers::inv_pi;
  using ::Kokkos::numbers::inv_pi_v;
  using ::Kokkos::numbers::inv_sqrt3;
  using ::Kokkos::numbers::inv_sqrt3_v;
  using ::Kokkos::numbers::inv_sqrtpi;
  using ::Kokkos::numbers::inv_sqrtpi_v;
  using ::Kokkos::numbers::ln10;
  using ::Kokkos::numbers::ln10_v;
  using ::Kokkos::numbers::ln2;
  using ::Kokkos::numbers::ln2_v;
  using ::Kokkos::numbers::log10e;
  using ::Kokkos::numbers::log10e_v;
  using ::Kokkos::numbers::log2e;
  using ::Kokkos::numbers::log2e_v;
  using ::Kokkos::numbers::phi;
  using ::Kokkos::numbers::phi_v;
  using ::Kokkos::numbers::pi;
  using ::Kokkos::numbers::pi_v;
  using ::Kokkos::numbers::sqrt2;
  using ::Kokkos::numbers::sqrt2_v;
  using ::Kokkos::numbers::sqrt3;
  using ::Kokkos::numbers::sqrt3_v;
  }  // namespace numbers

  // profiling
  namespace Profiling {
  using ::Kokkos::Profiling::createProfileSection;
  using ::Kokkos::Profiling::destroyProfileSection;
  using ::Kokkos::Profiling::KokkosPDeviceInfo;
  using ::Kokkos::Profiling::markEvent;
  using ::Kokkos::Profiling::popRegion;
  using ::Kokkos::Profiling::pushRegion;
  using ::Kokkos::Profiling::SpaceHandle;
  using ::Kokkos::Profiling::startSection;
  using ::Kokkos::Profiling::stopSection;
  namespace Experimental {
  using ::Kokkos::Profiling::Experimental::set_begin_parallel_for_callback;
  using ::Kokkos::Profiling::Experimental::set_begin_parallel_reduce_callback;
  using ::Kokkos::Profiling::Experimental::set_begin_parallel_scan_callback;
  using ::Kokkos::Profiling::Experimental::set_create_profile_section_callback;
  using ::Kokkos::Profiling::Experimental::set_destroy_profile_section_callback;
  using ::Kokkos::Profiling::Experimental::set_end_deep_copy_callback;
  using ::Kokkos::Profiling::Experimental::set_end_parallel_for_callback;
  using ::Kokkos::Profiling::Experimental::set_end_parallel_reduce_callback;
  using ::Kokkos::Profiling::Experimental::set_end_parallel_scan_callback;
  using ::Kokkos::Profiling::Experimental::set_finalize_callback;
  using ::Kokkos::Profiling::Experimental::set_pop_region_callback;
  using ::Kokkos::Profiling::Experimental::set_profile_event_callback;
  using ::Kokkos::Profiling::Experimental::set_push_region_callback;
  using ::Kokkos::Profiling::Experimental::set_start_profile_section_callback;
  using ::Kokkos::Profiling::Experimental::set_stop_profile_section_callback;
  }  // namespace Experimental
  }  // namespace Profiling

  // tools
  namespace Tools {
  namespace Experimental {

  using ::Kokkos::Tools::Experimental::begin_context;
  using ::Kokkos::Tools::Experimental::CandidateValueType;
  using ::Kokkos::Tools::Experimental::declare_input_type;
  using ::Kokkos::Tools::Experimental::declare_output_type;
  using ::Kokkos::Tools::Experimental::get_new_context_id;
  using ::Kokkos::Tools::Experimental::make_candidate_set;
  using ::Kokkos::Tools::Experimental::make_categorical_tuner;
  using ::Kokkos::Tools::Experimental::make_variable_value;
  using ::Kokkos::Tools::Experimental::request_output_values;
  using ::Kokkos::Tools::Experimental::set_input_values;
  using ::Kokkos::Tools::Experimental::SetOrRange;
  using ::Kokkos::Tools::Experimental::StatisticalCategory;
  using ::Kokkos::Tools::Experimental::ValueType;

  using ::Kokkos::Tools::Experimental::CategoricalTuner;
  using ::Kokkos::Tools::Experimental::contextBeginFunction;
  using ::Kokkos::Tools::Experimental::contextEndFunction;
  using ::Kokkos::Tools::Experimental::declare_optimization_goal;
  using ::Kokkos::Tools::Experimental::device_id;
  using ::Kokkos::Tools::Experimental::device_id_root;
  using ::Kokkos::Tools::Experimental::DeviceType;
  using ::Kokkos::Tools::Experimental::devicetype_from_uint32t;
  using ::Kokkos::Tools::Experimental::end_context;
  using ::Kokkos::Tools::Experimental::EventSet;
  using ::Kokkos::Tools::Experimental::ExecutionSpaceIdentifier;
  using ::Kokkos::Tools::Experimental::ExtendableTunerMixin;
  using ::Kokkos::Tools::Experimental::get_callbacks;
  using ::Kokkos::Tools::Experimental::get_current_context_id;
  using ::Kokkos::Tools::Experimental::get_new_variable_id;
  using ::Kokkos::Tools::Experimental::have_tuning_tool;
  using ::Kokkos::Tools::Experimental::identifier_from_devid;
  using ::Kokkos::Tools::Experimental::inputTypeDeclarationFunction;
  using ::Kokkos::Tools::Experimental::int_for_synchronization_reason;
  using ::Kokkos::Tools::Experimental::make_candidate_range;
  using ::Kokkos::Tools::Experimental::
      make_multidimensional_sparse_tuning_problem;
  using ::Kokkos::Tools::Experimental::MDRangeTuner;
  using ::Kokkos::Tools::Experimental::MultidimensionalSparseTuningProblem;
  using ::Kokkos::Tools::Experimental::num_avail_bits;
  using ::Kokkos::Tools::Experimental::num_device_bits;
  using ::Kokkos::Tools::Experimental::num_instance_bits;
  using ::Kokkos::Tools::Experimental::num_type_bits;
  using ::Kokkos::Tools::Experimental::NumReservedDeviceIDs;
  using ::Kokkos::Tools::Experimental::OptimizationGoal;
  using ::Kokkos::Tools::Experimental::optimizationGoalDeclarationFunction;
  using ::Kokkos::Tools::Experimental::outputTypeDeclarationFunction;
  using ::Kokkos::Tools::Experimental::pause_tools;
  using ::Kokkos::Tools::Experimental::provideToolProgrammingInterfaceFunction;
  using ::Kokkos::Tools::Experimental::RangePolicyOccupancyTuner;
  using ::Kokkos::Tools::Experimental::requestToolSettingsFunction;
  using ::Kokkos::Tools::Experimental::requestValueFunction;
  using ::Kokkos::Tools::Experimental::resume_tools;
  using ::Kokkos::Tools::Experimental::set_allocate_data_callback;
  using ::Kokkos::Tools::Experimental::set_begin_context_callback;
  using ::Kokkos::Tools::Experimental::set_begin_deep_copy_callback;
  using ::Kokkos::Tools::Experimental::set_begin_fence_callback;
  using ::Kokkos::Tools::Experimental::set_begin_parallel_for_callback;
  using ::Kokkos::Tools::Experimental::set_begin_parallel_reduce_callback;
  using ::Kokkos::Tools::Experimental::set_begin_parallel_scan_callback;
  using ::Kokkos::Tools::Experimental::set_callbacks;
  using ::Kokkos::Tools::Experimental::set_create_profile_section_callback;
  using ::Kokkos::Tools::Experimental::set_deallocate_data_callback;
  using ::Kokkos::Tools::Experimental::set_declare_input_type_callback;
  using ::Kokkos::Tools::Experimental::set_declare_metadata_callback;
  using ::Kokkos::Tools::Experimental::set_declare_optimization_goal_callback;
  using ::Kokkos::Tools::Experimental::set_declare_output_type_callback;
  using ::Kokkos::Tools::Experimental::set_destroy_profile_section_callback;
  using ::Kokkos::Tools::Experimental::set_dual_view_modify_callback;
  using ::Kokkos::Tools::Experimental::set_dual_view_sync_callback;
  using ::Kokkos::Tools::Experimental::set_end_context_callback;
  using ::Kokkos::Tools::Experimental::set_end_deep_copy_callback;
  using ::Kokkos::Tools::Experimental::set_end_fence_callback;
  using ::Kokkos::Tools::Experimental::set_end_parallel_for_callback;
  using ::Kokkos::Tools::Experimental::set_end_parallel_reduce_callback;
  using ::Kokkos::Tools::Experimental::set_end_parallel_scan_callback;
  using ::Kokkos::Tools::Experimental::set_finalize_callback;
  using ::Kokkos::Tools::Experimental::set_init_callback;
  using ::Kokkos::Tools::Experimental::set_parse_args_callback;
  using ::Kokkos::Tools::Experimental::set_pop_region_callback;
  using ::Kokkos::Tools::Experimental::set_print_help_callback;
  using ::Kokkos::Tools::Experimental::set_profile_event_callback;
  using ::Kokkos::Tools::Experimental::
      set_provide_tool_programming_interface_callback;
  using ::Kokkos::Tools::Experimental::set_push_region_callback;
  using ::Kokkos::Tools::Experimental::set_request_output_values_callback;
  using ::Kokkos::Tools::Experimental::set_request_tool_settings_callback;
  using ::Kokkos::Tools::Experimental::set_start_profile_section_callback;
  using ::Kokkos::Tools::Experimental::set_stop_profile_section_callback;
  using ::Kokkos::Tools::Experimental::SingleDimensionalRangeTuner;
  using ::Kokkos::Tools::Experimental::SpecialSynchronizationCases;
  using ::Kokkos::Tools::Experimental::TeamSizeTuner;
  using ::Kokkos::Tools::Experimental::toolInvokedFenceFunction;
  using ::Kokkos::Tools::Experimental::ToolProgrammingInterface;
  using ::Kokkos::Tools::Experimental::ToolSettings;
  using ::Kokkos::Tools::Experimental::TuningString;
  using ::Kokkos::Tools::Experimental::ValueRange;
  using ::Kokkos::Tools::Experimental::ValueSet;
  using ::Kokkos::Tools::Experimental::VariableInfo;
  using ::Kokkos::Tools::Experimental::VariableValue;
  }  // namespace Experimental
  using ::Kokkos::Tools::allocateData;
  using ::Kokkos::Tools::allocateDataFunction;
  using ::Kokkos::Tools::beginDeepCopy;
  using ::Kokkos::Tools::beginDeepCopyFunction;
  using ::Kokkos::Tools::beginFence;
  using ::Kokkos::Tools::beginFenceFunction;
  using ::Kokkos::Tools::beginFunction;
  using ::Kokkos::Tools::beginParallelFor;
  using ::Kokkos::Tools::beginParallelReduce;
  using ::Kokkos::Tools::beginParallelScan;
  using ::Kokkos::Tools::createProfileSection;
  using ::Kokkos::Tools::createProfileSectionFunction;
  using ::Kokkos::Tools::deallocateData;
  using ::Kokkos::Tools::deallocateDataFunction;
  using ::Kokkos::Tools::declareMetadata;
  using ::Kokkos::Tools::declareMetadataFunction;
  using ::Kokkos::Tools::destroyProfileSection;
  using ::Kokkos::Tools::destroyProfileSectionFunction;
  using ::Kokkos::Tools::dualViewModifyFunction;
  using ::Kokkos::Tools::dualViewSyncFunction;
  using ::Kokkos::Tools::endDeepCopy;
  using ::Kokkos::Tools::endDeepCopyFunction;
  using ::Kokkos::Tools::endFence;
  using ::Kokkos::Tools::endFenceFunction;
  using ::Kokkos::Tools::endFunction;
  using ::Kokkos::Tools::endParallelFor;
  using ::Kokkos::Tools::endParallelReduce;
  using ::Kokkos::Tools::endParallelScan;
  using ::Kokkos::Tools::finalize;
  using ::Kokkos::Tools::finalizeFunction;
  using ::Kokkos::Tools::InitArguments;
  using ::Kokkos::Tools::initFunction;
  using ::Kokkos::Tools::initialize;
  using ::Kokkos::Tools::make_space_handle;
  using ::Kokkos::Tools::markEvent;
  using ::Kokkos::Tools::modifyDualView;
  using ::Kokkos::Tools::parseArgs;
  using ::Kokkos::Tools::parseArgsFunction;
  using ::Kokkos::Tools::popFunction;
  using ::Kokkos::Tools::popRegion;
  using ::Kokkos::Tools::printHelp;
  using ::Kokkos::Tools::printHelpFunction;
  using ::Kokkos::Tools::profileEventFunction;
  using ::Kokkos::Tools::profileLibraryLoaded;
  using ::Kokkos::Tools::pushFunction;
  using ::Kokkos::Tools::pushRegion;
  using ::Kokkos::Tools::SpaceHandle;
  using ::Kokkos::Tools::startProfileSectionFunction;
  using ::Kokkos::Tools::startSection;
  using ::Kokkos::Tools::stopProfileSectionFunction;
  using ::Kokkos::Tools::stopSection;
  using ::Kokkos::Tools::syncDualView;
  }  // namespace Tools

  // Crs
  using ::Kokkos::count_and_fill_crs;
  using ::Kokkos::CountAndFill;
  using ::Kokkos::CountAndFillBase;
  using ::Kokkos::Crs;
  using ::Kokkos::get_crs_row_map_from_counts;
  using ::Kokkos::get_crs_transpose_counts;
  using ::Kokkos::transpose_crs;

  // mdspan
  using ::Kokkos::default_accessor;
  using ::Kokkos::dextents;
  using ::Kokkos::dynamic_extent;
  using ::Kokkos::extents;
  using ::Kokkos::full_extent;
  using ::Kokkos::full_extent_t;
  using ::Kokkos::layout_left;
  using ::Kokkos::layout_right;
  using ::Kokkos::layout_stride;
  using ::Kokkos::mdspan;
  using ::Kokkos::mdspan_non_standard;
  using ::Kokkos::mdspan_non_standard_tag;
  using ::Kokkos::strided_slice;
  using ::Kokkos::submdspan;
  using ::Kokkos::submdspan_extents;
  using ::Kokkos::submdspan_mapping_result;
  namespace Experimental {
  using ::Kokkos::Experimental::dims;
  using ::Kokkos::Experimental::layout_left_padded;
  using ::Kokkos::Experimental::layout_right_padded;
  }  // namespace Experimental

  // UniqueToken
  namespace Experimental {
  using ::Kokkos::Experimental::AcquireTeamUniqueToken;
  using ::Kokkos::Experimental::AcquireUniqueToken;
  using ::Kokkos::Experimental::UniqueToken;
  using ::Kokkos::Experimental::UniqueTokenScope;
  }  // namespace Experimental

  // operators
  using ::Kokkos::operator+;
  using ::Kokkos::operator-;
  using ::Kokkos::operator*;
  using ::Kokkos::operator/;
  using ::Kokkos::operator<<;
  using ::Kokkos::operator>>;
  using ::Kokkos::operator==;
  using ::Kokkos::operator!=;
  using ::Kokkos::operator<;
  using ::Kokkos::operator<=;
  using ::Kokkos::operator>;
  using ::Kokkos::operator>=;
  namespace Experimental {
  using ::Kokkos::Experimental::operator&;
  using ::Kokkos::Experimental::operator==;
  using ::Kokkos::Experimental::operator|;
  }  // namespace Experimental

  }  // namespace Kokkos
}
