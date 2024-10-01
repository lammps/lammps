//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

/// \file Kokkos_DualView.hpp
/// \brief Declaration and definition of Kokkos::DualView.
///
/// This header file declares and defines Kokkos::DualView and its
/// related nonmember functions.

#ifndef KOKKOS_DUALVIEW_HPP
#define KOKKOS_DUALVIEW_HPP
#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_DUALVIEW
#endif

#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Error.hpp>

namespace Kokkos {

/* \class DualView
 * \brief Container to manage mirroring a Kokkos::View that lives
 *   in device memory with a Kokkos::View that lives in host memory.
 *
 * This class provides capabilities to manage data which exists in two
 * memory spaces at the same time.  It keeps views of the same layout
 * on two memory spaces as well as modified flags for both
 * allocations.  Users are responsible for setting the modified flags
 * manually if they change the data in either memory space, by calling
 * the sync() method templated on the device where they modified the
 * data.  Users may synchronize data by calling the modify() function,
 * templated on the device towards which they want to synchronize
 * (i.e., the target of the one-way copy operation).
 *
 * The DualView class also provides convenience methods such as
 * realloc, resize and capacity which call the appropriate methods of
 * the underlying Kokkos::View objects.
 *
 * The four template arguments are the same as those of Kokkos::View.
 * (Please refer to that class' documentation for a detailed
 * description.)
 *
 *   \tparam DataType The type of the entries stored in the container.
 *
 *   \tparam Layout The array's layout in memory.
 *
 *   \tparam Device The Kokkos Device type.  If its memory space is
 *     not the same as the host's memory space, then DualView will
 *     contain two separate Views: one in device memory, and one in
 *     host memory.  Otherwise, DualView will only store one View.
 *
 *   \tparam MemoryTraits (optional) The user's intended memory access
 *     behavior.  Please see the documentation of Kokkos::View for
 *     examples.  The default suffices for most users.
 */

namespace Impl {

#ifdef KOKKOS_ENABLE_CUDA

inline const Kokkos::Cuda& get_cuda_space(const Kokkos::Cuda& in) { return in; }

inline const Kokkos::Cuda& get_cuda_space() {
  return *Kokkos::Impl::cuda_get_deep_copy_space();
}

template <typename NonCudaExecSpace>
inline const Kokkos::Cuda& get_cuda_space(const NonCudaExecSpace&) {
  return get_cuda_space();
}

#endif  // KOKKOS_ENABLE_CUDA

}  // namespace Impl

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
template <class DataType, class Arg1Type = void, class Arg2Type = void,
          class Arg3Type = void>
class DualView;
#else
template <class DataType, class... Properties>
class DualView;
#endif

template <class>
struct is_dual_view : public std::false_type {};

template <class DT, class... DP>
struct is_dual_view<DualView<DT, DP...>> : public std::true_type {};

template <class DT, class... DP>
struct is_dual_view<const DualView<DT, DP...>> : public std::true_type {};

template <class T>
inline constexpr bool is_dual_view_v = is_dual_view<T>::value;

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
template <class DataType, class Arg1Type, class Arg2Type, class Arg3Type>
class DualView : public ViewTraits<DataType, Arg1Type, Arg2Type, Arg3Type> {
  template <class, class, class, class>
#else
template <class DataType, class... Properties>
class DualView : public ViewTraits<DataType, Properties...> {
  template <class, class...>
#endif
  friend class DualView;

 public:
  //! \name Typedefs for device types and various Kokkos::View specializations.
  //@{
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
  using traits = ViewTraits<DataType, Arg1Type, Arg2Type, Arg3Type>;
#else
  using traits      = ViewTraits<DataType, Properties...>;
#endif

  //! The Kokkos Host Device type;
  using host_mirror_space = typename traits::host_mirror_space;

  //! The type of a Kokkos::View on the device.
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
  using t_dev = View<typename traits::data_type, Arg1Type, Arg2Type, Arg3Type>;
#else
  using t_dev       = View<typename traits::data_type, Properties...>;
#endif

  /// \typedef t_host
  /// \brief The type of a Kokkos::View host mirror of \c t_dev.
  using t_host = typename t_dev::HostMirror;

  //! The type of a const View on the device.
  //! The type of a Kokkos::View on the device.
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
  using t_dev_const =
      View<typename traits::const_data_type, Arg1Type, Arg2Type, Arg3Type>;
#else
  using t_dev_const = View<typename traits::const_data_type, Properties...>;
#endif

  /// \typedef t_host_const
  /// \brief The type of a const View host mirror of \c t_dev_const.
  using t_host_const = typename t_dev_const::HostMirror;

  //! The type of a const, random-access View on the device.
  using t_dev_const_randomread =
      View<typename traits::const_data_type, typename traits::array_layout,
           typename traits::device_type,
           Kokkos::MemoryTraits<Kokkos::RandomAccess>>;

  /// \typedef t_host_const_randomread
  /// \brief The type of a const, random-access View host mirror of
  ///   \c t_dev_const_randomread.
  using t_host_const_randomread = typename t_dev_const_randomread::HostMirror;

  //! The type of an unmanaged View on the device.
  using t_dev_um =
      View<typename traits::data_type, typename traits::array_layout,
           typename traits::device_type, MemoryUnmanaged>;

  //! The type of an unmanaged View host mirror of \c t_dev_um.
  using t_host_um =
      View<typename t_host::data_type, typename t_host::array_layout,
           typename t_host::device_type, MemoryUnmanaged>;

  //! The type of a const unmanaged View on the device.
  using t_dev_const_um =
      View<typename traits::const_data_type, typename traits::array_layout,
           typename traits::device_type, MemoryUnmanaged>;

  //! The type of a const unmanaged View host mirror of \c t_dev_const_um.
  using t_host_const_um =
      View<typename t_host::const_data_type, typename t_host::array_layout,
           typename t_host::device_type, MemoryUnmanaged>;

  //! The type of a const, random-access View on the device.
  using t_dev_const_randomread_um =
      View<typename t_host::const_data_type, typename t_host::array_layout,
           typename t_host::device_type,
           Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>;

  /// \typedef t_host_const_randomread
  /// \brief The type of a const, random-access View host mirror of
  ///   \c t_dev_const_randomread.
  using t_host_const_randomread_um =
      typename t_dev_const_randomread_um::HostMirror;

  //@}
  //! \name Counters to keep track of changes ("modified" flags)
  //@{

 protected:
  // modified_flags[0] -> host
  // modified_flags[1] -> device
  using t_modified_flags = View<unsigned int[2], LayoutLeft, Kokkos::HostSpace>;
  t_modified_flags modified_flags;

 public:
  //@}

  // Moved this specifically after modified_flags to resolve an alignment issue
  // on MSVC/NVCC
  //! \name The two View instances.
  //@{
  t_dev d_view;
  t_host h_view;
  //@}

  //! \name Constructors
  //@{

  /// \brief Empty constructor.
  ///
  /// Both device and host View objects are constructed using their
  /// default constructors.  The "modified" flags are both initialized
  /// to "unmodified."
  DualView() = default;

  /// \brief Constructor that allocates View objects on both host and device.
  ///
  /// This constructor works like the analogous constructor of View.
  /// The first argument is a string label, which is entirely for your
  /// benefit.  (Different DualView objects may have the same label if
  /// you like.)  The arguments that follow are the dimensions of the
  /// View objects.  For example, if the View has three dimensions,
  /// the first three integer arguments will be nonzero, and you may
  /// omit the integer arguments that follow.
  DualView(const std::string& label,
           const size_t n0 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
           const size_t n1 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
           const size_t n2 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
           const size_t n3 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
           const size_t n4 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
           const size_t n5 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
           const size_t n6 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
           const size_t n7 = KOKKOS_IMPL_CTOR_DEFAULT_ARG)
      : modified_flags(
            Kokkos::view_alloc(typename t_modified_flags::execution_space{},
                               "DualView::modified_flags")),
        d_view(label, n0, n1, n2, n3, n4, n5, n6, n7),
        h_view(create_mirror_view(d_view))  // without UVM, host View mirrors
  {}

  /// \brief Constructor that allocates View objects on both host and device.
  ///
  /// This constructor works like the analogous constructor of View.
  /// The first arguments are wrapped up in a ViewCtor class, this allows
  /// for a label, without initializing, and all of the other things that can
  /// be wrapped up in a Ctor class.
  /// The arguments that follow are the dimensions of the
  /// View objects.  For example, if the View has three dimensions,
  /// the first three integer arguments will be nonzero, and you may
  /// omit the integer arguments that follow.
  template <class... P>
  DualView(const Impl::ViewCtorProp<P...>& arg_prop,
           std::enable_if_t<!Impl::ViewCtorProp<P...>::has_pointer,
                            size_t> const n0 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
           const size_t n1                   = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
           const size_t n2                   = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
           const size_t n3                   = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
           const size_t n4                   = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
           const size_t n5                   = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
           const size_t n6                   = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
           const size_t n7                   = KOKKOS_IMPL_CTOR_DEFAULT_ARG)
      : modified_flags(t_modified_flags("DualView::modified_flags")),
        d_view(arg_prop, n0, n1, n2, n3, n4, n5, n6, n7) {
    // without UVM, host View mirrors
    if constexpr (Kokkos::Impl::has_type<Impl::WithoutInitializing_t,
                                         P...>::value)
      h_view = Kokkos::create_mirror_view(Kokkos::WithoutInitializing, d_view);
    else
      h_view = Kokkos::create_mirror_view(d_view);
  }

  //! Copy constructor (shallow copy)
  template <typename DT, typename... DP>
  DualView(const DualView<DT, DP...>& src)
      : modified_flags(src.modified_flags),
        d_view(src.d_view),
        h_view(src.h_view) {}

  //! Subview constructor
  template <class DT, class... DP, class Arg0, class... Args>
  DualView(const DualView<DT, DP...>& src, const Arg0& arg0, Args... args)
      : modified_flags(src.modified_flags),
        d_view(Kokkos::subview(src.d_view, arg0, args...)),
        h_view(Kokkos::subview(src.h_view, arg0, args...)) {}

  /// \brief Create DualView from existing device and host View objects.
  ///
  /// This constructor assumes that the device and host View objects
  /// are synchronized.  You, the caller, are responsible for making
  /// sure this is the case before calling this constructor.  After
  /// this constructor returns, you may use DualView's sync() and
  /// modify() methods to ensure synchronization of the View objects.
  ///
  /// \param d_view_ Device View
  /// \param h_view_ Host View (must have type t_host = t_dev::HostMirror)
  DualView(const t_dev& d_view_, const t_host& h_view_)
      : modified_flags(t_modified_flags("DualView::modified_flags")),
        d_view(d_view_),
        h_view(h_view_) {
    if (int(d_view.rank) != int(h_view.rank) ||
        d_view.extent(0) != h_view.extent(0) ||
        d_view.extent(1) != h_view.extent(1) ||
        d_view.extent(2) != h_view.extent(2) ||
        d_view.extent(3) != h_view.extent(3) ||
        d_view.extent(4) != h_view.extent(4) ||
        d_view.extent(5) != h_view.extent(5) ||
        d_view.extent(6) != h_view.extent(6) ||
        d_view.extent(7) != h_view.extent(7) ||
        d_view.stride_0() != h_view.stride_0() ||
        d_view.stride_1() != h_view.stride_1() ||
        d_view.stride_2() != h_view.stride_2() ||
        d_view.stride_3() != h_view.stride_3() ||
        d_view.stride_4() != h_view.stride_4() ||
        d_view.stride_5() != h_view.stride_5() ||
        d_view.stride_6() != h_view.stride_6() ||
        d_view.stride_7() != h_view.stride_7() ||
        d_view.span() != h_view.span()) {
      Kokkos::Impl::throw_runtime_exception(
          "DualView constructed with incompatible views");
    }
  }
  // does the DualView have only one device
  struct impl_dualview_is_single_device {
    enum : bool {
      value = std::is_same<typename t_dev::device_type,
                           typename t_host::device_type>::value
    };
  };

  // does the given device match the device of t_dev?
  template <typename Device>
  struct impl_device_matches_tdev_device {
    enum : bool {
      value = std::is_same<typename t_dev::device_type, Device>::value
    };
  };
  // does the given device match the device of t_host?
  template <typename Device>
  struct impl_device_matches_thost_device {
    enum : bool {
      value = std::is_same<typename t_host::device_type, Device>::value
    };
  };

  // does the given device match the execution space of t_host?
  template <typename Device>
  struct impl_device_matches_thost_exec {
    enum : bool {
      value = std::is_same<typename t_host::execution_space, Device>::value
    };
  };

  // does the given device match the execution space of t_dev?
  template <typename Device>
  struct impl_device_matches_tdev_exec {
    enum : bool {
      value = std::is_same<typename t_dev::execution_space, Device>::value
    };
  };

  // does the given device's memory space match the memory space of t_dev?
  template <typename Device>
  struct impl_device_matches_tdev_memory_space {
    enum : bool {
      value = std::is_same<typename t_dev::memory_space,
                           typename Device::memory_space>::value
    };
  };

  //@}
  //! \name Methods for synchronizing, marking as modified, and getting Views.
  //@{

  /// \brief Return a View on a specific device \c Device.
  ///
  /// Please don't be afraid of the nested if_c expressions in the return
  /// value's type.  That just tells the method what the return type
  /// should be: t_dev if the \c Device template parameter matches
  /// this DualView's device type, else t_host.
  ///
  /// For example, suppose you create a DualView on Cuda, like this:
  /// \code
  ///   using dual_view_type =
  ///       Kokkos::DualView<float, Kokkos::LayoutRight, Kokkos::Cuda>;
  ///   dual_view_type DV ("my dual view", 100);
  /// \endcode
  /// If you want to get the CUDA device View, do this:
  /// \code
  ///   typename dual_view_type::t_dev cudaView = DV.view<Kokkos::Cuda> ();
  /// \endcode
  /// and if you want to get the host mirror of that View, do this:
  /// \code
  ///   using host_device_type = typename Kokkos::HostSpace::execution_space;
  ///   typename dual_view_type::t_host hostView = DV.view<host_device_type> ();
  /// \endcode
  template <class Device>
  KOKKOS_INLINE_FUNCTION const typename std::conditional_t<
      impl_device_matches_tdev_device<Device>::value, t_dev,
      typename std::conditional_t<
          impl_device_matches_thost_device<Device>::value, t_host,
          typename std::conditional_t<
              impl_device_matches_thost_exec<Device>::value, t_host,
              typename std::conditional_t<
                  impl_device_matches_tdev_exec<Device>::value, t_dev,
                  typename std::conditional_t<
                      impl_device_matches_tdev_memory_space<Device>::value,
                      t_dev, t_host>>>>>
  view() const {
    constexpr bool device_is_memspace =
        std::is_same<Device, typename Device::memory_space>::value;
    constexpr bool device_is_execspace =
        std::is_same<Device, typename Device::execution_space>::value;
    constexpr bool device_exec_is_t_dev_exec =
        std::is_same<typename Device::execution_space,
                     typename t_dev::execution_space>::value;
    constexpr bool device_mem_is_t_dev_mem =
        std::is_same<typename Device::memory_space,
                     typename t_dev::memory_space>::value;
    constexpr bool device_exec_is_t_host_exec =
        std::is_same<typename Device::execution_space,
                     typename t_host::execution_space>::value;
    constexpr bool device_mem_is_t_host_mem =
        std::is_same<typename Device::memory_space,
                     typename t_host::memory_space>::value;
    constexpr bool device_is_t_host_device =
        std::is_same<typename Device::execution_space,
                     typename t_host::device_type>::value;
    constexpr bool device_is_t_dev_device =
        std::is_same<typename Device::memory_space,
                     typename t_host::device_type>::value;

    static_assert(
        device_is_t_dev_device || device_is_t_host_device ||
            (device_is_memspace &&
             (device_mem_is_t_dev_mem || device_mem_is_t_host_mem)) ||
            (device_is_execspace &&
             (device_exec_is_t_dev_exec || device_exec_is_t_host_exec)) ||
            ((!device_is_execspace && !device_is_memspace) &&
             ((device_mem_is_t_dev_mem || device_mem_is_t_host_mem) ||
              (device_exec_is_t_dev_exec || device_exec_is_t_host_exec))),
        "Template parameter to .view() must exactly match one of the "
        "DualView's device types or one of the execution or memory spaces");

    return Impl::if_c<std::is_same<typename t_dev::memory_space,
                                   typename Device::memory_space>::value,
                      t_dev, t_host>::select(d_view, h_view);
  }

  KOKKOS_INLINE_FUNCTION
  t_host view_host() const { return h_view; }

  KOKKOS_INLINE_FUNCTION
  t_dev view_device() const { return d_view; }

  KOKKOS_INLINE_FUNCTION constexpr bool is_allocated() const {
    return (d_view.is_allocated() && h_view.is_allocated());
  }

  template <class Device>
  static int get_device_side() {
    constexpr bool device_is_memspace =
        std::is_same<Device, typename Device::memory_space>::value;
    constexpr bool device_is_execspace =
        std::is_same<Device, typename Device::execution_space>::value;
    constexpr bool device_exec_is_t_dev_exec =
        std::is_same<typename Device::execution_space,
                     typename t_dev::execution_space>::value;
    constexpr bool device_mem_is_t_dev_mem =
        std::is_same<typename Device::memory_space,
                     typename t_dev::memory_space>::value;
    constexpr bool device_exec_is_t_host_exec =
        std::is_same<typename Device::execution_space,
                     typename t_host::execution_space>::value;
    constexpr bool device_mem_is_t_host_mem =
        std::is_same<typename Device::memory_space,
                     typename t_host::memory_space>::value;
    constexpr bool device_is_t_host_device =
        std::is_same<typename Device::execution_space,
                     typename t_host::device_type>::value;
    constexpr bool device_is_t_dev_device =
        std::is_same<typename Device::memory_space,
                     typename t_host::device_type>::value;

    static_assert(
        device_is_t_dev_device || device_is_t_host_device ||
            (device_is_memspace &&
             (device_mem_is_t_dev_mem || device_mem_is_t_host_mem)) ||
            (device_is_execspace &&
             (device_exec_is_t_dev_exec || device_exec_is_t_host_exec)) ||
            ((!device_is_execspace && !device_is_memspace) &&
             ((device_mem_is_t_dev_mem || device_mem_is_t_host_mem) ||
              (device_exec_is_t_dev_exec || device_exec_is_t_host_exec))),
        "Template parameter to .sync() must exactly match one of the "
        "DualView's device types or one of the execution or memory spaces");

    int dev = -1;
    if (device_is_t_dev_device)
      dev = 1;
    else if (device_is_t_host_device)
      dev = 0;
    else {
      if (device_is_memspace) {
        if (device_mem_is_t_dev_mem) dev = 1;
        if (device_mem_is_t_host_mem) dev = 0;
        if (device_mem_is_t_host_mem && device_mem_is_t_dev_mem) dev = -1;
      }
      if (device_is_execspace) {
        if (device_exec_is_t_dev_exec) dev = 1;
        if (device_exec_is_t_host_exec) dev = 0;
        if (device_exec_is_t_host_exec && device_exec_is_t_dev_exec) dev = -1;
      }
      if (!device_is_execspace && !device_is_memspace) {
        if (device_mem_is_t_dev_mem) dev = 1;
        if (device_mem_is_t_host_mem) dev = 0;
        if (device_mem_is_t_host_mem && device_mem_is_t_dev_mem) dev = -1;
        if (device_exec_is_t_dev_exec) dev = 1;
        if (device_exec_is_t_host_exec) dev = 0;
        if (device_exec_is_t_host_exec && device_exec_is_t_dev_exec) dev = -1;
      }
    }
    return dev;
  }
  static constexpr const int view_header_size = 128;
  void impl_report_host_sync() const noexcept {
    if (Kokkos::Tools::Experimental::get_callbacks().sync_dual_view !=
        nullptr) {
      Kokkos::Tools::syncDualView(
          h_view.label(),
          reinterpret_cast<void*>(reinterpret_cast<uintptr_t>(h_view.data()) -
                                  view_header_size),
          false);
    }
  }
  void impl_report_device_sync() const noexcept {
    if (Kokkos::Tools::Experimental::get_callbacks().sync_dual_view !=
        nullptr) {
      Kokkos::Tools::syncDualView(
          d_view.label(),
          reinterpret_cast<void*>(reinterpret_cast<uintptr_t>(d_view.data()) -
                                  view_header_size),
          true);
    }
  }

  /// \brief Update data on device or host only if data in the other
  ///   space has been marked as modified.
  ///
  /// If \c Device is the same as this DualView's device type, then
  /// copy data from host to device.  Otherwise, copy data from device
  /// to host.  In either case, only copy if the source of the copy
  /// has been modified.
  ///
  /// This is a one-way synchronization only.  If the target of the
  /// copy has been modified, this operation will discard those
  /// modifications.  It will also reset both device and host modified
  /// flags.
  ///
  /// \note This method doesn't know on its own whether you modified
  ///   the data in either View.  You must manually mark modified data
  ///   as modified, by calling the modify() method with the
  ///   appropriate template parameter.
  // deliberately passing args by cref as they're used multiple times
  template <class Device, class... Args>
  void sync_impl(std::true_type, Args const&... args) {
    if (modified_flags.data() == nullptr) return;

    int dev = get_device_side<Device>();

    if (dev == 1) {  // if Device is the same as DualView's device type
      if ((modified_flags(0) > 0) && (modified_flags(0) >= modified_flags(1))) {
#ifdef KOKKOS_ENABLE_CUDA
        if (std::is_same<typename t_dev::memory_space,
                         Kokkos::CudaUVMSpace>::value) {
          if (d_view.data() == h_view.data())
            Kokkos::Impl::cuda_prefetch_pointer(
                Impl::get_cuda_space(args...), d_view.data(),
                sizeof(typename t_dev::value_type) * d_view.span(), true);
        }
#endif

        deep_copy(args..., d_view, h_view);
        modified_flags(0) = modified_flags(1) = 0;
        impl_report_device_sync();
      }
    }
    if (dev == 0) {  // hopefully Device is the same as DualView's host type
      if ((modified_flags(1) > 0) && (modified_flags(1) >= modified_flags(0))) {
#ifdef KOKKOS_ENABLE_CUDA
        if (std::is_same<typename t_dev::memory_space,
                         Kokkos::CudaUVMSpace>::value) {
          if (d_view.data() == h_view.data())
            Kokkos::Impl::cuda_prefetch_pointer(
                Impl::get_cuda_space(args...), d_view.data(),
                sizeof(typename t_dev::value_type) * d_view.span(), false);
        }
#endif

        deep_copy(args..., h_view, d_view);
        modified_flags(0) = modified_flags(1) = 0;
        impl_report_host_sync();
      }
    }
    if constexpr (std::is_same<typename t_host::memory_space,
                               typename t_dev::memory_space>::value) {
      typename t_dev::execution_space().fence(
          "Kokkos::DualView<>::sync: fence after syncing DualView");
      typename t_host::execution_space().fence(
          "Kokkos::DualView<>::sync: fence after syncing DualView");
    }
  }

  template <class Device>
  void sync(const std::enable_if_t<
                (std::is_same<typename traits::data_type,
                              typename traits::non_const_data_type>::value) ||
                    (std::is_same<Device, int>::value),
                int>& = 0) {
    sync_impl<Device>(std::true_type{});
  }

  template <class Device, class ExecutionSpace>
  void sync(const ExecutionSpace& exec,
            const std::enable_if_t<
                (std::is_same<typename traits::data_type,
                              typename traits::non_const_data_type>::value) ||
                    (std::is_same<Device, int>::value),
                int>& = 0) {
    sync_impl<Device>(std::true_type{}, exec);
  }

  // deliberately passing args by cref as they're used multiple times
  template <class Device, class... Args>
  void sync_impl(std::false_type, Args const&...) {
    if (modified_flags.data() == nullptr) return;

    int dev = get_device_side<Device>();

    if (dev == 1) {  // if Device is the same as DualView's device type
      if ((modified_flags(0) > 0) && (modified_flags(0) >= modified_flags(1))) {
        Impl::throw_runtime_exception(
            "Calling sync on a DualView with a const datatype.");
      }
      impl_report_device_sync();
    }
    if (dev == 0) {  // hopefully Device is the same as DualView's host type
      if ((modified_flags(1) > 0) && (modified_flags(1) >= modified_flags(0))) {
        Impl::throw_runtime_exception(
            "Calling sync on a DualView with a const datatype.");
      }
      impl_report_host_sync();
    }
  }

  template <class Device>
  void sync(const std::enable_if_t<
                (!std::is_same<typename traits::data_type,
                               typename traits::non_const_data_type>::value) ||
                    (std::is_same<Device, int>::value),
                int>& = 0) {
    sync_impl<Device>(std::false_type{});
  }
  template <class Device, class ExecutionSpace>
  void sync(const ExecutionSpace& exec,
            const std::enable_if_t<
                (!std::is_same<typename traits::data_type,
                               typename traits::non_const_data_type>::value) ||
                    (std::is_same<Device, int>::value),
                int>& = 0) {
    sync_impl<Device>(std::false_type{}, exec);
  }

  // deliberately passing args by cref as they're used multiple times
  template <typename... Args>
  void sync_host_impl(Args const&... args) {
    if (!std::is_same<typename traits::data_type,
                      typename traits::non_const_data_type>::value)
      Impl::throw_runtime_exception(
          "Calling sync_host on a DualView with a const datatype.");
    if (modified_flags.data() == nullptr) return;
    if (modified_flags(1) > modified_flags(0)) {
#ifdef KOKKOS_ENABLE_CUDA
      if (std::is_same<typename t_dev::memory_space,
                       Kokkos::CudaUVMSpace>::value) {
        if (d_view.data() == h_view.data())
          Kokkos::Impl::cuda_prefetch_pointer(
              Impl::get_cuda_space(args...), d_view.data(),
              sizeof(typename t_dev::value_type) * d_view.span(), false);
      }
#endif

      deep_copy(args..., h_view, d_view);
      modified_flags(1) = modified_flags(0) = 0;
      impl_report_host_sync();
    }
  }

  template <class ExecSpace>
  void sync_host(const ExecSpace& exec) {
    sync_host_impl(exec);
  }
  void sync_host() { sync_host_impl(); }

  // deliberately passing args by cref as they're used multiple times
  template <typename... Args>
  void sync_device_impl(Args const&... args) {
    if (!std::is_same<typename traits::data_type,
                      typename traits::non_const_data_type>::value)
      Impl::throw_runtime_exception(
          "Calling sync_device on a DualView with a const datatype.");
    if (modified_flags.data() == nullptr) return;
    if (modified_flags(0) > modified_flags(1)) {
#ifdef KOKKOS_ENABLE_CUDA
      if (std::is_same<typename t_dev::memory_space,
                       Kokkos::CudaUVMSpace>::value) {
        if (d_view.data() == h_view.data())
          Kokkos::Impl::cuda_prefetch_pointer(
              Impl::get_cuda_space(args...), d_view.data(),
              sizeof(typename t_dev::value_type) * d_view.span(), true);
      }
#endif

      deep_copy(args..., d_view, h_view);
      modified_flags(1) = modified_flags(0) = 0;
      impl_report_device_sync();
    }
  }

  template <class ExecSpace>
  void sync_device(const ExecSpace& exec) {
    sync_device_impl(exec);
  }
  void sync_device() { sync_device_impl(); }

  template <class Device>
  bool need_sync() const {
    if (modified_flags.data() == nullptr) return false;
    int dev = get_device_side<Device>();

    if (dev == 1) {  // if Device is the same as DualView's device type
      if ((modified_flags(0) > 0) && (modified_flags(0) >= modified_flags(1))) {
        return true;
      }
    }
    if (dev == 0) {  // hopefully Device is the same as DualView's host type
      if ((modified_flags(1) > 0) && (modified_flags(1) >= modified_flags(0))) {
        return true;
      }
    }
    return false;
  }

  inline bool need_sync_host() const {
    if (modified_flags.data() == nullptr) return false;
    return modified_flags(0) < modified_flags(1);
  }

  inline bool need_sync_device() const {
    if (modified_flags.data() == nullptr) return false;
    return modified_flags(1) < modified_flags(0);
  }
  void impl_report_device_modification() {
    if (Kokkos::Tools::Experimental::get_callbacks().modify_dual_view !=
        nullptr) {
      Kokkos::Tools::modifyDualView(
          d_view.label(),
          reinterpret_cast<void*>(reinterpret_cast<uintptr_t>(d_view.data()) -
                                  view_header_size),
          true);
    }
  }
  void impl_report_host_modification() {
    if (Kokkos::Tools::Experimental::get_callbacks().modify_dual_view !=
        nullptr) {
      Kokkos::Tools::modifyDualView(
          h_view.label(),
          reinterpret_cast<void*>(reinterpret_cast<uintptr_t>(h_view.data()) -
                                  view_header_size),
          false);
    }
  }
  /// \brief Mark data as modified on the given device \c Device.
  ///
  /// If \c Device is the same as this DualView's device type, then
  /// mark the device's data as modified.  Otherwise, mark the host's
  /// data as modified.
  template <class Device, class Dummy = DualView,
            std::enable_if_t<!Dummy::impl_dualview_is_single_device::value>* =
                nullptr>
  void modify() {
    if (modified_flags.data() == nullptr) {
      modified_flags = t_modified_flags("DualView::modified_flags");
    }

    int dev = get_device_side<Device>();

    if (dev == 1) {  // if Device is the same as DualView's device type
      // Increment the device's modified count.
      modified_flags(1) =
          (modified_flags(1) > modified_flags(0) ? modified_flags(1)
                                                 : modified_flags(0)) +
          1;
      impl_report_device_modification();
    }
    if (dev == 0) {  // hopefully Device is the same as DualView's host type
      // Increment the host's modified count.
      modified_flags(0) =
          (modified_flags(1) > modified_flags(0) ? modified_flags(1)
                                                 : modified_flags(0)) +
          1;
      impl_report_host_modification();
    }

#ifdef KOKKOS_ENABLE_DEBUG_DUALVIEW_MODIFY_CHECK
    if (modified_flags(0) && modified_flags(1)) {
      std::string msg = "Kokkos::DualView::modify ERROR: ";
      msg += "Concurrent modification of host and device views ";
      msg += "in DualView \"";
      msg += d_view.label();
      msg += "\"\n";
      Kokkos::abort(msg.c_str());
    }
#endif
  }

  template <
      class Device, class Dummy = DualView,
      std::enable_if_t<Dummy::impl_dualview_is_single_device::value>* = nullptr>
  void modify() {
    return;
  }

  template <class Dummy = DualView,
            std::enable_if_t<!Dummy::impl_dualview_is_single_device::value>* =
                nullptr>
  inline void modify_host() {
    if (modified_flags.data() != nullptr) {
      modified_flags(0) =
          (modified_flags(1) > modified_flags(0) ? modified_flags(1)
                                                 : modified_flags(0)) +
          1;
      impl_report_host_modification();
#ifdef KOKKOS_ENABLE_DEBUG_DUALVIEW_MODIFY_CHECK
      if (modified_flags(0) && modified_flags(1)) {
        std::string msg = "Kokkos::DualView::modify_host ERROR: ";
        msg += "Concurrent modification of host and device views ";
        msg += "in DualView \"";
        msg += d_view.label();
        msg += "\"\n";
        Kokkos::abort(msg.c_str());
      }
#endif
    }
  }

  template <
      class Dummy = DualView,
      std::enable_if_t<Dummy::impl_dualview_is_single_device::value>* = nullptr>
  inline void modify_host() {
    return;
  }

  template <class Dummy = DualView,
            std::enable_if_t<!Dummy::impl_dualview_is_single_device::value>* =
                nullptr>
  inline void modify_device() {
    if (modified_flags.data() != nullptr) {
      modified_flags(1) =
          (modified_flags(1) > modified_flags(0) ? modified_flags(1)
                                                 : modified_flags(0)) +
          1;
      impl_report_device_modification();
#ifdef KOKKOS_ENABLE_DEBUG_DUALVIEW_MODIFY_CHECK
      if (modified_flags(0) && modified_flags(1)) {
        std::string msg = "Kokkos::DualView::modify_device ERROR: ";
        msg += "Concurrent modification of host and device views ";
        msg += "in DualView \"";
        msg += d_view.label();
        msg += "\"\n";
        Kokkos::abort(msg.c_str());
      }
#endif
    }
  }

  template <
      class Dummy = DualView,
      std::enable_if_t<Dummy::impl_dualview_is_single_device::value>* = nullptr>
  inline void modify_device() {
    return;
  }

  inline void clear_sync_state() {
    if (modified_flags.data() != nullptr)
      modified_flags(1) = modified_flags(0) = 0;
  }

  //@}
  //! \name Methods for reallocating or resizing the View objects.
  //@{

  /// \brief Reallocate both View objects.
  ///
  /// This discards any existing contents of the objects, and resets
  /// their modified flags.  It does <i>not</i> copy the old contents
  /// of either View into the new View objects.
  template <class... ViewCtorArgs>
  void impl_realloc(const size_t n0, const size_t n1, const size_t n2,
                    const size_t n3, const size_t n4, const size_t n5,
                    const size_t n6, const size_t n7,
                    const Impl::ViewCtorProp<ViewCtorArgs...>& arg_prop) {
    using alloc_prop_input = Impl::ViewCtorProp<ViewCtorArgs...>;

    static_assert(!alloc_prop_input::has_label,
                  "The view constructor arguments passed to Kokkos::realloc "
                  "must not include a label!");
    static_assert(
        !alloc_prop_input::has_pointer,
        "The view constructor arguments passed to Kokkos::realloc must "
        "not include a pointer!");
    static_assert(
        !alloc_prop_input::has_memory_space,
        "The view constructor arguments passed to Kokkos::realloc must "
        "not include a memory space instance!");

    const size_t new_extents[8] = {n0, n1, n2, n3, n4, n5, n6, n7};
    const bool sizeMismatch =
        Impl::size_mismatch(h_view, h_view.rank_dynamic, new_extents);

    if (sizeMismatch) {
      ::Kokkos::realloc(arg_prop, d_view, n0, n1, n2, n3, n4, n5, n6, n7);
      if constexpr (alloc_prop_input::initialize) {
        h_view = create_mirror_view(typename t_host::memory_space(), d_view);
      } else {
        h_view = create_mirror_view(Kokkos::WithoutInitializing,
                                    typename t_host::memory_space(), d_view);
      }
    } else if constexpr (alloc_prop_input::initialize) {
      if constexpr (alloc_prop_input::has_execution_space) {
        const auto& exec_space =
            Impl::get_property<Impl::ExecutionSpaceTag>(arg_prop);
        ::Kokkos::deep_copy(exec_space, d_view, typename t_dev::value_type{});
      } else
        ::Kokkos::deep_copy(d_view, typename t_dev::value_type{});
    }

    /* Reset dirty flags */
    if (modified_flags.data() == nullptr) {
      modified_flags = t_modified_flags("DualView::modified_flags");
    } else
      modified_flags(1) = modified_flags(0) = 0;
  }

  template <class... ViewCtorArgs>
  void realloc(const Impl::ViewCtorProp<ViewCtorArgs...>& arg_prop,
               const size_t n0 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
               const size_t n1 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
               const size_t n2 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
               const size_t n3 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
               const size_t n4 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
               const size_t n5 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
               const size_t n6 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
               const size_t n7 = KOKKOS_IMPL_CTOR_DEFAULT_ARG) {
    impl_realloc(n0, n1, n2, n3, n4, n5, n6, n7, arg_prop);
  }

  void realloc(const size_t n0 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
               const size_t n1 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
               const size_t n2 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
               const size_t n3 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
               const size_t n4 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
               const size_t n5 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
               const size_t n6 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
               const size_t n7 = KOKKOS_IMPL_CTOR_DEFAULT_ARG) {
    impl_realloc(n0, n1, n2, n3, n4, n5, n6, n7, Impl::ViewCtorProp<>{});
  }

  template <typename I>
  std::enable_if_t<Impl::is_view_ctor_property<I>::value> realloc(
      const I& arg_prop, const size_t n0 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
      const size_t n1 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
      const size_t n2 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
      const size_t n3 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
      const size_t n4 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
      const size_t n5 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
      const size_t n6 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
      const size_t n7 = KOKKOS_IMPL_CTOR_DEFAULT_ARG) {
    impl_realloc(n0, n1, n2, n3, n4, n5, n6, n7, Kokkos::view_alloc(arg_prop));
  }

  /// \brief Resize both views, copying old contents into new if necessary.
  ///
  /// This method only copies the old contents into the new View
  /// objects for the device which was last marked as modified.
  template <class... ViewCtorArgs>
  void impl_resize(const Impl::ViewCtorProp<ViewCtorArgs...>& arg_prop,
                   const size_t n0, const size_t n1, const size_t n2,
                   const size_t n3, const size_t n4, const size_t n5,
                   const size_t n6, const size_t n7) {
    using alloc_prop_input = Impl::ViewCtorProp<ViewCtorArgs...>;

    static_assert(!alloc_prop_input::has_label,
                  "The view constructor arguments passed to Kokkos::resize "
                  "must not include a label!");
    static_assert(
        !alloc_prop_input::has_pointer,
        "The view constructor arguments passed to Kokkos::resize must "
        "not include a pointer!");
    static_assert(
        !alloc_prop_input::has_memory_space,
        "The view constructor arguments passed to Kokkos::resize must "
        "not include a memory space instance!");

    const size_t new_extents[8] = {n0, n1, n2, n3, n4, n5, n6, n7};
    const bool sizeMismatch =
        Impl::size_mismatch(h_view, h_view.rank_dynamic, new_extents);

    if (modified_flags.data() == nullptr) {
      modified_flags = t_modified_flags("DualView::modified_flags");
    }

    [[maybe_unused]] auto resize_on_device = [&](const auto& properties) {
      /* Resize on Device */
      if (sizeMismatch) {
        ::Kokkos::resize(properties, d_view, n0, n1, n2, n3, n4, n5, n6, n7);
        // this part of the lambda was relocated in a method as it contains a
        // `if constexpr`. In some cases, both branches were evaluated
        // leading to a compile error
        resync_host(properties);

        /* Mark Device copy as modified */
        ++modified_flags(1);
      }
    };

    [[maybe_unused]] auto resize_on_host = [&](const auto& properties) {
      /* Resize on Host */
      if (sizeMismatch) {
        ::Kokkos::resize(properties, h_view, n0, n1, n2, n3, n4, n5, n6, n7);
        // this part of the lambda was relocated in a method as it contains a
        // `if constexpr`. In some cases, both branches were evaluated
        // leading to a compile error
        resync_device(properties);

        /* Mark Host copy as modified */
        ++modified_flags(0);
      }
    };

    constexpr bool has_execution_space = alloc_prop_input::has_execution_space;

    if constexpr (has_execution_space) {
      using ExecSpace = typename alloc_prop_input::execution_space;
      const auto& exec_space =
          Impl::get_property<Impl::ExecutionSpaceTag>(arg_prop);
      constexpr bool exec_space_can_access_device =
          SpaceAccessibility<ExecSpace,
                             typename t_dev::memory_space>::accessible;
      constexpr bool exec_space_can_access_host =
          SpaceAccessibility<ExecSpace,
                             typename t_host::memory_space>::accessible;
      static_assert(exec_space_can_access_device || exec_space_can_access_host);
      if constexpr (exec_space_can_access_device) {
        sync<typename t_dev::memory_space>(exec_space);
        resize_on_device(arg_prop);
        return;
      }
      if constexpr (exec_space_can_access_host) {
        sync<typename t_host::memory_space>(exec_space);
        resize_on_host(arg_prop);
        return;
      }
    } else {
      if (modified_flags(1) >= modified_flags(0)) {
        resize_on_device(arg_prop);
      } else {
        resize_on_host(arg_prop);
      }
    }
  }

 private:
  // resync host mirror from device
  // this code was relocated from a lambda as it contains a `if constexpr`.
  // In some cases, both branches were evaluated, leading to a compile error
  template <class... ViewCtorArgs>
  inline void resync_host(Impl::ViewCtorProp<ViewCtorArgs...> const&) {
    using alloc_prop_input = Impl::ViewCtorProp<ViewCtorArgs...>;

    if constexpr (alloc_prop_input::initialize) {
      h_view = create_mirror_view(typename t_host::memory_space(), d_view);
    } else {
      h_view = create_mirror_view(Kokkos::WithoutInitializing,
                                  typename t_host::memory_space(), d_view);
    }
  }

  // resync device mirror from host
  // this code was relocated from a lambda as it contains a `if constexpr`
  // In some cases, both branches were evaluated leading to a compile error
  template <class... ViewCtorArgs>
  inline void resync_device(Impl::ViewCtorProp<ViewCtorArgs...> const&) {
    using alloc_prop_input = Impl::ViewCtorProp<ViewCtorArgs...>;

    if constexpr (alloc_prop_input::initialize) {
      d_view = create_mirror_view(typename t_dev::memory_space(), h_view);

    } else {
      d_view = create_mirror_view(Kokkos::WithoutInitializing,
                                  typename t_dev::memory_space(), h_view);
    }
  }

 public:
  void resize(const size_t n0 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
              const size_t n1 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
              const size_t n2 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
              const size_t n3 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
              const size_t n4 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
              const size_t n5 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
              const size_t n6 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
              const size_t n7 = KOKKOS_IMPL_CTOR_DEFAULT_ARG) {
    impl_resize(Impl::ViewCtorProp<>{}, n0, n1, n2, n3, n4, n5, n6, n7);
  }

  template <class... ViewCtorArgs>
  void resize(const Impl::ViewCtorProp<ViewCtorArgs...>& arg_prop,
              const size_t n0 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
              const size_t n1 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
              const size_t n2 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
              const size_t n3 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
              const size_t n4 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
              const size_t n5 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
              const size_t n6 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
              const size_t n7 = KOKKOS_IMPL_CTOR_DEFAULT_ARG) {
    impl_resize(arg_prop, n0, n1, n2, n3, n4, n5, n6, n7);
  }

  template <class I>
  std::enable_if_t<Impl::is_view_ctor_property<I>::value> resize(
      const I& arg_prop, const size_t n0 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
      const size_t n1 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
      const size_t n2 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
      const size_t n3 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
      const size_t n4 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
      const size_t n5 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
      const size_t n6 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
      const size_t n7 = KOKKOS_IMPL_CTOR_DEFAULT_ARG) {
    impl_resize(Kokkos::view_alloc(arg_prop), n0, n1, n2, n3, n4, n5, n6, n7);
  }

  //@}
  //! \name Methods for getting capacity, stride, or dimension(s).
  //@{

  //! The allocation size (same as Kokkos::View::span).
  KOKKOS_INLINE_FUNCTION constexpr size_t span() const { return d_view.span(); }

  KOKKOS_INLINE_FUNCTION bool span_is_contiguous() const {
    return d_view.span_is_contiguous();
  }

  //! Get stride(s) for each dimension.
  template <typename iType>
  void stride(iType* stride_) const {
    d_view.stride(stride_);
  }

  template <typename iType>
  KOKKOS_INLINE_FUNCTION constexpr std::enable_if_t<
      std::is_integral<iType>::value, size_t>
  extent(const iType& r) const {
    return d_view.extent(r);
  }

  template <typename iType>
  KOKKOS_INLINE_FUNCTION constexpr std::enable_if_t<
      std::is_integral<iType>::value, int>
  extent_int(const iType& r) const {
    return static_cast<int>(d_view.extent(r));
  }

  //@}
};

}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//
// Partial specializations of Kokkos::subview() for DualView objects.
//

namespace Kokkos {
namespace Impl {

template <class V>
struct V2DV;

template <class D, class... P>
struct V2DV<View<D, P...>> {
  using type = DualView<D, P...>;
};
} /* namespace Impl */

template <class DataType, class... Properties, class... Args>
auto subview(const DualView<DataType, Properties...>& src, Args&&... args) {
  // leverage Kokkos::View facilities to deduce the properties of the subview
  using deduce_subview_type =
      decltype(subview(std::declval<View<DataType, Properties...>>(),
                       std::forward<Args>(args)...));
  // map it back to dual view
  return typename Impl::V2DV<deduce_subview_type>::type(
      src, std::forward<Args>(args)...);
}

} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

//
// Partial specialization of Kokkos::deep_copy() for DualView objects.
//

template <class DT, class... DP, class ST, class... SP>
void deep_copy(DualView<DT, DP...>& dst, const DualView<ST, SP...>& src) {
  if (src.need_sync_device()) {
    deep_copy(dst.h_view, src.h_view);
    dst.modify_host();
  } else {
    deep_copy(dst.d_view, src.d_view);
    dst.modify_device();
  }
}

template <class ExecutionSpace, class DT, class... DP, class ST, class... SP>
void deep_copy(const ExecutionSpace& exec, DualView<DT, DP...>& dst,
               const DualView<ST, SP...>& src) {
  if (src.need_sync_device()) {
    deep_copy(exec, dst.h_view, src.h_view);
    dst.modify_host();
  } else {
    deep_copy(exec, dst.d_view, src.d_view);
    dst.modify_device();
  }
}

}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

//
// Non-member resize and realloc
//

template <class... Properties, class... Args>
void resize(DualView<Properties...>& dv, Args&&... args) noexcept(
    noexcept(dv.resize(std::forward<Args>(args)...))) {
  dv.resize(std::forward<Args>(args)...);
}

template <class... ViewCtorArgs, class... Properties, class... Args>
void resize(
    const Impl::ViewCtorProp<ViewCtorArgs...>& arg_prop,
    DualView<Properties...>& dv,
    Args&&... args) noexcept(noexcept(dv.resize(arg_prop,
                                                std::forward<Args>(args)...))) {
  dv.resize(arg_prop, std::forward<Args>(args)...);
}

template <class I, class... Properties, class... Args>
std::enable_if_t<Impl::is_view_ctor_property<I>::value> resize(
    const I& arg_prop, DualView<Properties...>& dv,
    Args&&... args) noexcept(noexcept(dv.resize(arg_prop,
                                                std::forward<Args>(args)...))) {
  dv.resize(arg_prop, std::forward<Args>(args)...);
}

template <class... ViewCtorArgs, class... Properties, class... Args>
void realloc(const Impl::ViewCtorProp<ViewCtorArgs...>& arg_prop,
             DualView<Properties...>& dv,
             Args&&... args) noexcept(noexcept(dv
                                                   .realloc(std::forward<Args>(
                                                       args)...))) {
  dv.realloc(arg_prop, std::forward<Args>(args)...);
}

template <class... Properties, class... Args>
void realloc(DualView<Properties...>& dv, Args&&... args) noexcept(
    noexcept(dv.realloc(std::forward<Args>(args)...))) {
  dv.realloc(std::forward<Args>(args)...);
}

template <class I, class... Properties, class... Args>
std::enable_if_t<Impl::is_view_ctor_property<I>::value> realloc(
    const I& arg_prop, DualView<Properties...>& dv,
    Args&&... args) noexcept(noexcept(dv.realloc(arg_prop,
                                                 std::forward<Args>(
                                                     args)...))) {
  dv.realloc(arg_prop, std::forward<Args>(args)...);
}

}  // end namespace Kokkos

#ifdef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_DUALVIEW
#undef KOKKOS_IMPL_PUBLIC_INCLUDE
#undef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_DUALVIEW
#endif
#endif
