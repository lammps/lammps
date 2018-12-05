/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 2.0
//              Copyright (2014) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Christian R. Trott (crtrott@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

/// \file Kokkos_DualView.hpp
/// \brief Declaration and definition of Kokkos::DualView.
///
/// This header file declares and defines Kokkos::DualView and its
/// related nonmember functions.

#ifndef KOKKOS_DUALVIEW_HPP
#define KOKKOS_DUALVIEW_HPP

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
template< class DataType ,
          class Arg1Type = void ,
          class Arg2Type = void ,
          class Arg3Type = void>
class DualView : public ViewTraits< DataType , Arg1Type , Arg2Type, Arg3Type >
{
template< class , class , class , class > friend class DualView ;
public:
  //! \name Typedefs for device types and various Kokkos::View specializations.
  //@{
  typedef ViewTraits< DataType , Arg1Type , Arg2Type, Arg3Type > traits ;

  //! The Kokkos Host Device type;
  typedef typename traits::host_mirror_space host_mirror_space ;

  //! The type of a Kokkos::View on the device.
  typedef View< typename traits::data_type ,
                Arg1Type ,
                Arg2Type ,
                Arg3Type > t_dev ;

  /// \typedef t_host
  /// \brief The type of a Kokkos::View host mirror of \c t_dev.
  typedef typename t_dev::HostMirror t_host ;

  //! The type of a const View on the device.
  //! The type of a Kokkos::View on the device.
  typedef View< typename traits::const_data_type ,
                Arg1Type ,
                Arg2Type ,
                Arg3Type > t_dev_const ;

  /// \typedef t_host_const
  /// \brief The type of a const View host mirror of \c t_dev_const.
  typedef typename t_dev_const::HostMirror t_host_const;

  //! The type of a const, random-access View on the device.
  typedef View< typename traits::const_data_type ,
                typename traits::array_layout ,
                typename traits::device_type ,
                Kokkos::MemoryTraits<Kokkos::RandomAccess> > t_dev_const_randomread ;

  /// \typedef t_host_const_randomread
  /// \brief The type of a const, random-access View host mirror of
  ///   \c t_dev_const_randomread.
  typedef typename t_dev_const_randomread::HostMirror t_host_const_randomread;

  //! The type of an unmanaged View on the device.
  typedef View< typename traits::data_type ,
                typename traits::array_layout ,
                typename traits::device_type ,
                MemoryUnmanaged> t_dev_um;

  //! The type of an unmanaged View host mirror of \c t_dev_um.
  typedef View< typename t_host::data_type ,
                typename t_host::array_layout ,
                typename t_host::device_type ,
                MemoryUnmanaged> t_host_um;

  //! The type of a const unmanaged View on the device.
  typedef View< typename traits::const_data_type ,
                typename traits::array_layout ,
                typename traits::device_type ,
                MemoryUnmanaged> t_dev_const_um;

  //! The type of a const unmanaged View host mirror of \c t_dev_const_um.
  typedef View<typename t_host::const_data_type,
               typename t_host::array_layout,
               typename t_host::device_type,
               MemoryUnmanaged> t_host_const_um;

  //! The type of a const, random-access View on the device.
  typedef View< typename t_host::const_data_type ,
                typename t_host::array_layout ,
                typename t_host::device_type ,
                Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess> > t_dev_const_randomread_um ;

  /// \typedef t_host_const_randomread
  /// \brief The type of a const, random-access View host mirror of
  ///   \c t_dev_const_randomread.
  typedef typename t_dev_const_randomread::HostMirror t_host_const_randomread_um;

  //@}
  //! \name The two View instances.
  //@{

  t_dev d_view;
  t_host h_view;

  //@}
  //! \name Counters to keep track of changes ("modified" flags)
  //@{

#ifndef KOKKOS_ENABLE_DEPRECATED_CODE
protected:
  // modified_flags[0] -> host
  // modified_flags[1] -> device
  typedef View<unsigned int[2],LayoutLeft,Kokkos::HostSpace> t_modified_flags;
  t_modified_flags modified_flags;

public:
#else
  typedef View<unsigned int[2],LayoutLeft,typename t_host::execution_space> t_modified_flags;
  typedef View<unsigned int,LayoutLeft,typename t_host::execution_space> t_modified_flag;
  t_modified_flags modified_flags;
  t_modified_flag modified_host,modified_device;
#endif

  //@}
  //! \name Constructors
  //@{

  /// \brief Empty constructor.
  ///
  /// Both device and host View objects are constructed using their
  /// default constructors.  The "modified" flags are both initialized
  /// to "unmodified."
#ifndef KOKKOS_ENABLE_DEPRECATED_CODE
  DualView () = default;
#else
  DualView ():modified_flags (t_modified_flags("DualView::modified_flags")) {
    modified_host = t_modified_flag(modified_flags,0);
    modified_device = t_modified_flag(modified_flags,1);
  }
#endif

  /// \brief Constructor that allocates View objects on both host and device.
  ///
  /// This constructor works like the analogous constructor of View.
  /// The first argument is a string label, which is entirely for your
  /// benefit.  (Different DualView objects may have the same label if
  /// you like.)  The arguments that follow are the dimensions of the
  /// View objects.  For example, if the View has three dimensions,
  /// the first three integer arguments will be nonzero, and you may
  /// omit the integer arguments that follow.
  DualView (const std::string& label,
            const size_t n0 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
            const size_t n1 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
            const size_t n2 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
            const size_t n3 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
            const size_t n4 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
            const size_t n5 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
            const size_t n6 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
            const size_t n7 = KOKKOS_IMPL_CTOR_DEFAULT_ARG)
    : d_view (label, n0, n1, n2, n3, n4, n5, n6, n7)
    , h_view (create_mirror_view (d_view)) // without UVM, host View mirrors
    , modified_flags (t_modified_flags("DualView::modified_flags"))
  {
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
    modified_host = t_modified_flag(modified_flags,0);
    modified_device = t_modified_flag(modified_flags,1);
#endif
  }

  //! Copy constructor (shallow copy)
  template<class SS, class LS, class DS, class MS>
  DualView (const DualView<SS,LS,DS,MS>& src) :
    d_view (src.d_view),
    h_view (src.h_view),
    modified_flags (src.modified_flags)
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
    , modified_host(src.modified_host)
    , modified_device(src.modified_device)
#endif
  {}

  //! Subview constructor
  template< class SD, class S1 , class S2 , class S3
          , class Arg0 , class ... Args >
  DualView( const DualView<SD,S1,S2,S3> & src
          , const Arg0 & arg0
          , Args ... args
          )
    : d_view( Kokkos::subview( src.d_view , arg0 , args ... ) )
    , h_view( Kokkos::subview( src.h_view , arg0 , args ... ) )
    , modified_flags (src.modified_flags)
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
    , modified_host(src.modified_host)
    , modified_device(src.modified_device)
#endif
    {}

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
  DualView (const t_dev& d_view_, const t_host& h_view_) :
    d_view (d_view_),
    h_view (h_view_),
    modified_flags (t_modified_flags("DualView::modified_flags"))
  {
    if ( int(d_view.rank)     != int(h_view.rank) ||
         d_view.extent(0) != h_view.extent(0) ||
         d_view.extent(1) != h_view.extent(1) ||
         d_view.extent(2) != h_view.extent(2) ||
         d_view.extent(3) != h_view.extent(3) ||
         d_view.extent(4) != h_view.extent(4) ||
         d_view.extent(5) != h_view.extent(5) ||
         d_view.extent(6) != h_view.extent(6) ||
         d_view.extent(7) != h_view.extent(7) ||
         d_view.stride_0()    != h_view.stride_0() ||
         d_view.stride_1()    != h_view.stride_1() ||
         d_view.stride_2()    != h_view.stride_2() ||
         d_view.stride_3()    != h_view.stride_3() ||
         d_view.stride_4()    != h_view.stride_4() ||
         d_view.stride_5()    != h_view.stride_5() ||
         d_view.stride_6()    != h_view.stride_6() ||
         d_view.stride_7()    != h_view.stride_7() ||
         d_view.span()        != h_view.span() ) {
      Kokkos::Impl::throw_runtime_exception("DualView constructed with incompatible views");
    }
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
    modified_host = t_modified_flag(modified_flags,0);
    modified_device = t_modified_flag(modified_flags,1);
#endif
  }

  //@}
  //! \name Methods for synchronizing, marking as modified, and getting Views.
  //@{

  /// \brief Return a View on a specific device \c Device.
  ///
  /// Please don't be afraid of the if_c expression in the return
  /// value's type.  That just tells the method what the return type
  /// should be: t_dev if the \c Device template parameter matches
  /// this DualView's device type, else t_host.
  ///
  /// For example, suppose you create a DualView on Cuda, like this:
  /// \code
  /// typedef Kokkos::DualView<float, Kokkos::LayoutRight, Kokkos::Cuda> dual_view_type;
  /// dual_view_type DV ("my dual view", 100);
  /// \endcode
  /// If you want to get the CUDA device View, do this:
  /// \code
  /// typename dual_view_type::t_dev cudaView = DV.view<Kokkos::Cuda> ();
  /// \endcode
  /// and if you want to get the host mirror of that View, do this:
  /// \code
  /// typedef typename Kokkos::HostSpace::execution_space host_device_type;
  /// typename dual_view_type::t_host hostView = DV.view<host_device_type> ();
  /// \endcode
  template< class Device >
  KOKKOS_INLINE_FUNCTION
  const typename Impl::if_c<
    std::is_same<typename t_dev::memory_space,
                          typename Device::memory_space>::value,
    t_dev,
    t_host>::type& view () const
  {
    #ifndef KOKKOS_ENABLE_DEPRECATED_CODE
    constexpr bool device_is_memspace  = std::is_same<Device,typename Device::memory_space>::value;
    constexpr bool device_is_execspace = std::is_same<Device,typename Device::execution_space>::value;
    constexpr bool device_exec_is_t_dev_exec  = std::is_same<typename Device::execution_space,typename t_dev::execution_space>::value;
    constexpr bool device_mem_is_t_dev_mem    = std::is_same<typename Device::memory_space,typename t_dev::memory_space>::value;
    constexpr bool device_exec_is_t_host_exec  = std::is_same<typename Device::execution_space,typename t_host::execution_space>::value;
    constexpr bool device_mem_is_t_host_mem    = std::is_same<typename Device::memory_space,typename t_host::memory_space>::value;
    constexpr bool device_is_t_host_device  = std::is_same<typename Device::execution_space,typename t_host::device_type>::value;
    constexpr bool device_is_t_dev_device    = std::is_same<typename Device::memory_space,typename t_host::device_type>::value;

    static_assert(
        device_is_t_dev_device || device_is_t_host_device ||
        (device_is_memspace  && (device_mem_is_t_dev_mem   || device_mem_is_t_host_mem) ) ||
        (device_is_execspace && (device_exec_is_t_dev_exec || device_exec_is_t_host_exec) ) ||
        (
          (!device_is_execspace && !device_is_memspace) && (
            (device_mem_is_t_dev_mem   || device_mem_is_t_host_mem)  ||
            (device_exec_is_t_dev_exec || device_exec_is_t_host_exec)
          )
        )
        ,
        "Template parameter to .view() must exactly match one of the DualView's device types or one of the execution or memory spaces");
    #endif

    return Impl::if_c<
      std::is_same<
        typename t_dev::memory_space,
        typename Device::memory_space>::value,
      t_dev,
      t_host >::select (d_view , h_view);
  }

  KOKKOS_INLINE_FUNCTION
  t_host view_host() const {
    return h_view;
  }

  KOKKOS_INLINE_FUNCTION
  t_dev view_device() const {
    return d_view;
  }

  template<class Device>
  static int get_device_side() {
    constexpr bool device_is_memspace  = std::is_same<Device,typename Device::memory_space>::value;
    constexpr bool device_is_execspace = std::is_same<Device,typename Device::execution_space>::value;
    constexpr bool device_exec_is_t_dev_exec  = std::is_same<typename Device::execution_space,typename t_dev::execution_space>::value;
    constexpr bool device_mem_is_t_dev_mem    = std::is_same<typename Device::memory_space,typename t_dev::memory_space>::value;
    constexpr bool device_exec_is_t_host_exec  = std::is_same<typename Device::execution_space,typename t_host::execution_space>::value;
    constexpr bool device_mem_is_t_host_mem    = std::is_same<typename Device::memory_space,typename t_host::memory_space>::value;
    constexpr bool device_is_t_host_device  = std::is_same<typename Device::execution_space,typename t_host::device_type>::value;
    constexpr bool device_is_t_dev_device    = std::is_same<typename Device::memory_space,typename t_host::device_type>::value;

    #ifndef KOKKOS_ENABLE_DEPRECATED_CODE
    static_assert(
        device_is_t_dev_device || device_is_t_host_device ||
        (device_is_memspace  && (device_mem_is_t_dev_mem   || device_mem_is_t_host_mem) ) ||
        (device_is_execspace && (device_exec_is_t_dev_exec || device_exec_is_t_host_exec) ) ||
        (
          (!device_is_execspace && !device_is_memspace) && (
            (device_mem_is_t_dev_mem   || device_mem_is_t_host_mem)  ||
            (device_exec_is_t_dev_exec || device_exec_is_t_host_exec)
          )
        )
        ,
        "Template parameter to .sync() must exactly match one of the DualView's device types or one of the execution or memory spaces");
    #endif

    #ifndef KOKKOS_ENABLE_DEPRECATED_CODE
    int dev = -1;
    #else
    int dev = 0;
    #endif
    if(device_is_t_dev_device) dev = 1;
    else if(device_is_t_host_device) dev = 0;
    else {
      if(device_is_memspace) {
        if(device_mem_is_t_dev_mem) dev = 1;
        if(device_mem_is_t_host_mem) dev = 0;
        if(device_mem_is_t_host_mem && device_mem_is_t_dev_mem) dev = -1;
      }
      if(device_is_execspace) {
        if(device_exec_is_t_dev_exec) dev = 1;
        if(device_exec_is_t_host_exec) dev = 0;
        if(device_exec_is_t_host_exec && device_exec_is_t_dev_exec) dev = -1;
      }
      if(!device_is_execspace && !device_is_memspace) {
        if(device_mem_is_t_dev_mem) dev = 1;
        if(device_mem_is_t_host_mem) dev = 0;
        if(device_mem_is_t_host_mem && device_mem_is_t_dev_mem) dev = -1;
        if(device_exec_is_t_dev_exec) dev = 1;
        if(device_exec_is_t_host_exec) dev = 0;
        if(device_exec_is_t_host_exec && device_exec_is_t_dev_exec) dev = -1;
      }
    }
    return dev;
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
  template<class Device>
  void sync( const typename Impl::enable_if<
        ( std::is_same< typename traits::data_type , typename traits::non_const_data_type>::value) ||
        ( std::is_same< Device , int>::value)
        , int >::type& = 0)
  {
    if(modified_flags.data()==NULL) return;

    int dev = get_device_side<Device>();

    if (dev == 1) { // if Device is the same as DualView's device type
      if ((modified_flags(0) > 0) && (modified_flags(0) >= modified_flags(1))) {
        deep_copy (d_view, h_view);
        modified_flags(0) = modified_flags(1) = 0;
      }
    }
    if (dev == 0) { // hopefully Device is the same as DualView's host type
      if ((modified_flags(1) > 0) && (modified_flags(1) >= modified_flags(0))) {
        deep_copy (h_view, d_view);
        modified_flags(0) = modified_flags(1) = 0;
      }
    }
    if(std::is_same<typename t_host::memory_space,typename t_dev::memory_space>::value) {
      t_dev::execution_space::fence();
      t_host::execution_space::fence();
    }
  }

  template<class Device>
  void sync ( const typename Impl::enable_if<
      ( ! std::is_same< typename traits::data_type , typename traits::non_const_data_type>::value ) ||
      ( std::is_same< Device , int>::value)
      , int >::type& = 0 )
  {
    if(modified_flags.data()==NULL) return;

    int dev = get_device_side<Device>();

    if (dev == 1) { // if Device is the same as DualView's device type
      if ((modified_flags(0) > 0) && (modified_flags(0) >= modified_flags(1))) {
        Impl::throw_runtime_exception("Calling sync on a DualView with a const datatype.");
      }
    }
    if (dev == 0){ // hopefully Device is the same as DualView's host type
      if ((modified_flags(1) > 0) && (modified_flags(1) >= modified_flags(0))) {
        Impl::throw_runtime_exception("Calling sync on a DualView with a const datatype.");
      }
    }
  }

  void sync_host() {
    if( ! std::is_same< typename traits::data_type , typename traits::non_const_data_type>::value )
      Impl::throw_runtime_exception("Calling sync_host on a DualView with a const datatype.");
    if(modified_flags.data()==NULL) return;
    if(modified_flags(1) > modified_flags(0)) {
      deep_copy (h_view, d_view);
      modified_flags(1) = modified_flags(0) = 0;
    }
  }

  void sync_device() {
    if( ! std::is_same< typename traits::data_type , typename traits::non_const_data_type>::value )
      Impl::throw_runtime_exception("Calling sync_device on a DualView with a const datatype.");
    if(modified_flags.data()==NULL) return;
    if(modified_flags(0) > modified_flags(1)) {
      deep_copy (d_view, h_view);
      modified_flags(1) = modified_flags(0) = 0;
    }
  }

  template<class Device>
  bool need_sync() const
  {
    if(modified_flags.data()==NULL) return false;
    int dev = get_device_side<Device>();

    if (dev == 1) { // if Device is the same as DualView's device type
      if ((modified_flags(0) > 0) && (modified_flags(0) >= modified_flags(1))) {
        return true;
      }
    }
    if (dev == 0){ // hopefully Device is the same as DualView's host type
      if ((modified_flags(1) > 0) && (modified_flags(1) >= modified_flags(0))) {
        return true;
      }
    }
    return false;
  }

  inline bool need_sync_host() const {
    if(modified_flags.data()==NULL) return false;
    return modified_flags(0)<modified_flags(1);
  }

  inline bool need_sync_device() const {
    if(modified_flags.data()==NULL) return false;
    return modified_flags(1)<modified_flags(0);
  }

  /// \brief Mark data as modified on the given device \c Device.
  ///
  /// If \c Device is the same as this DualView's device type, then
  /// mark the device's data as modified.  Otherwise, mark the host's
  /// data as modified.
  template<class Device>
  void modify () {
    if(modified_flags.data()==NULL) return;
    int dev = get_device_side<Device>();

    if (dev == 1) { // if Device is the same as DualView's device type
      // Increment the device's modified count.
      modified_flags(1) = (modified_flags(1) > modified_flags(0) ?
                            modified_flags(1) : modified_flags(0)) + 1;
    }
    if (dev == 0) { // hopefully Device is the same as DualView's host type
      // Increment the host's modified count.
      modified_flags(0) = (modified_flags(1) > modified_flags(0) ?
                          modified_flags(1) : modified_flags(0))  + 1;
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

  inline void modify_host() {
    if(modified_flags.data()!=NULL) {
      modified_flags(0) = (modified_flags(1) > modified_flags(0) ?
          modified_flags(1) : modified_flags(0))  + 1;
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

  inline void modify_device() {
    if(modified_flags.data()!=NULL) {
      modified_flags(1) = (modified_flags(1) > modified_flags(0) ?
          modified_flags(1) : modified_flags(0))  + 1;
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

  inline void clear_sync_state() {
    if(modified_flags.data()!=NULL) 
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
  void realloc( const size_t n0 = KOKKOS_IMPL_CTOR_DEFAULT_ARG ,
           const size_t n1 = KOKKOS_IMPL_CTOR_DEFAULT_ARG ,
           const size_t n2 = KOKKOS_IMPL_CTOR_DEFAULT_ARG ,
           const size_t n3 = KOKKOS_IMPL_CTOR_DEFAULT_ARG ,
           const size_t n4 = KOKKOS_IMPL_CTOR_DEFAULT_ARG ,
           const size_t n5 = KOKKOS_IMPL_CTOR_DEFAULT_ARG ,
           const size_t n6 = KOKKOS_IMPL_CTOR_DEFAULT_ARG ,
           const size_t n7 = KOKKOS_IMPL_CTOR_DEFAULT_ARG ) {
    ::Kokkos::realloc(d_view,n0,n1,n2,n3,n4,n5,n6,n7);
     h_view = create_mirror_view( d_view );

     /* Reset dirty flags */
     if(modified_flags.data()==NULL) {
       modified_flags = t_modified_flags("DualView::modified_flags");
     } else
       modified_flags(1) = modified_flags(0) = 0;
  }

  /// \brief Resize both views, copying old contents into new if necessary.
  ///
  /// This method only copies the old contents into the new View
  /// objects for the device which was last marked as modified.
  void resize( const size_t n0 = KOKKOS_IMPL_CTOR_DEFAULT_ARG ,
           const size_t n1 = KOKKOS_IMPL_CTOR_DEFAULT_ARG ,
           const size_t n2 = KOKKOS_IMPL_CTOR_DEFAULT_ARG ,
           const size_t n3 = KOKKOS_IMPL_CTOR_DEFAULT_ARG ,
           const size_t n4 = KOKKOS_IMPL_CTOR_DEFAULT_ARG ,
           const size_t n5 = KOKKOS_IMPL_CTOR_DEFAULT_ARG ,
           const size_t n6 = KOKKOS_IMPL_CTOR_DEFAULT_ARG ,
           const size_t n7 = KOKKOS_IMPL_CTOR_DEFAULT_ARG ) {
   if(modified_flags.data()==NULL) {
     modified_flags = t_modified_flags("DualView::modified_flags");
   }
   if(modified_flags(1) >= modified_flags(0)) {
     /* Resize on Device */
     ::Kokkos::resize(d_view,n0,n1,n2,n3,n4,n5,n6,n7);
     h_view = create_mirror_view( d_view );

     /* Mark Device copy as modified */
     modified_flags(1) = modified_flags(1)+1;

   } else {
     /* Realloc on Device */

     ::Kokkos::realloc(d_view,n0,n1,n2,n3,n4,n5,n6,n7);

     const bool sizeMismatch = ( h_view.extent(0) != n0 ) ||
         ( h_view.extent(1) != n1 ) ||
         ( h_view.extent(2) != n2 ) ||
         ( h_view.extent(3) != n3 ) ||
         ( h_view.extent(4) != n4 ) ||
         ( h_view.extent(5) != n5 ) ||
         ( h_view.extent(6) != n6 ) ||
         ( h_view.extent(7) != n7 );
     if ( sizeMismatch )
       ::Kokkos::resize(h_view,n0,n1,n2,n3,n4,n5,n6,n7);

     t_host temp_view = create_mirror_view( d_view );

     /* Remap on Host */
     Kokkos::deep_copy( temp_view , h_view );

     h_view = temp_view;

     d_view = create_mirror_view( typename t_dev::execution_space(), h_view );

     /* Mark Host copy as modified */
     modified_flags(0) = modified_flags(0)+1;
   }
  }

  //@}
  //! \name Methods for getting capacity, stride, or dimension(s).
  //@{

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
  //! The allocation size (same as Kokkos::View::capacity).
  size_t capacity() const {
    return d_view.span();
  }
#endif

  //! The allocation size (same as Kokkos::View::span).
  KOKKOS_INLINE_FUNCTION constexpr size_t span() const {
    return d_view.span();
  }

  KOKKOS_INLINE_FUNCTION bool span_is_contiguous() const { 
    return d_view.span_is_contiguous(); 
  }

  //! Get stride(s) for each dimension.
  template< typename iType>
  void stride(iType* stride_) const {
    d_view.stride(stride_);
  }

  template< typename iType >
   KOKKOS_INLINE_FUNCTION constexpr
   typename std::enable_if< std::is_integral<iType>::value , size_t >::type
   extent( const iType & r ) const
     { return d_view.extent(r); }

   template< typename iType >
   KOKKOS_INLINE_FUNCTION constexpr
   typename std::enable_if< std::is_integral<iType>::value , int >::type
   extent_int( const iType & r ) const
     { return static_cast<int>(d_view.extent(r)); }

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
  /*  Deprecate all 'dimension' functions in favor of
   *  ISO/C++ vocabulary 'extent'.
   */

  /* \brief return size of dimension 0 */
  size_t dimension_0() const {return d_view.extent(0);}
  /* \brief return size of dimension 1 */
  size_t dimension_1() const {return d_view.extent(1);}
  /* \brief return size of dimension 2 */
  size_t dimension_2() const {return d_view.extent(2);}
  /* \brief return size of dimension 3 */
  size_t dimension_3() const {return d_view.extent(3);}
  /* \brief return size of dimension 4 */
  size_t dimension_4() const {return d_view.extent(4);}
  /* \brief return size of dimension 5 */
  size_t dimension_5() const {return d_view.extent(5);}
  /* \brief return size of dimension 6 */
  size_t dimension_6() const {return d_view.extent(6);}
  /* \brief return size of dimension 7 */
  size_t dimension_7() const {return d_view.extent(7);}
#endif

  //@}
};

} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//
// Partial specializations of Kokkos::subview() for DualView objects.
//

namespace Kokkos {
namespace Impl {

template< class D, class A1, class A2, class A3, class ... Args >
struct DualViewSubview {

  typedef typename Kokkos::Impl::ViewMapping
    < void
    , Kokkos::ViewTraits< D, A1, A2, A3 >
    , Args ...
    >::traits_type dst_traits ;

  typedef Kokkos::DualView
    < typename dst_traits::data_type
    , typename dst_traits::array_layout
    , typename dst_traits::device_type
    , typename dst_traits::memory_traits
    > type ;
};

} /* namespace Impl */


template< class D , class A1 , class A2 , class A3 , class ... Args >
typename Impl::DualViewSubview<D,A1,A2,A3,Args...>::type
subview( const DualView<D,A1,A2,A3> & src , Args ... args )
{
  return typename
    Impl::DualViewSubview<D,A1,A2,A3,Args...>::type( src , args ... );
}

} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

//
// Partial specialization of Kokkos::deep_copy() for DualView objects.
//

template< class DT , class DL , class DD , class DM ,
          class ST , class SL , class SD , class SM >
void
deep_copy (DualView<DT,DL,DD,DM> dst, // trust me, this must not be a reference
           const DualView<ST,SL,SD,SM>& src )
{
  if(src.modified_flags.data()==NULL || dst.modified_flags.data()==NULL) {
    return deep_copy(dst.d_view, src.d_view);
  }
  if (src.modified_flags(1) >= src.modified_flags(0)) {
    deep_copy (dst.d_view, src.d_view);
    dst.template modify<typename DualView<DT,DL,DD,DM>::device_type> ();
  } else {
    deep_copy (dst.h_view, src.h_view);
    dst.template modify<typename DualView<DT,DL,DD,DM>::host_mirror_space> ();
  }
}

template< class ExecutionSpace ,
          class DT , class DL , class DD , class DM ,
          class ST , class SL , class SD , class SM >
void
deep_copy (const ExecutionSpace& exec ,
           DualView<DT,DL,DD,DM> dst, // trust me, this must not be a reference
           const DualView<ST,SL,SD,SM>& src )
{
  if(src.modified_flags.data()==NULL || dst.modified_flags.data()==NULL) {
    return deep_copy(exec, dst.d_view, src.d_view);
  }
  if (src.modified_flags(1) >= src.modified_flags(0)) {
    deep_copy (exec, dst.d_view, src.d_view);
    dst.template modify<typename DualView<DT,DL,DD,DM>::device_type> ();
  } else {
    deep_copy (exec, dst.h_view, src.h_view);
    dst.template modify<typename DualView<DT,DL,DD,DM>::host_mirror_space> ();
  }
}

} // namespace Kokkos

#endif

