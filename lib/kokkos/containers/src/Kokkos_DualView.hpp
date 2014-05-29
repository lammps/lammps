/*
//@HEADER
// ************************************************************************
//
//                             Kokkos
//         Manycore Performance-Portable Multidimensional Arrays
//
//              Copyright (2012) Sandia Corporation
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
// Questions?  Contact  H. Carter Edwards (hcedwar@sandia.gov)
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

#include <Kokkos_View.hpp>

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
template< class T , class L , class D, class M = MemoryManaged>
class DualView {
public:
  //! \name Typedefs for device types and various Kokkos::View specializations.
  //@{

  //! The Kokkos Device type; same as the \c Device template parameter.
  typedef D device_type;
  //! The host mirror Kokkos Device type.
  typedef typename D::host_mirror_device_type host_mirror_device_type;

  //! The type of a Kokkos::View on the device.
  typedef Kokkos::View<T,L,D,M> t_dev ;

  /// \typedef t_host
  /// \brief The type of a Kokkos::View host mirror of \c t_dev.
#if defined( CUDA_VERSION ) && ( 6000 <= CUDA_VERSION ) && defined(KOKKOS_USE_UVM)
  typedef t_dev t_host;
#else
  typedef typename t_dev::HostMirror t_host ;
#endif

  //! The type of a const View on the device.
  typedef Kokkos::View<typename t_dev::const_data_type,L,D,M> t_dev_const;

  /// \typedef t_host_const
  /// \brief The type of a const View host mirror of \c t_dev_const.
#if defined( CUDA_VERSION ) && ( 6000 <= CUDA_VERSION ) && defined(KOKKOS_USE_UVM)
  typedef t_dev_const t_host_const;
#else
  typedef typename t_dev_const::HostMirror t_host_const;
#endif

  //! The type of a const, random-access View on the device.
  typedef Kokkos::View<typename t_dev::const_data_type,L,D,Kokkos::MemoryRandomAccess> t_dev_const_randomread ;

  /// \typedef t_host_const_randomread
  /// \brief The type of a const, random-access View host mirror of
  ///   \c t_dev_const_randomread.
#if defined( CUDA_VERSION ) && ( 6000 <= CUDA_VERSION ) && defined(KOKKOS_USE_UVM)
  typedef t_dev_const_randomread t_host_const_randomread;
#else
  typedef typename t_dev_const_randomread::HostMirror t_host_const_randomread;
#endif

  //! The type of an unmanaged View on the device.
  typedef Kokkos::View<T, L, D, Kokkos::MemoryUnmanaged> t_dev_um;
  //! The type of an unmanaged View host mirror of \c t_dev_um.
  typedef Kokkos::View<typename t_host::data_type,
                       typename t_host::array_layout,
                       typename t_host::device_type,
                       Kokkos::MemoryUnmanaged> t_host_um;

  //! The type of a const unmanaged View on the device.
  typedef Kokkos::View<typename t_dev::const_data_type, L, D,
                       Kokkos::MemoryUnmanaged> t_dev_const_um;
  //! The type of a const unmanaged View host mirror of \c t_dev_const_um.
  typedef Kokkos::View<typename t_host::const_data_type,
                       typename t_host::array_layout,
                       typename t_host::device_type,
                       Kokkos::MemoryUnmanaged> t_host_const_um;
  //@}
  //! \name The same typedefs as a View for scalar, data, and value types.
  //@{

  typedef typename t_dev::value_type value_type;
  typedef typename t_dev::const_value_type const_value_type;
  typedef typename t_dev::non_const_value_type non_const_value_type;

  //@}
  //! \name The two View instances.
  //@{

  t_dev d_view;
  t_host h_view;

  //@}
  //! \name Counters to keep track of changes ("modified" flags)
  //@{

  View<unsigned int,LayoutLeft,host_mirror_device_type> modified_device;
  View<unsigned int,LayoutLeft,host_mirror_device_type> modified_host;

  //@}
  //! \name Constructors
  //@{

  /// \brief Empty constructor.
  ///
  /// Both device and host View objects are constructed using their
  /// default constructors.  The "modified" flags are both initialized
  /// to "unmodified."
  DualView () :
    modified_device (View<unsigned int,LayoutLeft,host_mirror_device_type> ("DualView::modified_device")),
    modified_host (View<unsigned int,LayoutLeft,host_mirror_device_type> ("DualView::modified_host"))
  {}

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
            const size_t n0 = 0,
            const size_t n1 = 0,
            const size_t n2 = 0,
            const size_t n3 = 0,
            const size_t n4 = 0,
            const size_t n5 = 0,
            const size_t n6 = 0,
            const size_t n7 = 0)
    : d_view (label, n0, n1, n2, n3, n4, n5, n6, n7)
#if defined( CUDA_VERSION ) && ( 6000 <= CUDA_VERSION ) && defined(KOKKOS_USE_UVM)
    , h_view (d_view) // with UVM, host View is _always_ a shallow copy
#else
    , h_view (create_mirror_view (d_view)) // without UVM, host View mirrors
#endif
    , modified_device (View<unsigned int,LayoutLeft,host_mirror_device_type> ("DualView::modified_device"))
    , modified_host (View<unsigned int,LayoutLeft,host_mirror_device_type> ("DualView::modified_host"))
  {}

  //! Copy constructor (shallow copy)
  template<class SS, class LS, class DS, class MS>
  DualView (const DualView<SS,LS,DS,MS>& src) :
    d_view (src.d_view),
    h_view (src.h_view),
    modified_device (src.modified_device),
    modified_host (src.modified_host)
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
    modified_device (View<unsigned int,LayoutLeft,host_mirror_device_type> ("DualView::modified_device")),
    modified_host (View<unsigned int,LayoutLeft,host_mirror_device_type> ("DualView::modified_host"))
  {
    Impl::assert_shapes_are_equal (d_view.shape (), h_view.shape ());
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
  /// typedef typename Kokkos::Cuda::host_mirror_device_type host_device_type;
  /// typename dual_view_type::t_host hostView = DV.view<host_device_type> ();
  /// \endcode
  template< class Device >
  const typename Kokkos::Impl::if_c<
    Kokkos::Impl::is_same<typename t_dev::memory_space,
                          typename Device::memory_space>::value,
    t_dev,
    t_host>::type view () const
  {
    return Kokkos::Impl::if_c<
      Kokkos::Impl::is_same<
        typename t_dev::memory_space,
        typename Device::memory_space>::value,
      t_dev,
      t_host >::select (d_view , h_view);
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
  void sync () {
    const unsigned int dev =
      Kokkos::Impl::if_c<
        Kokkos::Impl::is_same<
          typename t_dev::memory_space,
          typename Device::memory_space>::value ,
        unsigned int,
        unsigned int>::select (1, 0);

    if (dev) { // if Device is the same as DualView's device type
      if ((modified_host () > 0) && (modified_host () >= modified_device ())) {
        Kokkos::deep_copy (d_view, h_view);
        modified_host() = modified_device() = 0;
      }
    } else { // hopefully Device is the same as DualView's host type
      if ((modified_device () > 0) && (modified_device () >= modified_host ())) {
        Kokkos::deep_copy (h_view, d_view);
        modified_host() = modified_device() = 0;
      }
    }
  }

  /// \brief Mark data as modified on the given device \c Device.
  ///
  /// If \c Device is the same as this DualView's device type, then
  /// mark the device's data as modified.  Otherwise, mark the host's
  /// data as modified.
  template<class Device>
  void modify () {
    const unsigned int dev =
      Kokkos::Impl::if_c<
        Kokkos::Impl::is_same<
          typename t_dev::memory_space,
          typename Device::memory_space>::value,
        unsigned int,
        unsigned int>::select (1, 0);

    if (dev) { // if Device is the same as DualView's device type
      // Increment the device's modified count.
      modified_device () = (modified_device () > modified_host () ?
                            modified_device () : modified_host ()) + 1;
    } else { // hopefully Device is the same as DualView's host type
      // Increment the host's modified count.
      modified_host () = (modified_device () > modified_host () ?
                          modified_device () : modified_host ())  + 1;
    }
  }

  //@}
  //! \name Methods for reallocating or resizing the View objects.
  //@{

  /// \brief Reallocate both View objects.
  ///
  /// This discards any existing contents of the objects, and resets
  /// their modified flags.  It does <i>not</i> copy the old contents
  /// of either View into the new View objects.
  void realloc( const size_t n0 = 0 ,
           const size_t n1 = 0 ,
           const size_t n2 = 0 ,
           const size_t n3 = 0 ,
           const size_t n4 = 0 ,
           const size_t n5 = 0 ,
           const size_t n6 = 0 ,
           const size_t n7 = 0 ) {
     Kokkos::realloc(d_view,n0,n1,n2,n3,n4,n5,n6,n7);
#if defined( CUDA_VERSION ) && ( 6000 <= CUDA_VERSION ) && defined(KOKKOS_USE_UVM)
     h_view = d_view ;
#else
     h_view = create_mirror_view( d_view );
#endif
     /* Reset dirty flags */
     modified_device() = modified_host() = 0;
  }

  /// \brief Resize both views, copying old contents into new if necessary.
  ///
  /// This method only copies the old contents into the new View
  /// objects for the device which was last marked as modified.
  void resize( const size_t n0 = 0 ,
           const size_t n1 = 0 ,
           const size_t n2 = 0 ,
           const size_t n3 = 0 ,
           const size_t n4 = 0 ,
           const size_t n5 = 0 ,
           const size_t n6 = 0 ,
           const size_t n7 = 0 ) {
   if(modified_device() >= modified_host()) {
     /* Resize on Device */
     Kokkos::resize(d_view,n0,n1,n2,n3,n4,n5,n6,n7);
#if defined( CUDA_VERSION ) && ( 6000 <= CUDA_VERSION ) && defined(KOKKOS_USE_UVM)
     h_view = d_view ;
#else
     h_view = create_mirror_view( d_view );
#endif

     /* Mark Device copy as modified */
     modified_device() = modified_device()+1;

   } else {
     /* Realloc on Device */

     Kokkos::realloc(d_view,n0,n1,n2,n3,n4,n5,n6,n7);
#if defined( CUDA_VERSION ) && ( 6000 <= CUDA_VERSION ) && defined(KOKKOS_USE_UVM)
     t_host temp_view = d_view ;
#else
     t_host temp_view = create_mirror_view( d_view );
#endif

     /* Remap on Host */
     Kokkos::Impl::ViewRemap< t_host , t_host >( temp_view , h_view );
     h_view = temp_view;

     /* Mark Host copy as modified */
     modified_host() = modified_host()+1;
   }
  }

  //@}
  //! \name Methods for getting capacity, stride, or dimension(s).
  //@{

  //! The allocation size (same as Kokkos::View::capacity).
  size_t capacity() const {
    return d_view.capacity();
  }

  //! Get stride(s) for each dimension.
  template< typename iType>
  void stride(iType* stride_) const {
    d_view.stride(stride_);
  }

  /* \brief return size of dimension 0 */
  size_t dimension_0() const {return d_view.dimension_0();}
  /* \brief return size of dimension 1 */
  size_t dimension_1() const {return d_view.dimension_1();}
  /* \brief return size of dimension 2 */
  size_t dimension_2() const {return d_view.dimension_2();}
  /* \brief return size of dimension 3 */
  size_t dimension_3() const {return d_view.dimension_3();}
  /* \brief return size of dimension 4 */
  size_t dimension_4() const {return d_view.dimension_4();}
  /* \brief return size of dimension 5 */
  size_t dimension_5() const {return d_view.dimension_5();}
  /* \brief return size of dimension 6 */
  size_t dimension_6() const {return d_view.dimension_6();}
  /* \brief return size of dimension 7 */
  size_t dimension_7() const {return d_view.dimension_7();}

  //@}
};

//
// Partial specializations of Kokkos::subview() for DualView objects.
//

template< class DstViewType ,
          class T , class L , class D , class M ,
          class ArgType0 >
DstViewType
subview( const DualView<T,L,D,M> & src ,
         const ArgType0 & arg0 )
{
  DstViewType sub_view;
  sub_view.d_view = subview<typename DstViewType::t_dev>(src.d_view,arg0);
  sub_view.h_view = subview<typename DstViewType::t_host>(src.h_view,arg0);
  sub_view.modified_device = src.modified_device;
  sub_view.modified_host = src.modified_host;
  return sub_view;
}


template< class DstViewType ,
          class T , class L , class D , class M ,
          class ArgType0 , class ArgType1 >
DstViewType
subview( const DualView<T,L,D,M> & src ,
         const ArgType0 & arg0 ,
         const ArgType1 & arg1 )
{
  DstViewType sub_view;
  sub_view.d_view = subview<typename DstViewType::t_dev>(src.d_view,arg0,arg1);
  sub_view.h_view = subview<typename DstViewType::t_host>(src.h_view,arg0,arg1);
  sub_view.modified_device = src.modified_device;
  sub_view.modified_host = src.modified_host;
  return sub_view;
}

template< class DstViewType ,
          class T , class L , class D , class M ,
          class ArgType0 , class ArgType1 , class ArgType2 >
DstViewType
subview( const DualView<T,L,D,M> & src ,
         const ArgType0 & arg0 ,
         const ArgType1 & arg1 ,
         const ArgType2 & arg2 )
{
  DstViewType sub_view;
  sub_view.d_view = subview<typename DstViewType::t_dev>(src.d_view,arg0,arg1,arg2);
  sub_view.h_view = subview<typename DstViewType::t_host>(src.h_view,arg0,arg1,arg2);
  sub_view.modified_device = src.modified_device;
  sub_view.modified_host = src.modified_host;
  return sub_view;
}

template< class DstViewType ,
          class T , class L , class D , class M ,
          class ArgType0 , class ArgType1 , class ArgType2 , class ArgType3 >
DstViewType
subview( const DualView<T,L,D,M> & src ,
         const ArgType0 & arg0 ,
         const ArgType1 & arg1 ,
         const ArgType2 & arg2 ,
         const ArgType3 & arg3 )
{
  DstViewType sub_view;
  sub_view.d_view = subview<typename DstViewType::t_dev>(src.d_view,arg0,arg1,arg2,arg3);
  sub_view.h_view = subview<typename DstViewType::t_host>(src.h_view,arg0,arg1,arg2,arg3);
  sub_view.modified_device = src.modified_device;
  sub_view.modified_host = src.modified_host;
  return sub_view;
}

template< class DstViewType ,
          class T , class L , class D , class M ,
          class ArgType0 , class ArgType1 , class ArgType2 , class ArgType3 ,
          class ArgType4 >
DstViewType
subview( const DualView<T,L,D,M> & src ,
         const ArgType0 & arg0 ,
         const ArgType1 & arg1 ,
         const ArgType2 & arg2 ,
         const ArgType3 & arg3 ,
         const ArgType4 & arg4 )
{
  DstViewType sub_view;
  sub_view.d_view = subview<typename DstViewType::t_dev>(src.d_view,arg0,arg1,arg2,arg3,arg4);
  sub_view.h_view = subview<typename DstViewType::t_host>(src.h_view,arg0,arg1,arg2,arg3,arg4);
  sub_view.modified_device = src.modified_device;
  sub_view.modified_host = src.modified_host;
  return sub_view;
}

template< class DstViewType ,
          class T , class L , class D , class M ,
          class ArgType0 , class ArgType1 , class ArgType2 , class ArgType3 ,
          class ArgType4 , class ArgType5 >
DstViewType
subview( const DualView<T,L,D,M> & src ,
         const ArgType0 & arg0 ,
         const ArgType1 & arg1 ,
         const ArgType2 & arg2 ,
         const ArgType3 & arg3 ,
         const ArgType4 & arg4 ,
         const ArgType5 & arg5 )
{
  DstViewType sub_view;
  sub_view.d_view = subview<typename DstViewType::t_dev>(src.d_view,arg0,arg1,arg2,arg3,arg4,arg5);
  sub_view.h_view = subview<typename DstViewType::t_host>(src.h_view,arg0,arg1,arg2,arg3,arg4,arg5);
  sub_view.modified_device = src.modified_device;
  sub_view.modified_host = src.modified_host;
  return sub_view;
}

template< class DstViewType ,
          class T , class L , class D , class M ,
          class ArgType0 , class ArgType1 , class ArgType2 , class ArgType3 ,
          class ArgType4 , class ArgType5 , class ArgType6 >
DstViewType
subview( const DualView<T,L,D,M> & src ,
         const ArgType0 & arg0 ,
         const ArgType1 & arg1 ,
         const ArgType2 & arg2 ,
         const ArgType3 & arg3 ,
         const ArgType4 & arg4 ,
         const ArgType5 & arg5 ,
         const ArgType6 & arg6 )
{
  DstViewType sub_view;
  sub_view.d_view = subview<typename DstViewType::t_dev>(src.d_view,arg0,arg1,arg2,arg3,arg4,arg5,arg6);
  sub_view.h_view = subview<typename DstViewType::t_host>(src.h_view,arg0,arg1,arg2,arg3,arg4,arg5,arg6);
  sub_view.modified_device = src.modified_device;
  sub_view.modified_host = src.modified_host;
  return sub_view;
}

template< class DstViewType ,
          class T , class L , class D , class M ,
          class ArgType0 , class ArgType1 , class ArgType2 , class ArgType3 ,
          class ArgType4 , class ArgType5 , class ArgType6 , class ArgType7 >
DstViewType
subview( const DualView<T,L,D,M> & src ,
         const ArgType0 & arg0 ,
         const ArgType1 & arg1 ,
         const ArgType2 & arg2 ,
         const ArgType3 & arg3 ,
         const ArgType4 & arg4 ,
         const ArgType5 & arg5 ,
         const ArgType6 & arg6 ,
         const ArgType7 & arg7 )
{
  DstViewType sub_view;
  sub_view.d_view = subview<typename DstViewType::t_dev>(src.d_view,arg0,arg1,arg2,arg3,arg4,arg5,arg6,arg7);
  sub_view.h_view = subview<typename DstViewType::t_host>(src.h_view,arg0,arg1,arg2,arg3,arg4,arg5,arg6,arg7);
  sub_view.modified_device = src.modified_device;
  sub_view.modified_host = src.modified_host;
  return sub_view;
}

//
// Partial specialization of Kokkos::deep_copy() for DualView objects.
//

template< class DT , class DL , class DD , class DM ,
          class ST , class SL , class SD , class SM >
void
deep_copy (DualView<DT,DL,DD,DM> dst, // trust me, this must not be a reference
           const DualView<ST,SL,SD,SM>& src)
{
  if (src.modified_device () >= src.modified_host ()) {
    Kokkos::deep_copy (dst.d_view, src.d_view);
    dst.template modify<typename DualView<DT,DL,DD,DM>::device_type> ();
  } else {
    Kokkos::deep_copy (dst.h_view, src.h_view);
    dst.template modify<typename DualView<DT,DL,DD,DM>::host_mirror_device_type> ();
  }
}

} // namespace Kokkos

#endif
