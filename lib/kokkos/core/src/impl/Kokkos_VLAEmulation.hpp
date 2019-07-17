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

#ifndef KOKKOS_IMPL_VLAEMULATION_HPP
#define KOKKOS_IMPL_VLAEMULATION_HPP

#include <Kokkos_Macros.hpp>
#if defined( KOKKOS_ENABLE_TASKDAG )


#include <Kokkos_Core_fwd.hpp>

#include <impl/Kokkos_Error.hpp> // KOKKOS_EXPECTS

#include <type_traits> // std::is_abstract<>, ...

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template <
  class Derived,
  class VLAValueType,
  class EntryCountType = int32_t
>
struct ObjectWithVLAEmulation;

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

/** @brief Attorney to enable private CRTP inheritance from ObjectWithVLAEmulation
 */
struct VLAEmulationAccess {
private:

  template <class, class, class>
  friend struct ObjectWithVLAEmulation;

  template <class Derived, class VLAValueType, class EntryCountType>
  KOKKOS_FORCEINLINE_FUNCTION
  static constexpr Derived*
  _cast_to_derived(ObjectWithVLAEmulation<Derived, VLAValueType, EntryCountType>* base) noexcept
  {
    return static_cast<Derived*>(base);
  }

  template <class Derived, class VLAValueType, class EntryCountType>
  KOKKOS_FORCEINLINE_FUNCTION
  static constexpr Derived const*
  _cast_to_derived(ObjectWithVLAEmulation<Derived, VLAValueType, EntryCountType> const* base) noexcept
  {
    return static_cast<Derived const*>(base);
  }

};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

/** \brief A CRTP base class for a type that includes a variable-length array by allocation
 *
 *  The storage for the derived type must be allocated manually and the objects
 *  (both derived type and VLA objects) must be constructed with placement new.
 *  Obviously, this can't be done for objects on the stack.
 *
 *  Note: Though most uses of this currently delete the copy and move constructor
 *  in the `Derived` type, this type is intended to have value semantics.
 *
 *  \todo @documentation elaborate on implications of value semantics for this class template
 *
 */
template <
  class Derived,
  class VLAValueType,
  class EntryCountType /* = int32_t */
>
struct ObjectWithVLAEmulation {
public:

  using object_type = Derived;
  using vla_value_type = VLAValueType;
  using vla_entry_count_type = EntryCountType;

  using iterator = VLAValueType*;
  using const_iterator = typename std::add_const<VLAValueType>::type*;


  // TODO @tasking @minor DSH require that Derived be marked final? (note that std::is_final is C++14)
  // TODO @tasking @minor DSH delete non-placement operator new for Derived type?

private:

  vla_entry_count_type m_num_entries;

  // CRTP boilerplate

  KOKKOS_FORCEINLINE_FUNCTION
  /* KOKKOS_CONSTEXPR_14 */
  Derived* _this() noexcept { return VLAEmulationAccess::_cast_to_derived(this); }

  KOKKOS_FORCEINLINE_FUNCTION
  /* KOKKOS_CONSTEXPR_14 */
  Derived const* _this() const noexcept { return VLAEmulationAccess::_cast_to_derived(this); }

  // Note: can't be constexpr because of reinterpret_cast
  KOKKOS_FORCEINLINE_FUNCTION
  /* KOKKOS_CONSTEXPR_14 */
  vla_value_type* _vla_pointer() noexcept {
    // The data starts right after the aligned storage of Derived
    return reinterpret_cast<vla_value_type*>(_this() + 1);
  }

  // Note: can't be constexpr because of reinterpret_cast
  KOKKOS_FORCEINLINE_FUNCTION
  /* KOKKOS_CONSTEXPR_14 */
  vla_value_type const* _vla_pointer() const noexcept {
    // The data starts right after the aligned storage of Derived
    return reinterpret_cast<vla_value_type const*>(_this() + 1);
  }

public:

  KOKKOS_INLINE_FUNCTION
  static /* KOKKOS_CONSTEXPR_14 */ size_t
  required_allocation_size(vla_entry_count_type num_vla_entries) {
    KOKKOS_EXPECTS(num_vla_entries >= 0);
    return sizeof(Derived) + num_vla_entries * sizeof(VLAValueType);
  }

  //----------------------------------------------------------------------------
  // <editor-fold desc="Constructors, destructor, and assignment"> {{{2

  // TODO @tasking @optimization DSH specialization for trivially constructible VLAValueType?
  // TODO @tasking @minor DSH SFINAE-out this constructor for non-default contructible vla_value_types
  KOKKOS_INLINE_FUNCTION
  explicit
  ObjectWithVLAEmulation(vla_entry_count_type num_entries)
    noexcept(noexcept(vla_value_type()))
    : m_num_entries(num_entries)
  {
    // Note: We can't do this at class scope because it unnecessarily requires
    // object_type to be a complete type
    static_assert(
      alignof(object_type) >= alignof(vla_value_type),
      "Can't append emulated variable length array of type with greater alignment than"
      "  the type to which the VLA is being appended"
    );

    // Note: We can't do this at class scope because it unnecessarily requires
    // vla_value_type to be a complete type
    static_assert(
      not std::is_abstract<vla_value_type>::value,
      "Can't use abstract type with VLA emulation"
    );

    KOKKOS_EXPECTS(num_entries >= 0);
    for(vla_entry_count_type i = 0; i < m_num_entries; ++i) {
      new (_vla_pointer() + i) vla_value_type();
    }
  }

  KOKKOS_INLINE_FUNCTION
  ~ObjectWithVLAEmulation()
    noexcept(noexcept(std::declval<vla_value_type>().~vla_value_type()))
  {
    for(auto&& value : *this) { value.~vla_value_type(); }
  }

  // TODO @tasking @new_feature DSH constrained analogs for move and copy ctors and assignment ops
  // TODO @tasking @new_feature DSH forwarding in_place constructor
  // TODO @tasking @new_feature DSH initializer_list constructor?

  // </editor-fold> end Constructors, destructor, and assignment }}}2
  //----------------------------------------------------------------------------


  KOKKOS_INLINE_FUNCTION
  constexpr EntryCountType n_vla_entries() const noexcept { return m_num_entries; }


  //----------------------------------------------------------------------------
  // <editor-fold desc="Accessing the object and the VLA values"> {{{2

  KOKKOS_INLINE_FUNCTION
  object_type& object() & { return static_cast<Derived&>(*this); }

  KOKKOS_INLINE_FUNCTION
  object_type const& object() const & { return static_cast<Derived const&>(*this); }

  KOKKOS_INLINE_FUNCTION
  object_type&& object() && { return static_cast<Derived&&>(*this); }


  KOKKOS_INLINE_FUNCTION
  vla_value_type& vla_value_at(vla_entry_count_type n) &
  {
    KOKKOS_EXPECTS(n < n_vla_entries());
    return _vla_pointer()[n];
  }

  KOKKOS_INLINE_FUNCTION
  vla_value_type const& vla_value_at(vla_entry_count_type n) const &
  {
    KOKKOS_EXPECTS(n < n_vla_entries());
    return _vla_pointer()[n];
  }

  KOKKOS_INLINE_FUNCTION
  vla_value_type& vla_value_at(vla_entry_count_type n) &&
  {
    KOKKOS_EXPECTS(n < n_vla_entries());
    return _vla_pointer()[n];
  }

  // </editor-fold> end Accessing the object and the VLA values }}}2
  //----------------------------------------------------------------------------


  //----------------------------------------------------------------------------
  // <editor-fold desc="Iterators"> {{{2

  KOKKOS_INLINE_FUNCTION
  iterator begin() noexcept { return _vla_pointer(); }

  KOKKOS_INLINE_FUNCTION
  const_iterator begin() const noexcept { return _vla_pointer(); }

  KOKKOS_INLINE_FUNCTION
  const_iterator cbegin() noexcept { return _vla_pointer(); }

  KOKKOS_INLINE_FUNCTION
  iterator end() noexcept { return _vla_pointer() + m_num_entries; }

  KOKKOS_INLINE_FUNCTION
  const_iterator end() const noexcept { return _vla_pointer() + m_num_entries; }

  KOKKOS_INLINE_FUNCTION
  const_iterator cend() noexcept { return _vla_pointer() + m_num_entries; }

  // </editor-fold> end Iterators }}}2
  //----------------------------------------------------------------------------

};

} /* namespace Impl */
} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #if defined( KOKKOS_ENABLE_TASKDAG ) */
#endif /* #ifndef KOKKOS_IMPL_VLAEMULATION_HPP */

