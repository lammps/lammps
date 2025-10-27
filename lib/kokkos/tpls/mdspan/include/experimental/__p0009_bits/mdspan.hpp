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

#pragma once

#include "default_accessor.hpp"
#include "layout_right.hpp"
#include "extents.hpp"
#include "trait_backports.hpp"
#include "compressed_pair.hpp"

namespace MDSPAN_IMPL_STANDARD_NAMESPACE {
template <
  class ElementType,
  class Extents,
  class LayoutPolicy = layout_right,
  class AccessorPolicy = default_accessor<ElementType>
>
class mdspan
{
private:
  static_assert(detail::impl_is_extents_v<Extents>,
                MDSPAN_IMPL_STANDARD_NAMESPACE_STRING "::mdspan's Extents template parameter must be a specialization of " MDSPAN_IMPL_STANDARD_NAMESPACE_STRING "::extents.");
  static_assert(std::is_same<ElementType, typename AccessorPolicy::element_type>::value,
                MDSPAN_IMPL_STANDARD_NAMESPACE_STRING "::mdspan's ElementType template parameter must be the same as its AccessorPolicy::element_type.");

  // Workaround for non-deducibility of the index sequence template parameter if it's given at the top level
  template <class>
  struct deduction_workaround;

  template <size_t... Idxs>
  struct deduction_workaround<std::index_sequence<Idxs...>>
  {
    MDSPAN_FORCE_INLINE_FUNCTION static constexpr
    size_t size(mdspan const& self) noexcept {
      return MDSPAN_IMPL_FOLD_TIMES_RIGHT((self.mapping_ref().extents().extent(Idxs)), /* * ... * */ size_t(1));
    }
    MDSPAN_FORCE_INLINE_FUNCTION static constexpr
    bool empty(mdspan const& self) noexcept {
      return (self.rank()>0) && MDSPAN_IMPL_FOLD_OR((self.mapping_ref().extents().extent(Idxs)==index_type(0)));
    }
    template <class ReferenceType, class SizeType, size_t N>
    MDSPAN_FORCE_INLINE_FUNCTION static constexpr
    ReferenceType callop(mdspan const& self, const std::array<SizeType, N>& indices) noexcept {
      return self.accessor_ref().access(self.ptr_ref(), self.mapping_ref()(indices[Idxs]...));
    }
#ifdef __cpp_lib_span
    template <class ReferenceType, class SizeType, size_t N>
    MDSPAN_FORCE_INLINE_FUNCTION static constexpr
    ReferenceType callop(mdspan const& self, const std::span<SizeType, N>& indices) noexcept {
      return self.accessor_ref().access(self.ptr_ref(), self.mapping_ref()(indices[Idxs]...));
    }
#endif
  };

public:

  //--------------------------------------------------------------------------------
  // Domain and codomain types

  using extents_type = Extents;
  using layout_type = LayoutPolicy;
  using accessor_type = AccessorPolicy;
  using mapping_type = typename layout_type::template mapping<extents_type>;
  using element_type = ElementType;
  using value_type = std::remove_cv_t<element_type>;
  using index_type = typename extents_type::index_type;
  using size_type = typename extents_type::size_type;
  using rank_type = typename extents_type::rank_type;
  using data_handle_type = typename accessor_type::data_handle_type;
  using reference = typename accessor_type::reference;

  MDSPAN_INLINE_FUNCTION static constexpr size_t rank() noexcept { return extents_type::rank(); }
  MDSPAN_INLINE_FUNCTION static constexpr size_t rank_dynamic() noexcept { return extents_type::rank_dynamic(); }
  MDSPAN_INLINE_FUNCTION static constexpr size_t static_extent(size_t r) noexcept { return extents_type::static_extent(r); }
  MDSPAN_INLINE_FUNCTION constexpr index_type extent(size_t r) const noexcept { return mapping_ref().extents().extent(r); };

private:

  // Can't use defaulted parameter in the deduction_workaround template because of a bug in MSVC warning C4348.
  using deduction_workaround_impl = deduction_workaround<std::make_index_sequence<extents_type::rank()>>;

  using map_acc_pair_t = detail::impl_compressed_pair<mapping_type, accessor_type>;

public:

  //--------------------------------------------------------------------------------
  // [mdspan.basic.cons], mdspan constructors, assignment, and destructor

#if !MDSPAN_HAS_CXX_20
  MDSPAN_INLINE_FUNCTION_DEFAULTED constexpr mdspan() = default;
#else
  MDSPAN_INLINE_FUNCTION_DEFAULTED constexpr mdspan()
    requires(
       // nvhpc has a bug where using just rank_dynamic() here doesn't work ...
       (extents_type::rank_dynamic() > 0) &&
       MDSPAN_IMPL_TRAIT(std::is_default_constructible, data_handle_type) &&
       MDSPAN_IMPL_TRAIT(std::is_default_constructible, mapping_type) &&
       MDSPAN_IMPL_TRAIT(std::is_default_constructible, accessor_type)
     ) = default;
#endif
  MDSPAN_INLINE_FUNCTION_DEFAULTED constexpr mdspan(const mdspan&) = default;
  MDSPAN_INLINE_FUNCTION_DEFAULTED constexpr mdspan(mdspan&&) = default;

  MDSPAN_TEMPLATE_REQUIRES(
    class... SizeTypes,
    /* requires */ (
      ((sizeof...(SizeTypes) == rank()) || (sizeof...(SizeTypes) == rank_dynamic())) &&
      (detail::are_valid_indices<index_type, SizeTypes...>()) &&
      MDSPAN_IMPL_TRAIT(std::is_constructible, mapping_type, extents_type) &&
      MDSPAN_IMPL_TRAIT(std::is_default_constructible, accessor_type)
    )
  )
  MDSPAN_INLINE_FUNCTION
  explicit constexpr mdspan(data_handle_type p, SizeTypes... dynamic_extents)
    // TODO @proposal-bug shouldn't I be allowed to do `move(p)` here?
    : m_members(std::move(p), map_acc_pair_t(mapping_type(extents_type(static_cast<index_type>(std::move(dynamic_extents))...)), accessor_type()))
  { }

  MDSPAN_TEMPLATE_REQUIRES(
    class SizeType, size_t N,
    /* requires */ (
      MDSPAN_IMPL_TRAIT(std::is_convertible, const SizeType&, index_type) &&
      MDSPAN_IMPL_TRAIT(std::is_nothrow_constructible, index_type, const SizeType&) &&
      ((N == rank()) || (N == rank_dynamic())) &&
      MDSPAN_IMPL_TRAIT(std::is_constructible, mapping_type, extents_type) &&
      MDSPAN_IMPL_TRAIT(std::is_default_constructible, accessor_type)
    )
  )
  MDSPAN_CONDITIONAL_EXPLICIT(N != rank_dynamic())
  MDSPAN_INLINE_FUNCTION
  constexpr mdspan(data_handle_type p, const std::array<SizeType, N>& dynamic_extents)
    : m_members(std::move(p), map_acc_pair_t(mapping_type(extents_type(dynamic_extents)), accessor_type()))
  { }

#ifdef __cpp_lib_span
  MDSPAN_TEMPLATE_REQUIRES(
    class SizeType, size_t N,
    /* requires */ (
      MDSPAN_IMPL_TRAIT(std::is_convertible, const SizeType&, index_type) &&
      MDSPAN_IMPL_TRAIT(std::is_nothrow_constructible, index_type, const SizeType&) &&
      ((N == rank()) || (N == rank_dynamic())) &&
      MDSPAN_IMPL_TRAIT(std::is_constructible, mapping_type, extents_type) &&
      MDSPAN_IMPL_TRAIT(std::is_default_constructible, accessor_type)
    )
  )
  MDSPAN_CONDITIONAL_EXPLICIT(N != rank_dynamic())
  MDSPAN_INLINE_FUNCTION
  constexpr mdspan(data_handle_type p, std::span<SizeType, N> dynamic_extents)
    : m_members(std::move(p), map_acc_pair_t(mapping_type(extents_type(as_const(dynamic_extents))), accessor_type()))
  { }
#endif

  MDSPAN_FUNCTION_REQUIRES(
    (MDSPAN_INLINE_FUNCTION constexpr),
    mdspan, (data_handle_type p, const extents_type& exts), ,
    /* requires */ (MDSPAN_IMPL_TRAIT(std::is_default_constructible, accessor_type) &&
                    MDSPAN_IMPL_TRAIT(std::is_constructible, mapping_type, const extents_type&))
  ) : m_members(std::move(p), map_acc_pair_t(mapping_type(exts), accessor_type()))
  { }

  MDSPAN_FUNCTION_REQUIRES(
    (MDSPAN_INLINE_FUNCTION constexpr),
    mdspan, (data_handle_type p, const mapping_type& m), ,
    /* requires */ (MDSPAN_IMPL_TRAIT(std::is_default_constructible, accessor_type))
  ) : m_members(std::move(p), map_acc_pair_t(m, accessor_type()))
  { }

  MDSPAN_INLINE_FUNCTION
  constexpr mdspan(data_handle_type p, const mapping_type& m, const accessor_type& a)
    : m_members(std::move(p), map_acc_pair_t(m, a))
  { }

  MDSPAN_TEMPLATE_REQUIRES(
    class OtherElementType, class OtherExtents, class OtherLayoutPolicy, class OtherAccessor,
    /* requires */ (
      MDSPAN_IMPL_TRAIT(std::is_constructible, mapping_type, const typename OtherLayoutPolicy::template mapping<OtherExtents>&) &&
      MDSPAN_IMPL_TRAIT(std::is_constructible, accessor_type, const OtherAccessor&)
    )
  )
  MDSPAN_CONDITIONAL_EXPLICIT(
    !MDSPAN_IMPL_TRAIT(std::is_convertible, const typename OtherLayoutPolicy::template mapping<OtherExtents>&, mapping_type) ||
    !MDSPAN_IMPL_TRAIT(std::is_convertible, const OtherAccessor&, accessor_type)
  )
  MDSPAN_INLINE_FUNCTION
  constexpr mdspan(const mdspan<OtherElementType, OtherExtents, OtherLayoutPolicy, OtherAccessor>& other)
    : m_members(other.ptr_ref(), map_acc_pair_t(other.mapping_ref(), other.accessor_ref()))
  {
      static_assert(MDSPAN_IMPL_TRAIT(std::is_constructible, data_handle_type, typename OtherAccessor::data_handle_type),"Incompatible data_handle_type for mdspan construction");
      static_assert(MDSPAN_IMPL_TRAIT(std::is_constructible, extents_type, OtherExtents),"Incompatible extents for mdspan construction");
      /*
       * TODO: Check precondition
       * For each rank index r of extents_type, static_extent(r) == dynamic_extent || static_extent(r) == other.extent(r) is true.
       */
  }

  /* Might need this on NVIDIA?
  MDSPAN_INLINE_FUNCTION_DEFAULTED
  ~mdspan() = default;
  */

  MDSPAN_INLINE_FUNCTION_DEFAULTED MDSPAN_IMPL_CONSTEXPR_14_DEFAULTED mdspan& operator=(const mdspan&) = default;
  MDSPAN_INLINE_FUNCTION_DEFAULTED MDSPAN_IMPL_CONSTEXPR_14_DEFAULTED mdspan& operator=(mdspan&&) = default;


  //--------------------------------------------------------------------------------
  // [mdspan.basic.mapping], mdspan mapping domain multidimensional index to access codomain element

  #if MDSPAN_USE_BRACKET_OPERATOR
  MDSPAN_TEMPLATE_REQUIRES(
    class... SizeTypes,
    /* requires */ (
      MDSPAN_IMPL_FOLD_AND(MDSPAN_IMPL_TRAIT(std::is_convertible, SizeTypes, index_type) /* && ... */) &&
      MDSPAN_IMPL_FOLD_AND(MDSPAN_IMPL_TRAIT(std::is_nothrow_constructible, index_type, SizeTypes) /* && ... */) &&
      (rank() == sizeof...(SizeTypes))
    )
  )
  MDSPAN_FORCE_INLINE_FUNCTION
  constexpr reference operator[](SizeTypes... indices) const
  {
    return accessor_ref().access(ptr_ref(), mapping_ref()(static_cast<index_type>(std::move(indices))...));
  }
  #endif

  MDSPAN_TEMPLATE_REQUIRES(
    class SizeType,
    /* requires */ (
      MDSPAN_IMPL_TRAIT(std::is_convertible, const SizeType&, index_type) &&
      MDSPAN_IMPL_TRAIT(std::is_nothrow_constructible, index_type, const SizeType&)
    )
  )
  MDSPAN_FORCE_INLINE_FUNCTION
  constexpr reference operator[](const std::array< SizeType, rank()>& indices) const
  {
    return deduction_workaround_impl::template callop<reference>(*this, indices);
  }

  #ifdef __cpp_lib_span
  MDSPAN_TEMPLATE_REQUIRES(
    class SizeType,
    /* requires */ (
      MDSPAN_IMPL_TRAIT(std::is_convertible, const SizeType&, index_type) &&
      MDSPAN_IMPL_TRAIT(std::is_nothrow_constructible, index_type, const SizeType&)
    )
  )
  MDSPAN_FORCE_INLINE_FUNCTION
  constexpr reference operator[](std::span<SizeType, rank()> indices) const
  {
    return deduction_workaround_impl::template callop<reference>(*this, indices);
  }
  #endif // __cpp_lib_span

  #if !MDSPAN_USE_BRACKET_OPERATOR
  MDSPAN_TEMPLATE_REQUIRES(
    class Index,
    /* requires */ (
      MDSPAN_IMPL_TRAIT(std::is_convertible, Index, index_type) &&
      MDSPAN_IMPL_TRAIT(std::is_nothrow_constructible, index_type, Index) &&
      extents_type::rank() == 1
    )
  )
  MDSPAN_FORCE_INLINE_FUNCTION
  constexpr reference operator[](Index idx) const
  {
    return accessor_ref().access(ptr_ref(), mapping_ref()(static_cast<index_type>(std::move(idx))));
  }
  #endif

  #if MDSPAN_USE_PAREN_OPERATOR
  MDSPAN_TEMPLATE_REQUIRES(
    class... SizeTypes,
    /* requires */ (
      extents_type::rank() == sizeof...(SizeTypes) &&
      (detail::are_valid_indices<index_type, SizeTypes...>())
    )
  )
  MDSPAN_FORCE_INLINE_FUNCTION
  constexpr reference operator()(SizeTypes... indices) const
  {
    return accessor_ref().access(ptr_ref(), mapping_ref()(static_cast<index_type>(std::move(indices))...));
  }

  MDSPAN_TEMPLATE_REQUIRES(
    class SizeType,
    /* requires */ (
      MDSPAN_IMPL_TRAIT(std::is_convertible, const SizeType&, index_type) &&
      MDSPAN_IMPL_TRAIT(std::is_nothrow_constructible, index_type, const SizeType&)
    )
  )
  MDSPAN_FORCE_INLINE_FUNCTION
  constexpr reference operator()(const std::array<SizeType, rank()>& indices) const
  {
    return deduction_workaround_impl::template callop<reference>(*this, indices);
  }

  #ifdef __cpp_lib_span
  MDSPAN_TEMPLATE_REQUIRES(
    class SizeType,
    /* requires */ (
      MDSPAN_IMPL_TRAIT(std::is_convertible, const SizeType&, index_type) &&
      MDSPAN_IMPL_TRAIT(std::is_nothrow_constructible, index_type, const SizeType&)
    )
  )
  MDSPAN_FORCE_INLINE_FUNCTION
  constexpr reference operator()(std::span<SizeType, rank()> indices) const
  {
    return deduction_workaround_impl::template callop<reference>(*this, indices);
  }
  #endif // __cpp_lib_span
  #endif // MDSPAN_USE_PAREN_OPERATOR

  MDSPAN_INLINE_FUNCTION constexpr size_type size() const noexcept {
    return static_cast<size_type>(deduction_workaround_impl::size(*this));
  };

  MDSPAN_INLINE_FUNCTION constexpr bool empty() const noexcept {
    return deduction_workaround_impl::empty(*this);
  };

  MDSPAN_INLINE_FUNCTION
  friend constexpr void swap(mdspan& x, mdspan& y) noexcept {
    // can't call the std::swap inside on HIP
    #if !defined(MDSPAN_IMPL_HAS_HIP) && !defined(MDSPAN_IMPL_HAS_CUDA)
    using std::swap;
    swap(x.ptr_ref(), y.ptr_ref());
    swap(x.mapping_ref(), y.mapping_ref());
    swap(x.accessor_ref(), y.accessor_ref());
    #else
    mdspan tmp = y;
    y = x;
    x = tmp;
    #endif
  }

  //--------------------------------------------------------------------------------
  // [mdspan.basic.domobs], mdspan observers of the domain multidimensional index space


  MDSPAN_INLINE_FUNCTION constexpr const extents_type& extents() const noexcept { return mapping_ref().extents(); };
  MDSPAN_INLINE_FUNCTION constexpr const data_handle_type& data_handle() const noexcept { return ptr_ref(); };
  MDSPAN_INLINE_FUNCTION constexpr const mapping_type& mapping() const noexcept { return mapping_ref(); };
  MDSPAN_INLINE_FUNCTION constexpr const accessor_type& accessor() const noexcept { return accessor_ref(); };

  //--------------------------------------------------------------------------------
  // [mdspan.basic.obs], mdspan observers of the mapping

  MDSPAN_INLINE_FUNCTION static constexpr bool is_always_unique() { return mapping_type::is_always_unique(); };
  MDSPAN_INLINE_FUNCTION static constexpr bool is_always_exhaustive() { return mapping_type::is_always_exhaustive(); };
  MDSPAN_INLINE_FUNCTION static constexpr bool is_always_strided() { return mapping_type::is_always_strided(); };

  MDSPAN_INLINE_FUNCTION constexpr bool is_unique() const { return mapping_ref().is_unique(); };
  MDSPAN_INLINE_FUNCTION constexpr bool is_exhaustive() const { return mapping_ref().is_exhaustive(); };
  MDSPAN_INLINE_FUNCTION constexpr bool is_strided() const { return mapping_ref().is_strided(); };
  MDSPAN_INLINE_FUNCTION constexpr index_type stride(size_t r) const { return mapping_ref().stride(r); };

private:

  detail::impl_compressed_pair<data_handle_type, map_acc_pair_t> m_members{};

  MDSPAN_FORCE_INLINE_FUNCTION MDSPAN_IMPL_CONSTEXPR_14 data_handle_type& ptr_ref() noexcept { return m_members.first(); }
  MDSPAN_FORCE_INLINE_FUNCTION constexpr data_handle_type const& ptr_ref() const noexcept { return m_members.first(); }
  MDSPAN_FORCE_INLINE_FUNCTION MDSPAN_IMPL_CONSTEXPR_14 mapping_type& mapping_ref() noexcept { return m_members.second().first(); }
  MDSPAN_FORCE_INLINE_FUNCTION constexpr mapping_type const& mapping_ref() const noexcept { return m_members.second().first(); }
  MDSPAN_FORCE_INLINE_FUNCTION MDSPAN_IMPL_CONSTEXPR_14 accessor_type& accessor_ref() noexcept { return m_members.second().second(); }
  MDSPAN_FORCE_INLINE_FUNCTION constexpr accessor_type const& accessor_ref() const noexcept { return m_members.second().second(); }

  template <class, class, class, class>
  friend class mdspan;

};

#if defined(MDSPAN_IMPL_USE_CLASS_TEMPLATE_ARGUMENT_DEDUCTION)
MDSPAN_TEMPLATE_REQUIRES(
  class ElementType, class... SizeTypes,
  /* requires */ MDSPAN_IMPL_FOLD_AND(MDSPAN_IMPL_TRAIT(std::is_convertible, SizeTypes, size_t) /* && ... */) &&
  (sizeof...(SizeTypes) > 0)
)
MDSPAN_DEDUCTION_GUIDE explicit mdspan(ElementType*, SizeTypes...)
  -> mdspan<ElementType, ::MDSPAN_IMPL_STANDARD_NAMESPACE::dextents<size_t, sizeof...(SizeTypes)>>;

MDSPAN_TEMPLATE_REQUIRES(
  class Pointer,
  (MDSPAN_IMPL_TRAIT(std::is_pointer, std::remove_reference_t<Pointer>))
)
MDSPAN_DEDUCTION_GUIDE mdspan(Pointer&&) -> mdspan<std::remove_pointer_t<std::remove_reference_t<Pointer>>, extents<size_t>>;

MDSPAN_TEMPLATE_REQUIRES(
  class CArray,
  (MDSPAN_IMPL_TRAIT(std::is_array, CArray) && (std::rank_v<CArray> == 1))
)
MDSPAN_DEDUCTION_GUIDE mdspan(CArray&) -> mdspan<std::remove_all_extents_t<CArray>, extents<size_t, ::std::extent_v<CArray,0>>>;

template <class ElementType, class SizeType, size_t N>
MDSPAN_DEDUCTION_GUIDE mdspan(ElementType*, const ::std::array<SizeType, N>&)
  -> mdspan<ElementType, ::MDSPAN_IMPL_STANDARD_NAMESPACE::dextents<size_t, N>>;

#ifdef __cpp_lib_span
template <class ElementType, class SizeType, size_t N>
MDSPAN_DEDUCTION_GUIDE mdspan(ElementType*, ::std::span<SizeType, N>)
  -> mdspan<ElementType, ::MDSPAN_IMPL_STANDARD_NAMESPACE::dextents<size_t, N>>;
#endif

// This one is necessary because all the constructors take `data_handle_type`s, not
// `ElementType*`s, and `data_handle_type` is taken from `accessor_type::data_handle_type`, which
// seems to throw off automatic deduction guides.
template <class ElementType, class SizeType, size_t... ExtentsPack>
MDSPAN_DEDUCTION_GUIDE mdspan(ElementType*, const extents<SizeType, ExtentsPack...>&)
  -> mdspan<ElementType, ::MDSPAN_IMPL_STANDARD_NAMESPACE::extents<SizeType, ExtentsPack...>>;

template <class ElementType, class MappingType>
MDSPAN_DEDUCTION_GUIDE mdspan(ElementType*, const MappingType&)
  -> mdspan<ElementType, typename MappingType::extents_type, typename MappingType::layout_type>;

template <class MappingType, class AccessorType>
MDSPAN_DEDUCTION_GUIDE mdspan(const typename AccessorType::data_handle_type, const MappingType&, const AccessorType&)
  -> mdspan<typename AccessorType::element_type, typename MappingType::extents_type, typename MappingType::layout_type, AccessorType>;
#endif

} // end namespace MDSPAN_IMPL_STANDARD_NAMESPACE
