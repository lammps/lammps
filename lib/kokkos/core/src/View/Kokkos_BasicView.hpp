// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
static_assert(false,
              "Including non-public Kokkos header files is not allowed.");
#endif

#ifndef KOKKOS_BASIC_VIEW_HPP
#define KOKKOS_BASIC_VIEW_HPP
#include <Kokkos_Macros.hpp>
#include <impl/Kokkos_InitializeFinalize.hpp>
#include <impl/Kokkos_Utilities.hpp>
#include <impl/Kokkos_SharedAlloc.hpp>
#include <View/Kokkos_ViewAlloc.hpp>
#include <View/Kokkos_ViewAccessPreconditionsCheck.hpp>
#include <View/Kokkos_ViewCtor.hpp>
#include <View/Kokkos_ViewTraits.hpp>
#include <View/MDSpan/Kokkos_MDSpan_Header.hpp>
#include <View/MDSpan/Kokkos_MDSpan_Accessor.hpp>
#include <View/MDSpan/Kokkos_MDSpan_Layout.hpp>

#include <optional>
#include <type_traits>

#if defined(KOKKOS_ENABLE_DEBUG_BOUNDS_CHECK)

#define KOKKOS_IMPL_BASICVIEW_OPERATOR_VERIFY(...)                             \
  if constexpr (Impl::IsReferenceCountedDataHandle<data_handle_type>::value) { \
    Kokkos::Impl::runtime_check_memory_access_violation<memory_space>(         \
        m_ptr.tracker());                                                      \
    Kokkos::Impl::view_verify_operator_bounds(                                 \
        m_ptr.tracker(), m_map.extents(), m_ptr.get(), __VA_ARGS__);           \
  } else {                                                                     \
    Kokkos::Impl::runtime_check_memory_access_violation<memory_space>(         \
        Kokkos::Impl::SharedAllocationTracker());                              \
    Kokkos::Impl::view_verify_operator_bounds(                                 \
        Kokkos::Impl::SharedAllocationTracker(), m_map.extents(), m_ptr,       \
        __VA_ARGS__);                                                          \
  }

#else

#define KOKKOS_IMPL_BASICVIEW_OPERATOR_VERIFY(...)                             \
  if constexpr (Impl::IsReferenceCountedDataHandle<data_handle_type>::value) { \
    Kokkos::Impl::runtime_check_memory_access_violation<memory_space>(         \
        m_ptr.tracker());                                                      \
  } else {                                                                     \
    Kokkos::Impl::runtime_check_memory_access_violation<memory_space>(         \
        Kokkos::Impl::SharedAllocationTracker());                              \
  }
#endif

#define KOKKOS_IMPL_NO_UNIQUE_ADDRESS MDSPAN_IMPL_NO_UNIQUE_ADDRESS
namespace Kokkos::Impl {

constexpr inline struct SubViewCtorTag {
  explicit SubViewCtorTag() = default;
} subview_ctor_tag{};

template <class T>
struct KokkosSliceToMDSpanSliceImpl {
  using type = T;
  KOKKOS_FUNCTION
  static constexpr decltype(auto) transform(const T &s) { return s; }
};

template <>
struct KokkosSliceToMDSpanSliceImpl<Kokkos::ALL_t> {
  using type = full_extent_t;
  KOKKOS_FUNCTION
  static constexpr decltype(auto) transform(Kokkos::ALL_t) {
    return full_extent;
  }
};

template <class T>
using kokkos_slice_to_mdspan_slice =
    typename KokkosSliceToMDSpanSliceImpl<T>::type;

template <class T>
KOKKOS_INLINE_FUNCTION constexpr decltype(auto)
transform_kokkos_slice_to_mdspan_slice(const T &s) {
  return KokkosSliceToMDSpanSliceImpl<T>::transform(s);
}

// Default implementation for computing allocation size (in #of elements)
// from mapping and accessor.
// This is used to figure out how much data to allocate, with this particular
// overload suitable for situations where the allocation element type is
// the element type of the View (something which is not true for example for
// Sacado).
template <class MappingType, class AccessorType>
KOKKOS_INLINE_FUNCTION constexpr size_t
allocation_size_from_mapping_and_accessor(const MappingType &map,
                                          const AccessorType &) {
  return map.required_span_size();
}

// Tag type to enable ADL for accessor_from_mapping_and_accessor_arg
// customization point
template <class AccessorType>
struct AccessorTypeTag {};

// Default implementation for creating an accessor from a mapping
// and an AccessorArg_t.
// In Sacado the accessor construction may require information from
// the mapping (specifically the span size) in some cases.
// Specifically it needs it if the elements of a FAD type are not
// consecutive but strided by the span size.
template <class AccessorType, class MappingType>
KOKKOS_INLINE_FUNCTION constexpr auto accessor_from_mapping_and_accessor_arg(
    const Kokkos::Impl::AccessorTypeTag<AccessorType> &, const MappingType &,
    const AccessorArg_t &arg) {
  return AccessorType(arg.value);
}

// FIXME spurious warnings like
// error: 'SR.14123' may be used uninitialized [-Werror=maybe-uninitialized]
#if defined(KOKKOS_COMPILER_GNU) && KOKKOS_COMPILER_GNU >= 1500
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#pragma GCC diagnostic ignored "-Wuninitialized"
#endif

// BasicView has to be in a different namespace than Impl;
// The reason for this is that if BasicView is in Impl, View (which derives from
// BasicView) can cause function resolution to ADL-find the Kokkos::Impl
// namespace, i.e., an unqualified call can find an internal Kokkos function.
// This was already exhibited in some Kokkos functions that were named the same
// in both Kokkos:: and Kokkos::Impl:: namespaces and caused an ambiguous call
namespace BV {
template <class ElementType, class Extents, class LayoutPolicy,
          class AccessorPolicy>
class BasicView {
 public:
  using mdspan_type =
      mdspan<ElementType, Extents, LayoutPolicy, AccessorPolicy>;
  using extents_type  = typename mdspan_type::extents_type;
  using layout_type   = typename mdspan_type::layout_type;
  using accessor_type = typename mdspan_type::accessor_type;
  using mapping_type  = typename mdspan_type::mapping_type;
  using element_type  = typename mdspan_type::element_type;
  using value_type    = typename mdspan_type::value_type;
  // FIXME: backwards compatibility, should be changed to the same as mdspan
  // index_type
  using index_type       = typename mdspan_type::size_type;
  using size_type        = typename mdspan_type::size_type;
  using rank_type        = typename mdspan_type::rank_type;
  using data_handle_type = typename mdspan_type::data_handle_type;
  using reference        = typename mdspan_type::reference;
  using memory_space     = typename accessor_type::memory_space;
  using execution_space  = typename memory_space::execution_space;

  KOKKOS_FUNCTION static constexpr rank_type rank() noexcept {
    return extents_type::rank();
  }
  KOKKOS_FUNCTION static constexpr rank_type rank_dynamic() noexcept {
    return extents_type::rank_dynamic();
  }
  KOKKOS_FUNCTION static constexpr size_t static_extent(rank_type r) noexcept {
    return extents_type::static_extent(r);
  }
  KOKKOS_FUNCTION constexpr index_type extent(rank_type r) const noexcept {
    return m_map.extents().extent(r);
  }

 protected:
  // These are pre-condition checks which are unconditionally (i.e. in release
  // mode) enabled in Kokkos::View 4.4
  template <class OtherMapping>
  KOKKOS_FUNCTION static constexpr void check_basic_view_constructibility(
      [[maybe_unused]] const OtherMapping &rhs) {
    using src_t                           = typename OtherMapping::layout_type;
    using dst_t                           = layout_type;
    [[maybe_unused]] constexpr size_t rnk = mdspan_type::rank();
    if constexpr (!std::is_same_v<src_t, dst_t>) {
      if constexpr (Impl::IsLayoutLeftPadded<dst_t>::value) {
        if constexpr (std::is_same_v<src_t, layout_stride>) {
          index_type stride = 1;
          for (size_t r = 0; r < rnk; r++) {
            if (rhs.stride(r) != stride)
              Kokkos::abort("View assignment must have compatible layouts");
            if constexpr (rnk > 1)
              stride *= (r == 0 ? rhs.stride(1) : rhs.extents().extent(r));
          }
        }
      }
      if constexpr (Impl::IsLayoutRightPadded<dst_t>::value) {
        if constexpr (std::is_same_v<src_t, layout_stride>) {
          index_type stride = 1;
          if constexpr (rnk > 0) {
            for (size_t r = rnk; r > 0; r--) {
              if (rhs.stride(r - 1) != stride)
                Kokkos::abort("View assignment must have compatible layouts");
              if constexpr (rnk > 1)
                stride *= (r == rnk ? rhs.stride(r - 2)
                                    : rhs.extents().extent(r - 1));
            }
          }
        }
      }
      if constexpr (std::is_same_v<dst_t, layout_left>) {
        if constexpr (std::is_same_v<src_t, layout_stride>) {
          index_type stride = 1;
          for (size_t r = 0; r < rnk; r++) {
            if (rhs.stride(r) != stride)
              Kokkos::abort("View assignment must have compatible layouts");
            stride *= rhs.extents().extent(r);
          }
        } else if constexpr (Impl::IsLayoutLeftPadded<src_t>::value &&
                             rnk > 1) {
          if (rhs.stride(1) != rhs.extents().extent(0))
            Kokkos::abort("View assignment must have compatible layouts");
        }
      }
      if constexpr (std::is_same_v<dst_t, layout_right>) {
        if constexpr (std::is_same_v<src_t, layout_stride>) {
          index_type stride = 1;
          if constexpr (rnk > 0) {
            for (size_t r = rnk; r > 0; r--) {
              if (rhs.stride(r - 1) != stride)
                Kokkos::abort("View assignment must have compatible layouts");
              stride *= rhs.extents().extent(r - 1);
            }
          }
        } else if constexpr (Impl::IsLayoutRightPadded<src_t>::value &&
                             rnk > 1) {
          if (rhs.stride(rnk - 2) != rhs.extents().extent(rnk - 1))
            Kokkos::abort("View assignment must have compatible layouts");
        }
      }
    }
  }

 public:
  KOKKOS_DEFAULTED_FUNCTION constexpr BasicView() = default;

  KOKKOS_FUNCTION constexpr BasicView(const mdspan_type &other)
      : m_ptr(other.data_handle()),
        m_map(other.mapping()),
        m_acc(other.accessor()) {}
  KOKKOS_FUNCTION constexpr BasicView(mdspan_type &&other)
      : m_ptr(std::move(other.data_handle())),
        m_map(std::move(other.mapping())),
        m_acc(std::move(other.accessor())) {}

  template <class OtherIndexType, size_t Size>
  // When doing C++20 we should switch to this, the conditional explicit we
  // can't do in 17
  //    requires(std::is_constructible_v<mdspan_type, data_handle_type,
  //                                     std::array<OtherIndexType, Size>>)
  //  explicit(Size != rank_dynamic())
  KOKKOS_FUNCTION constexpr BasicView(
      std::enable_if_t<
          std::is_constructible_v<mdspan_type, data_handle_type,
                                  std::array<OtherIndexType, Size>>,
          data_handle_type>
          p,
      const Array<OtherIndexType, Size> &exts)
      : m_ptr(std::move(p)), m_map(extents_type(exts)), m_acc{} {}

  KOKKOS_FUNCTION constexpr BasicView(data_handle_type p,
                                      const extents_type &exts)
    requires(std::is_default_constructible_v<accessor_type> &&
             std::is_constructible_v<mapping_type, const extents_type &>)
      : m_ptr(std::move(p)), m_map(exts), m_acc{} {}

  KOKKOS_FUNCTION constexpr BasicView(data_handle_type p, const mapping_type &m)
    requires(std::is_default_constructible_v<accessor_type>)
      : m_ptr(std::move(p)), m_map(m), m_acc{} {}

  KOKKOS_FUNCTION constexpr BasicView(data_handle_type p, const mapping_type &m,
                                      const accessor_type &a)
      : m_ptr(std::move(p)), m_map(m), m_acc(a) {}

  template <class OtherT, class OtherE, class OtherL, class OtherA>
    requires(std::is_constructible_v<mdspan_type,
                                     typename BasicView<OtherT, OtherE, OtherL,
                                                        OtherA>::mdspan_type>)
  explicit(
      !std::is_convertible_v<const typename OtherL::template mapping<OtherE> &,
                             mapping_type> ||
      !std::is_convertible_v<const OtherA &, accessor_type>)
      KOKKOS_INLINE_FUNCTION
      BasicView(const BasicView<OtherT, OtherE, OtherL, OtherA> &other)
      : m_ptr(other.m_ptr), m_map(other.m_map), m_acc(other.m_acc) {
    // Kokkos View precondition checks happen in release builds
    check_basic_view_constructibility(other.mapping());

    static_assert(
        std::is_constructible_v<data_handle_type,
                                const typename OtherA::data_handle_type &>,
        "Kokkos::View: incompatible data_handle_type for View construction");
    static_assert(std::is_constructible_v<extents_type, OtherE>,
                  "Kokkos::View: incompatible extents for View construction");
  }

  template <class OtherT, class OtherE, class OtherL, class OtherA>
    requires(std::is_constructible_v<mdspan_type,
                                     mdspan<OtherT, OtherE, OtherL, OtherA>>)
  explicit(
      !std::is_convertible_v<const typename OtherL::template mapping<OtherE> &,
                             mapping_type> ||
      !std::is_convertible_v<const OtherA &, accessor_type>)
      KOKKOS_INLINE_FUNCTION
      BasicView(const mdspan<OtherT, OtherE, OtherL, OtherA> &other,
                std::enable_if_t<
                    std::is_constructible_v<
                        mdspan_type, mdspan<OtherT, OtherE, OtherL, OtherA>>,
                    void *> = nullptr)
      : m_ptr(other.data_handle()),
        m_map(other.mapping()),
        m_acc(other.accessor()) {
    // Kokkos View precondition checks happen in release builds
    check_basic_view_constructibility(other.mapping());

    static_assert(
        std::is_constructible_v<data_handle_type,
                                const typename OtherA::data_handle_type &>,
        "Kokkos::View: incompatible data_handle_type for View construction");
    static_assert(std::is_constructible_v<extents_type, OtherE>,
                  "Kokkos::View: incompatible extents for View construction");
  }

  // Allocating constructors specific to BasicView
  ///
  /// Construct from a given mapping
  ///
  explicit constexpr BasicView(const std::string &label,
                               const mapping_type &mapping)
      : BasicView(view_alloc(label), mapping) {}

  ///
  /// Construct from a given extents
  ///
  explicit constexpr BasicView(const std::string &label,
                               const extents_type &ext)
      : BasicView(view_alloc(label), mapping_type{ext}) {}

 private:
  template <class... P>
  data_handle_type create_data_handle(
      const Impl::ViewCtorProp<P...> &arg_prop,
      const typename mdspan_type::mapping_type &arg_mapping,
      const typename mdspan_type::accessor_type &arg_accessor) {
    if (bool was_finalized = is_finalized();
        was_finalized || !is_initialized()) {
      std::stringstream ss;
      ss << "Kokkos ERROR: View ";
      constexpr bool has_label = Impl::ViewCtorProp<P...>::has_label;
      if constexpr (has_label) {
        auto const &lbl = Impl::get_property<Impl::LabelTag>(arg_prop);
        ss << "(label=\"" << lbl << "\") ";
      }
      ss << "is being constructed ";
      if (was_finalized) {
        ss << "after finalize() ";
      } else {
        ss << "before initialize() ";
      }
      ss << "has been called";
      auto const err = ss.str();
      abort(err.c_str());
    }
    using storage_value_type = typename data_handle_type::value_type;
    constexpr bool has_exec  = Impl::ViewCtorProp<P...>::has_execution_space;
    // Copy the input allocation properties with possibly defaulted properties
    // We need to split it in two to avoid MSVC compiler errors
    auto prop_copy_tmp =
        Impl::with_properties_if_unset(arg_prop, std::string{});
    auto prop_copy = Impl::with_properties_if_unset(
        prop_copy_tmp, memory_space{}, execution_space{});
    using alloc_prop = decltype(prop_copy);

    // get allocation size: may be different from
    // arg_mapping.required_span_size()
    size_t allocation_size =
        allocation_size_from_mapping_and_accessor(arg_mapping, arg_accessor);
    if constexpr (has_exec) {
      return data_handle_type{
          Impl::make_shared_allocation_record<storage_value_type>(
              allocation_size, Impl::get_property<Impl::LabelTag>(prop_copy),
              Impl::get_property<Impl::MemorySpaceTag>(prop_copy),
              std::make_optional(
                  Impl::get_property<Impl::ExecutionSpaceTag>(prop_copy)),
              std::bool_constant<alloc_prop::initialize>(),
              std::bool_constant<alloc_prop::sequential_host_init>())};
    } else {
      return data_handle_type{
          Impl::make_shared_allocation_record<storage_value_type>(
              allocation_size, Impl::get_property<Impl::LabelTag>(prop_copy),
              Impl::get_property<Impl::MemorySpaceTag>(prop_copy),
              std::optional<execution_space>{},
              std::bool_constant<alloc_prop::initialize>(),
              std::bool_constant<alloc_prop::sequential_host_init>())};
    }
  }

 public:
  // Ctors to pull out AccessorArg_t
  // Need also the other ones to keep the overload set consistent and all the
  // constraints mutually exclusive. Delegate to private ctors
  // We need to explicitly distinguish between the has_pointer and !has_pointer
  // versions since only the ones with a pointer can be marked host/device
  template <class... P>
  explicit BasicView(
      const Impl::ViewCtorProp<P...> &arg_prop,
      std::enable_if_t<!Impl::ViewCtorProp<P...>::has_pointer &&
                           Impl::ViewCtorProp<P...>::has_accessor_arg,
                       typename mdspan_type::mapping_type> const &arg_mapping)
      : BasicView(arg_prop, arg_mapping,
                  accessor_from_mapping_and_accessor_arg(
                      Impl::AccessorTypeTag<accessor_type>(), arg_mapping,
                      Impl::get_property<Impl::AccessorArgTag>(arg_prop))) {}

  template <class... P>
  KOKKOS_FUNCTION explicit BasicView(
      const Impl::ViewCtorProp<P...> &arg_prop,
      std::enable_if_t<Impl::ViewCtorProp<P...>::has_pointer &&
                           Impl::ViewCtorProp<P...>::has_accessor_arg,
                       typename mdspan_type::mapping_type> const &arg_mapping)
      : BasicView(arg_prop, arg_mapping,
                  accessor_from_mapping_and_accessor_arg(
                      Impl::AccessorTypeTag<accessor_type>(), arg_mapping,
                      Impl::get_property<Impl::AccessorArgTag>(arg_prop))) {}

  // Private Ctors coming from the case where we dealt with AccessorArg_t
  // We don't have public ctors which take ViewCtorProp and accessor.
  // I don't want to create the data handle above, because in a subsequent PR
  // the data handle creation itself will depend on the constructed accessor.
 private:
  template <class... P>
  explicit BasicView(
      const Impl::ViewCtorProp<P...> &arg_prop,
      std::enable_if_t<!Impl::ViewCtorProp<P...>::has_pointer &&
                           Impl::ViewCtorProp<P...>::has_accessor_arg,
                       typename mdspan_type::mapping_type> const &arg_mapping,
      const accessor_type &arg_accessor)
      : BasicView(create_data_handle(arg_prop, arg_mapping, arg_accessor),
                  arg_mapping, arg_accessor) {}

  template <class... P>
  KOKKOS_FUNCTION explicit BasicView(
      const Impl::ViewCtorProp<P...> &arg_prop,
      std::enable_if_t<ViewCtorProp<P...>::has_pointer &&
                           Impl::ViewCtorProp<P...>::has_accessor_arg,
                       typename mdspan_type::mapping_type> const &arg_mapping,
      const accessor_type &arg_accessor)
      : BasicView(
            data_handle_type{Impl::get_property<Impl::PointerTag>(arg_prop)},
            arg_mapping, arg_accessor) {}

  // Ctors taking CtorProp that don't have AccessorArg_t in it
 public:
  template <class... P>
  explicit BasicView(
      const Impl::ViewCtorProp<P...> &arg_prop,
      std::enable_if_t<!Impl::ViewCtorProp<P...>::has_pointer &&
                           !Impl::ViewCtorProp<P...>::has_accessor_arg,
                       typename mdspan_type::mapping_type> const &arg_mapping,
      const accessor_type &arg_accessor)
      : BasicView(create_data_handle(arg_prop, arg_mapping, arg_accessor),
                  arg_mapping, arg_accessor) {}

  template <class... P>
  KOKKOS_FUNCTION explicit BasicView(
      const Impl::ViewCtorProp<P...> &arg_prop,
      std::enable_if_t<ViewCtorProp<P...>::has_pointer &&
                           !Impl::ViewCtorProp<P...>::has_accessor_arg,
                       typename mdspan_type::mapping_type> const &arg_mapping,
      const accessor_type &arg_accessor)
      : BasicView(
            data_handle_type{Impl::get_property<Impl::PointerTag>(arg_prop)},
            arg_mapping, arg_accessor) {}

  template <class... P>
  explicit inline BasicView(
      const Impl::ViewCtorProp<P...> &arg_prop,
      std::enable_if_t<!Impl::ViewCtorProp<P...>::has_pointer &&
                           !Impl::ViewCtorProp<P...>::has_accessor_arg,
                       typename mdspan_type::mapping_type> const &arg_mapping)
      : BasicView(create_data_handle(arg_prop, arg_mapping, accessor_type{}),
                  arg_mapping) {}

  template <class... P>
  KOKKOS_FUNCTION explicit inline BasicView(
      const Impl::ViewCtorProp<P...> &arg_prop,
      std::enable_if_t<Impl::ViewCtorProp<P...>::has_pointer &&
                           !Impl::ViewCtorProp<P...>::has_accessor_arg,
                       typename mdspan_type::mapping_type> const &arg_mapping)
      : BasicView(
            data_handle_type{Impl::get_property<Impl::PointerTag>(arg_prop)},
            arg_mapping) {}

#ifdef KOKKOS_ENABLE_DEBUG_BOUNDS_CHECK
 private:
  template <class Arg>
  KOKKOS_INLINE_FUNCTION static bool is_in_bounds(size_t extent,
                                                  const Arg arg) {
    return static_cast<std::size_t>(arg) < extent;
  }

  template <class T>
  KOKKOS_INLINE_FUNCTION static bool is_in_bounds(size_t extent,
                                                  const std::pair<T, T> arg) {
    if constexpr (std::is_signed_v<T>) {
      return static_cast<std::size_t>(arg.second) <= extent && arg.first >= 0 &&
             arg.first <= arg.second;
    } else {
      return static_cast<std::size_t>(arg.second) <= extent &&
             arg.first <= arg.second;
    }
  }

  template <class T>
  KOKKOS_INLINE_FUNCTION static bool is_in_bounds(
      size_t extent, const Kokkos::pair<T, T> arg) {
    if constexpr (std::is_signed_v<T>) {
      return static_cast<std::size_t>(arg.second) <= extent && arg.first >= 0 &&
             arg.first <= arg.second;
    } else {
      return static_cast<std::size_t>(arg.second) <= extent &&
             arg.first <= arg.second;
    }
  }

  KOKKOS_INLINE_FUNCTION static bool is_in_bounds(size_t /*extent*/,
                                                  const Kokkos::ALL_t) {
    return true;
  }

  template <class RT, class... RP, size_t... Idx, class... Args>
  KOKKOS_INLINE_FUNCTION static bool subview_extents_valid(
      const BasicView<RT, RP...> &src_view, std::index_sequence<Idx...>,
      Args... args) {
    return (is_in_bounds(src_view.extent(Idx), args) && ... && true);
  }

  template <class Arg>
  static void append_error_message(std::stringstream &ss, size_t extent,
                                   const Arg arg) {
    ss << arg << " < " << extent;
  }

  template <class T>
  static void append_error_message(std::stringstream &ss, size_t extent,
                                   const std::pair<T, T> arg) {
    ss << arg.first << " <= " << arg.second << " <= " << extent;
  }

  template <class T>
  static void append_error_message(std::stringstream &ss, size_t extent,
                                   const Kokkos::pair<T, T> arg) {
    ss << arg.first << " <= " << arg.second << " <= " << extent;
  }

  static void append_error_message(std::stringstream &ss, size_t /*extent*/,
                                   const Kokkos::ALL_t) {
    ss << "Kokkos::ALL";
  }

  template <class RT, class... RP, size_t... Idx, class Arg0, class... Args>
  static std::stringstream generate_error_message(
      const BasicView<RT, RP...> &src_view, std::index_sequence<Idx...>,
      Arg0 arg0, Args... args) {
    std::stringstream ss;
    ss << "Kokkos::subview bounds error (";
    append_error_message(ss, src_view.extent(0), arg0);
    ((ss << ", ", append_error_message(ss, src_view.extent(Idx + 1), args)),
     ...);
    ss << ')';
    return ss;
  }
#endif

 protected:
  // Optimized subview constructor that calls submdspan_mapping directly
  // and initializes members in the constructor body
  template <class OtherElementType, class OtherExtents, class OtherLayoutPolicy,
            class OtherAccessorPolicy, class... SliceSpecifiers>
  KOKKOS_INLINE_FUNCTION BasicView(
      Impl::SubViewCtorTag,
      const BasicView<OtherElementType, OtherExtents, OtherLayoutPolicy,
                      OtherAccessorPolicy> &src_view,
      SliceSpecifiers... slices) {
    // Avoid calling submdspan to not create temporary mdspan objects.
    // Instead do what submdspan does: calling submdspan_mapping.
    const auto sub_mapping_result = submdspan_mapping(
        src_view.m_map,
        Impl::transform_kokkos_slice_to_mdspan_slice(slices)...);

    // Kokkos View precondition should happen in release build,
    // and before any checks happen inside mdspan mapping ctor itself.
    check_basic_view_constructibility(sub_mapping_result.mapping);

    // Explicit cast is needed because submdspan_mapping may return a different
    // layout type.
    using sub_accessor_t = typename OtherAccessorPolicy::offset_policy;
    m_ptr = src_view.m_acc.offset(src_view.m_ptr, sub_mapping_result.offset);
    m_map = mapping_type(sub_mapping_result.mapping);
    m_acc = sub_accessor_t(src_view.m_acc);

#ifdef KOKKOS_ENABLE_DEBUG_BOUNDS_CHECK
    bool valid = subview_extents_valid(
        src_view, std::make_index_sequence<sizeof...(SliceSpecifiers)>{},
        slices...);
    if (!valid) {
      KOKKOS_IF_ON_HOST(
          Kokkos::abort(
              generate_error_message(
                  src_view,
                  std::make_index_sequence<sizeof...(SliceSpecifiers) - 1>{},
                  slices...)
                  .str()
                  .c_str());)
      KOKKOS_IF_ON_DEVICE(Kokkos::abort("Kokkos::subview bounds error");)
    }
#endif
  }

 public:
  //----------------------------------------
  // Conversion to MDSpan
  template <class OtherElementType, class OtherExtents, class OtherLayoutPolicy,
            class OtherAccessor,
            typename = std::enable_if_t<
                std::is_assignable_v<mdspan<OtherElementType, OtherExtents,
                                            OtherLayoutPolicy, OtherAccessor>,
                                     mdspan_type>>>
  KOKKOS_INLINE_FUNCTION constexpr
  operator mdspan<OtherElementType, OtherExtents, OtherLayoutPolicy,
                  OtherAccessor>() const {
    return mdspan_type(m_ptr, m_map, m_acc);
  }

  // Here we use an overload instead of a default parameter as a workaround
  // to a potential compiler bug with clang 17. It may be present in other
  // compilers
  template <class OtherAccessorType = AccessorPolicy,
            typename                = std::enable_if_t<std::is_constructible_v<
                typename mdspan_type::data_handle_type,
                typename OtherAccessorType::data_handle_type>>>
  KOKKOS_INLINE_FUNCTION constexpr auto to_mdspan() const {
    using ret_mdspan_type =
        mdspan<typename mdspan_type::element_type,
               typename mdspan_type::extents_type,
               typename mdspan_type::layout_type, OtherAccessorType>;
    return ret_mdspan_type(
        static_cast<typename OtherAccessorType::data_handle_type>(
            data_handle()),
        mapping(), static_cast<OtherAccessorType>(accessor()));
  }

  template <
      class OtherAccessorType = AccessorPolicy,
      typename                = std::enable_if_t<std::is_assignable_v<
          data_handle_type, typename OtherAccessorType::data_handle_type>>>
  KOKKOS_INLINE_FUNCTION constexpr auto to_mdspan(
      const OtherAccessorType &other_accessor) const {
    using ret_mdspan_type =
        mdspan<element_type, extents_type, layout_type, OtherAccessorType>;
    return ret_mdspan_type(
        static_cast<typename OtherAccessorType::data_handle_type>(
            data_handle()),
        mapping(), other_accessor);
  }

  KOKKOS_FUNCTION void assign_data(element_type *ptr) { m_ptr = ptr; }

  // ========================= mdspan =================================

  // [mdspan.mdspan.members], members

// Introducing the C++20 and C++23 variants of the operators already
#ifndef KOKKOS_ENABLE_CXX20
  // C++23 only operator[]
  template <class... OtherIndexTypes>
    requires((std::is_convertible_v<OtherIndexTypes, index_type> && ...) &&
             (std::is_nothrow_constructible_v<index_type, OtherIndexTypes> &&
              ...) &&
             (sizeof...(OtherIndexTypes) == rank()))
  KOKKOS_FUNCTION constexpr reference operator[](
      OtherIndexTypes... indices) const {
    return m_acc.access(m_ptr,
                        m_map(static_cast<index_type>(std::move(indices))...));
  }
#endif

  // C++20 operator()
  template <class... OtherIndexTypes>
    requires((std::is_convertible_v<OtherIndexTypes, index_type> && ...) &&
             (std::is_nothrow_constructible_v<index_type, OtherIndexTypes> &&
              ...) &&
             (sizeof...(OtherIndexTypes) == rank()))
  KOKKOS_FUNCTION constexpr reference operator()(
      OtherIndexTypes... indices) const {
    KOKKOS_IMPL_BASICVIEW_OPERATOR_VERIFY(indices...);
    return m_acc.access(m_ptr,
                        m_map(static_cast<index_type>(std::move(indices))...));
  }

#undef KOKKOS_IMPL_BASICVIEW_OPERATOR_VERIFY

 private:
  // FIXME_CXX20: could use inline templated lambda in C++20 mode inside size()
  template <size_t... Idxs>
  KOKKOS_FUNCTION constexpr size_type size_impl(
      std::index_sequence<Idxs...>) const noexcept {
    // Note we restrict data_handle to be convertible to element_type* for now.
    // This is also different from mdspan: mdspan can NOT be legally in a state
    // where m_ptr is nullptr and the product of extents is non-zero
    // The default constructor of mdspan is constrained to dynamic_rank > 0
    // For View we do not have that constraint today
    if (data_handle() == nullptr) return 0u;
    return ((static_cast<size_type>(m_map.extents().extent(Idxs))) * ... *
            size_type(1));
  }

 public:
  KOKKOS_FUNCTION constexpr size_type size() const noexcept {
    return size_impl(std::make_index_sequence<rank()>());
  }

 private:
  // FIXME_CXX20: could use inline templated lambda in C++20 mode inside empty()
  template <size_t... Idxs>
  KOKKOS_FUNCTION constexpr bool empty_impl(
      std::index_sequence<Idxs...>) const noexcept {
    // Note we restrict data_handle to be convertible to element_type* for now.
    // This is also different from mdspan: mdspan can NOT be legally in a state
    // where m_ptr is nullptr and the product of extents is non-zero
    // The default constructor of mdspan is constrained to dynamic_rank > 0
    // For View we do not have that constraint today
    if (data_handle() == nullptr) return true;
    return (rank() > 0) &&
           ((m_map.extents().extent(Idxs) == index_type(0)) || ... || false);
  }

 public:
  [[nodiscard]] KOKKOS_FUNCTION constexpr bool empty() const noexcept {
    return empty_impl(std::make_index_sequence<rank()>());
  }

  KOKKOS_FUNCTION friend constexpr void swap(BasicView &x,
                                             BasicView &y) noexcept {
    kokkos_swap(x.m_ptr, y.m_ptr);
    kokkos_swap(x.m_map, y.m_map);
    kokkos_swap(x.m_acc, y.m_acc);
  }

  KOKKOS_FUNCTION constexpr const extents_type &extents() const noexcept {
    return m_map.extents();
  }
  KOKKOS_FUNCTION constexpr const data_handle_type &data_handle()
      const noexcept {
    return m_ptr;
  }
  KOKKOS_FUNCTION constexpr const mapping_type &mapping() const noexcept {
    return m_map;
  }
  KOKKOS_FUNCTION constexpr const accessor_type &accessor() const noexcept {
    return m_acc;
  }

  KOKKOS_FUNCTION static constexpr bool is_always_unique() noexcept {
    return mapping_type::is_always_unique();
  }
  KOKKOS_FUNCTION static constexpr bool is_always_exhaustive() noexcept {
    return mapping_type::is_always_exhaustive();
  }
  KOKKOS_FUNCTION static constexpr bool is_always_strided() noexcept {
    return mapping_type::is_always_strided();
  }

  KOKKOS_FUNCTION constexpr bool is_unique() const { return m_map.is_unique(); }
  KOKKOS_FUNCTION constexpr bool is_exhaustive() const {
    return m_map.is_exhaustive();
  }
  KOKKOS_FUNCTION constexpr bool is_strided() const {
    return m_map.is_strided();
  }
  KOKKOS_FUNCTION constexpr index_type stride(rank_type r) const {
    return m_map.stride(r);
  }

 protected:
#ifndef __NVCC__
  KOKKOS_IMPL_NO_UNIQUE_ADDRESS data_handle_type m_ptr{};
  KOKKOS_IMPL_NO_UNIQUE_ADDRESS mapping_type m_map{};
  KOKKOS_IMPL_NO_UNIQUE_ADDRESS accessor_type m_acc{};
#else
  data_handle_type m_ptr{};
  mapping_type m_map{};
  accessor_type m_acc{};
#endif

  template <class, class, class, class>
  friend class BasicView;
};

#if defined(KOKKOS_COMPILER_GNU) && KOKKOS_COMPILER_GNU >= 1500
#pragma GCC diagnostic pop
#endif

}  // namespace BV
}  // namespace Kokkos::Impl

#endif
