// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_IMPL_KOKKOS_GRAPHNODECTORPROPS_HPP
#define KOKKOS_IMPL_KOKKOS_GRAPHNODECTORPROPS_HPP

#include "impl/Kokkos_DeviceHandle.hpp"
#include "View/Kokkos_ViewCtor.hpp"

namespace Kokkos::Impl {

template <typename T>
concept ValidNodeProperty = is_view_label_v<T> || is_device_handle_v<T>;

// Transform a type into a valid property.
template <typename T>
struct NodeCtorPropTransform {
  using type = T;
};

template <Kokkos::Impl::ViewLabel T>
struct NodeCtorPropTransform<T> {
  using type = std::string;
};

// A single property.
template <ValidNodeProperty T>
struct NodeCtorProp {
  using value_type = typename NodeCtorPropTransform<T>::type;

  explicit NodeCtorProp(value_type value) : m_value(std::move(value)) {}

  value_type m_value;
};

// Aggregated properties.
template <typename... Props>
struct NodeCtorProps : public NodeCtorProp<Props>... {
  using properties_value_type_list_t =
      Kokkos::Impl::type_list<typename NodeCtorProp<Props>::value_type...>;

  static_assert(Kokkos::Impl::type_list_size_v<Kokkos::Impl::filter_type_list_t<
                        is_view_label, Kokkos::Impl::type_list<Props...>>> <= 1,
                "Only one label allowed.");
  static_assert(Kokkos::Impl::type_list_size_v<Kokkos::Impl::filter_type_list_t<
                        is_device_handle, Kokkos::Impl::type_list<Props...>>> <=
                    1,
                "Only one device handle allowed.");

  using uniform_type =
      NodeCtorProps<typename NodeCtorProp<Props>::value_type...>;

  template <typename T>
  static constexpr bool has =
      Kokkos::Impl::type_list_contains_v<T, properties_value_type_list_t>;

  NodeCtorProps() = default;

  // NOLINTBEGIN(modernize-type-traits)
  template <typename... Args>
    requires(std::constructible_from<NodeCtorProp<Props>, Args &&> && ...)
  // NOLINTEND(modernize-type-traits)
  explicit NodeCtorProps(Args&&... args)
      : NodeCtorProp<Props>{std::forward<Args>(args)}... {}
};

template <typename T>
struct is_node_props : public std::false_type {};

template <typename... Args>
struct is_node_props<NodeCtorProps<Args...>> : public std::true_type {};

template <typename T>
constexpr bool is_node_props_v = is_node_props<T>::value;

template <typename T>
concept NodeProperties = is_node_props_v<T>;

template <typename Prop, NodeProperties Props>
  requires Props::template
has<Prop> [[nodiscard]] constexpr decltype(auto) get_property(
    const Props& props) {
  return static_cast<const NodeCtorProp<Prop>&>(props).m_value;
}

template <typename Prop, NodeProperties Props>
  requires Props::template
has<Prop> [[nodiscard]] constexpr decltype(auto) extract_property(
    Props& props) {
  return std::move(static_cast<NodeCtorProp<Prop>&>(props).m_value);
}

struct WithProperty {
  template <typename... Props, typename Property>
    requires(sizeof...(Props) > 0)
  [[nodiscard]] static constexpr decltype(auto) set(
      NodeCtorProps<Props...> props, Property&& property) {
    using NewNodeCtorProps =
        typename NodeCtorProps<Props...,
                               std::remove_cvref_t<Property>>::uniform_type;
    return NewNodeCtorProps{extract_property<Props>(props)...,
                            std::forward<Property>(property)};
  }

  template <typename Property>
  [[nodiscard]] static constexpr decltype(auto) set(NodeCtorProps<>,
                                                    Property&& prop) {
    return typename NodeCtorProps<std::remove_cvref_t<Property>>::uniform_type{
        std::forward<Property>(prop)};
  }
};

template <NodeProperties Props>
[[nodiscard]] constexpr decltype(auto) with_properties_if_unset(
    Props node_props) noexcept {
  return node_props;
}

template <NodeProperties Props, typename Property, typename... Properties>
[[nodiscard]] constexpr decltype(auto) with_properties_if_unset(
    Props node_props, [[maybe_unused]] Property&& property,
    Properties&&... properties) {
  if constexpr (!Props::template has<typename Kokkos::Impl::NodeCtorProp<
                    std::remove_cvref_t<Property>>::value_type>) {
    return with_properties_if_unset(
        WithProperty::set(std::move(node_props),
                          std::forward<Property>(property)),
        std::forward<Properties>(properties)...);
  } else {
    return with_properties_if_unset(std::move(node_props),
                                    std::forward<Properties>(properties)...);
  }
}

}  // namespace Kokkos::Impl

namespace Kokkos::Experimental {
template <typename... Args>
[[nodiscard]] constexpr auto node_props(Args&&... args) {
  using return_t = typename Kokkos::Impl::NodeCtorProps<
      std::remove_cvref_t<Args>...>::uniform_type;
  return return_t{std::forward<Args>(args)...};
}

}  // namespace Kokkos::Experimental

#endif  // KOKKOS_IMPL_KOKKOS_GRAPHNODECTORPROPS_HPP
