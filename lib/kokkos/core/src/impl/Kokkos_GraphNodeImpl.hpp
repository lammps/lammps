// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_IMPL_GRAPHNODEIMPL_HPP
#define KOKKOS_IMPL_GRAPHNODEIMPL_HPP

#include <Kokkos_Macros.hpp>

#include <Kokkos_Core_fwd.hpp>
#include <Kokkos_Graph_fwd.hpp>

#include <impl/Kokkos_GraphImpl.hpp>
#include <impl/Kokkos_GraphNodeCustomization.hpp>

#include <impl/Kokkos_EBO.hpp>

#include <memory>

namespace Kokkos {
namespace Impl {

//==============================================================================
// <editor-fold desc="Fully type-erased GraphNodeImpl"> {{{1

// Base specialization for the case where both the kernel and the predecessor
// type information is type-erased
template <class ExecutionSpace>
struct GraphNodeImpl<ExecutionSpace, Kokkos::Experimental::TypeErasedTag,
                     Kokkos::Experimental::TypeErasedTag>
    : GraphNodeBackendSpecificDetails<ExecutionSpace>,
      InstanceStorage<ExecutionSpace> {
 public:
  using node_ref_t =
      Kokkos::Experimental::GraphNodeRef<ExecutionSpace,
                                         Kokkos::Experimental::TypeErasedTag,
                                         Kokkos::Experimental::TypeErasedTag>;

 protected:
  using implementation_base_t = GraphNodeBackendSpecificDetails<ExecutionSpace>;
  using execution_space_storage_base_t = InstanceStorage<ExecutionSpace>;

 public:
  virtual ~GraphNodeImpl() = default;

 protected:
  //----------------------------------------------------------------------------
  // <editor-fold desc="protected ctors and destructors"> {{{2

  explicit GraphNodeImpl(ExecutionSpace const& ex) noexcept
      : implementation_base_t(), execution_space_storage_base_t(ex) {}

  // </editor-fold> end protected ctors and destructors }}}2
  //----------------------------------------------------------------------------

 public:
  //----------------------------------------------------------------------------
  // <editor-fold desc="public(-ish) constructors"> {{{2

  template <class... Args>
  GraphNodeImpl(ExecutionSpace const& ex, _graph_node_is_root_ctor_tag,
                Args&&... args) noexcept
      : implementation_base_t(_graph_node_is_root_ctor_tag{}, (Args&&)args...),
        execution_space_storage_base_t(ex) {}

  // </editor-fold> end public(-ish) constructors }}}2
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  // <editor-fold desc="no other constructors"> {{{2

  GraphNodeImpl()                                = delete;
  GraphNodeImpl(GraphNodeImpl const&)            = delete;
  GraphNodeImpl(GraphNodeImpl&&)                 = delete;
  GraphNodeImpl& operator=(GraphNodeImpl const&) = delete;
  GraphNodeImpl& operator=(GraphNodeImpl&&)      = delete;

  // </editor-fold> end no other constructors }}}2
  //----------------------------------------------------------------------------

  ExecutionSpace const& execution_space_instance() const {
    return this->execution_space_storage_base_t::instance();
  }
};

// </editor-fold> end Fully type-erased GraphNodeImpl }}}1
//==============================================================================

//==============================================================================
// <editor-fold desc="Type-erased predecessor GraphNodeImpl"> {{{1

// Specialization for the case with the concrete type of the kernel, but the
// predecessor erased.
template <class ExecutionSpace, class Kernel>
struct GraphNodeImpl<ExecutionSpace, Kernel,
                     Kokkos::Experimental::TypeErasedTag>
    : GraphNodeImpl<ExecutionSpace, Kokkos::Experimental::TypeErasedTag,
                    Kokkos::Experimental::TypeErasedTag> {
 private:
  using base_t =
      GraphNodeImpl<ExecutionSpace, Kokkos::Experimental::TypeErasedTag,
                    Kokkos::Experimental::TypeErasedTag>;

 public:
  //----------------------------------------------------------------------------
  // <editor-fold desc="public member types"> {{{2

  using node_ref_t =
      Kokkos::Experimental::GraphNodeRef<ExecutionSpace, Kernel,
                                         Kokkos::Experimental::TypeErasedTag>;
  using kernel_type = Kernel;

  // </editor-fold> end public member types }}}2
  //----------------------------------------------------------------------------

 private:
  //----------------------------------------------------------------------------
  // <editor-fold desc="private data members"> {{{2

  Kernel m_kernel;

  // </editor-fold> end private data members }}}2
  //----------------------------------------------------------------------------

 public:
  //----------------------------------------------------------------------------
  // <editor-fold desc="Ctors, destructors, and assignment"> {{{2

  template <class KernelDeduced, class Tag,
            typename = std::enable_if_t<
                std::is_same_v<Tag, _graph_node_kernel_ctor_tag> ||
                std::is_same_v<Tag, _graph_node_capture_ctor_tag> ||
                std::is_same_v<Tag, _graph_node_host_ctor_tag>>>
  GraphNodeImpl(ExecutionSpace const& ex, Tag, KernelDeduced&& arg_kernel)
      : base_t(ex), m_kernel{(KernelDeduced&&)arg_kernel} {}

  template <class... Args>
  GraphNodeImpl(ExecutionSpace const& ex, _graph_node_is_root_ctor_tag,
                Args&&... args)
      : base_t(ex, _graph_node_is_root_ctor_tag{}, (Args&&)args...) {}

  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // <editor-fold desc="Rule of 6 for not copyable or movable"> {{{3

  // Not copyable or movable
  GraphNodeImpl()                                = delete;
  GraphNodeImpl(GraphNodeImpl const&)            = delete;
  GraphNodeImpl(GraphNodeImpl&&)                 = delete;
  GraphNodeImpl& operator=(GraphNodeImpl const&) = delete;
  GraphNodeImpl& operator=(GraphNodeImpl&&)      = delete;
  ~GraphNodeImpl() override                      = default;

  // </editor-fold> end Rule of 6 for not copyable or movable }}}3
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // </editor-fold> end Ctors, destructors, and assignment }}}2
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  // <editor-fold desc="member accessors"> {{{2

  // Reference qualified to prevent dangling reference to data member
  Kernel& get_kernel() & { return m_kernel; }
  Kernel const& get_kernel() const& { return m_kernel; }
  Kernel&& get_kernel() && = delete;

  // </editor-fold> end member accessors }}}2
  //----------------------------------------------------------------------------
};

// </editor-fold> end Type-erased predecessor GraphNodeImpl }}}1
//==============================================================================

//==============================================================================
// <editor-fold desc="Fully concrete GraphNodeImpl"> {{{1

// Specialization for the case where nothing is type-erased
template <class ExecutionSpace, class Kernel, class PredecessorRef>
struct GraphNodeImpl
    : GraphNodeImpl<ExecutionSpace, Kernel,
                    Kokkos::Experimental::TypeErasedTag>,
      GraphNodeBackendDetailsBeforeTypeErasure<ExecutionSpace, Kernel,
                                               PredecessorRef> {
 private:
  using base_t = GraphNodeImpl<ExecutionSpace, Kernel,
                               Kokkos::Experimental::TypeErasedTag>;
  using backend_details_base_t =
      GraphNodeBackendDetailsBeforeTypeErasure<ExecutionSpace, Kernel,
                                               PredecessorRef>;
  // The fully type-erased base type, for the destroy function
  using type_erased_base_t =
      GraphNodeImpl<ExecutionSpace, Kokkos::Experimental::TypeErasedTag,
                    Kokkos::Experimental::TypeErasedTag>;

 public:
  //----------------------------------------------------------------------------
  // <editor-fold desc="public data members"> {{{2

  using node_ref_t = Kokkos::Experimental::GraphNodeRef<ExecutionSpace, Kernel,
                                                        PredecessorRef>;

  // </editor-fold> end public data members }}}2
  //----------------------------------------------------------------------------

 private:
  //----------------------------------------------------------------------------
  // <editor-fold desc="private data members"> {{{2

  PredecessorRef m_predecessor_ref;

  // </editor-fold> end private data members }}}2
  //----------------------------------------------------------------------------

 public:
  //----------------------------------------------------------------------------
  // <editor-fold desc="Ctors, destructors, and assignment"> {{{2

  // Not copyable or movable
  GraphNodeImpl()                                = delete;
  GraphNodeImpl(GraphNodeImpl const&)            = delete;
  GraphNodeImpl(GraphNodeImpl&&)                 = delete;
  GraphNodeImpl& operator=(GraphNodeImpl const&) = delete;
  GraphNodeImpl& operator=(GraphNodeImpl&&)      = delete;
  ~GraphNodeImpl() override                      = default;

  // Normal kernel-and-predecessor or capture-and-predecessor constructor.
  template <class KernelDeduced, class PredecessorPtrDeduced, class Tag,
            typename = std::enable_if_t<
                std::is_same_v<Tag, _graph_node_kernel_ctor_tag> ||
                std::is_same_v<Tag, _graph_node_capture_ctor_tag> ||
                std::is_same_v<Tag, _graph_node_host_ctor_tag>>>
  GraphNodeImpl(ExecutionSpace const& ex, Tag, KernelDeduced&& arg_kernel,
                _graph_node_predecessor_ctor_tag,
                PredecessorPtrDeduced&& arg_predecessor)
      : base_t(ex, Tag{}, (KernelDeduced&&)arg_kernel),
        // The backend gets the ability to store (weak, non-owning) references
        // to the kernel in it's final resting place here if it wants. The
        // predecessor is already a pointer, so it doesn't matter that it isn't
        // already at its final address
        backend_details_base_t(ex, this->base_t::get_kernel(), arg_predecessor,
                               *this),
        m_predecessor_ref((PredecessorPtrDeduced&&)arg_predecessor) {}

  // Root-tagged constructor
  template <class... Args>
  GraphNodeImpl(ExecutionSpace const& ex, _graph_node_is_root_ctor_tag,
                Args&&... args)
      : base_t(ex, _graph_node_is_root_ctor_tag{}, (Args&&)args...),
        backend_details_base_t(ex, _graph_node_is_root_ctor_tag{}, *this),
        m_predecessor_ref() {}

  // </editor-fold> end Ctors, destructors, and assignment }}}2
  //------------------------------------------------------------------------------
};

// </editor-fold> end Fully concrete GraphNodeImpl }}}1
//==============================================================================
}  // end namespace Impl
}  // end namespace Kokkos

#endif  // KOKKOS_IMPL_GRAPHNODEIMPL_HPP
