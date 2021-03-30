/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
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
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
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

#ifndef KOKKOS_IMPL_GRAPHNODEIMPL_HPP
#define KOKKOS_IMPL_GRAPHNODEIMPL_HPP

#include <Kokkos_Macros.hpp>

#include <Kokkos_Core_fwd.hpp>
#include <Kokkos_Graph_fwd.hpp>

#include <impl/Kokkos_SimpleTaskScheduler.hpp>  // ExecutionSpaceInstanceStorage
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
      ExecutionSpaceInstanceStorage<ExecutionSpace> {
 public:
  using node_ref_t =
      Kokkos::Experimental::GraphNodeRef<ExecutionSpace,
                                         Kokkos::Experimental::TypeErasedTag,
                                         Kokkos::Experimental::TypeErasedTag>;

 protected:
  using implementation_base_t = GraphNodeBackendSpecificDetails<ExecutionSpace>;
  using execution_space_storage_base_t =
      ExecutionSpaceInstanceStorage<ExecutionSpace>;

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
      : implementation_base_t(_graph_node_is_root_ctor_tag{},
                              (Args &&) args...),
        execution_space_storage_base_t(ex) {}

  // </editor-fold> end public(-ish) constructors }}}2
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  // <editor-fold desc="no other constructors"> {{{2

  GraphNodeImpl()                     = delete;
  GraphNodeImpl(GraphNodeImpl const&) = delete;
  GraphNodeImpl(GraphNodeImpl&&)      = delete;
  GraphNodeImpl& operator=(GraphNodeImpl const&) = delete;
  GraphNodeImpl& operator=(GraphNodeImpl&&) = delete;

  // </editor-fold> end no other constructors }}}2
  //----------------------------------------------------------------------------

  ExecutionSpace const& execution_space_instance() const {
    return this->execution_space_storage_base_t::execution_space_instance();
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

  template <class KernelDeduced>
  GraphNodeImpl(ExecutionSpace const& ex, _graph_node_kernel_ctor_tag,
                KernelDeduced&& arg_kernel)
      : base_t(ex), m_kernel((KernelDeduced &&) arg_kernel) {}

  template <class... Args>
  GraphNodeImpl(ExecutionSpace const& ex, _graph_node_is_root_ctor_tag,
                Args&&... args)
      : base_t(ex, _graph_node_is_root_ctor_tag{}, (Args &&) args...) {}

  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // <editor-fold desc="Rule of 6 for not copyable or movable"> {{{3

  // Not copyable or movable
  GraphNodeImpl()                     = delete;
  GraphNodeImpl(GraphNodeImpl const&) = delete;
  GraphNodeImpl(GraphNodeImpl&&)      = delete;
  GraphNodeImpl& operator=(GraphNodeImpl const&) = delete;
  GraphNodeImpl& operator=(GraphNodeImpl&&) = delete;
  ~GraphNodeImpl() override                 = default;

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
  GraphNodeImpl()                     = delete;
  GraphNodeImpl(GraphNodeImpl const&) = delete;
  GraphNodeImpl(GraphNodeImpl&&)      = delete;
  GraphNodeImpl& operator=(GraphNodeImpl const&) = delete;
  GraphNodeImpl& operator=(GraphNodeImpl&&) = delete;
  ~GraphNodeImpl() override                 = default;

  // Normal kernel-and-predecessor constructor
  template <class KernelDeduced, class PredecessorPtrDeduced>
  GraphNodeImpl(ExecutionSpace const& ex, _graph_node_kernel_ctor_tag,
                KernelDeduced&& arg_kernel, _graph_node_predecessor_ctor_tag,
                PredecessorPtrDeduced&& arg_predecessor)
      : base_t(ex, _graph_node_kernel_ctor_tag{},
               (KernelDeduced &&) arg_kernel),
        // The backend gets the ability to store (weak, non-owning) references
        // to the kernel in it's final resting place here if it wants. The
        // predecessor is already a pointer, so it doesn't matter that it isn't
        // already at its final address
        backend_details_base_t(ex, this->base_t::get_kernel(), arg_predecessor,
                               *this),
        m_predecessor_ref((PredecessorPtrDeduced &&) arg_predecessor) {}

  // Root-tagged constructor
  template <class... Args>
  GraphNodeImpl(ExecutionSpace const& ex, _graph_node_is_root_ctor_tag,
                Args&&... args)
      : base_t(ex, _graph_node_is_root_ctor_tag{}, (Args &&) args...),
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
