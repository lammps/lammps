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

// This example simulates one timestep of an explicit
// finite-difference discretization of a time-dependent partial
// differential equation (PDE).  It shows how to take subviews of the
// mesh in order to represent particular boundaries or the interior of
// the mesh.

#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>
#include <cstdio>

using mesh_type = Kokkos::View<double***, Kokkos::LayoutRight>;

// These View types represent subviews of the mesh.  Some of the Views
// have layout LayoutStride, meaning that they have run-time "strides"
// in each dimension which may differ from that dimension.  For
// example, inner_mesh_type (which represents the interior of the
// mesh) has to skip over the boundaries when computing its stride;
// the dimensions of the interior mesh differ from these strides.  You
// may safely always use a LayoutStride layout when taking a subview
// of a LayoutRight or LayoutLeft subview, but strided accesses may
// cost a bit more, especially for 1-D Views.
using xz_plane_type   = Kokkos::View<double**, Kokkos::LayoutStride>;
using yz_plane_type   = Kokkos::View<double**, Kokkos::LayoutRight>;
using xy_plane_type   = Kokkos::View<double**, Kokkos::LayoutStride>;
using inner_mesh_type = Kokkos::View<double***, Kokkos::LayoutStride>;

// Functor to set all entries of a boundary of the mesh to a constant
// value.  The functor is templated on ViewType because different
// boundaries may have different layouts.
template <class ViewType>
struct set_boundary {
  ViewType a;
  double value;

  set_boundary(ViewType a_, double value_) : a(a_), value(value_) {}

  using size_type = typename ViewType::size_type;

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_type i) const {
    for (size_type j = 0; j < static_cast<size_type>(a.extent(1)); ++j) {
      a(i, j) = value;
    }
  }
};

// Functor to set all entries of a boundary of the mesh to a constant
// value.  The functor is templated on ViewType because different
// boundaries may have different layouts.
template <class ViewType>
struct set_inner {
  ViewType a;
  double value;

  set_inner(ViewType a_, double value_) : a(a_), value(value_) {}

  using size_type = typename ViewType::size_type;

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_type i) const {
    for (size_type j = 0; j < static_cast<size_type>(a.extent(1)); ++j) {
      for (size_type k = 0; k < static_cast<size_type>(a.extent(2)); ++k) {
        a(i, j, k) = value;
      }
    }
  }
};

// Update the interior of the mesh.  This simulates one timestep of a
// finite-difference method.
template <class ViewType>
struct update {
  ViewType a;
  const double dt;

  update(ViewType a_, const double dt_) : a(a_), dt(dt_) {}

  using size_type = typename ViewType::size_type;

  KOKKOS_INLINE_FUNCTION
  void operator()(size_type i) const {
    i++;
    for (size_type j = 1; j < static_cast<size_type>(a.extent(1) - 1); j++) {
      for (size_type k = 1; k < static_cast<size_type>(a.extent(2) - 1); k++) {
        a(i, j, k) += dt * (a(i, j, k + 1) - a(i, j, k - 1) + a(i, j + 1, k) -
                            a(i, j - 1, k) + a(i + 1, j, k) - a(i - 1, j, k));
      }
    }
  }
};

int main(int narg, char* arg[]) {
  using Kokkos::ALL;
  using Kokkos::pair;
  using Kokkos::parallel_for;
  using Kokkos::subview;
  using size_type = mesh_type::size_type;

  Kokkos::initialize(narg, arg);

  {
    // The number of mesh points along each dimension of the mesh, not
    // including boundaries.
    const size_type size = 100;

    // A is the full cubic 3-D mesh, including the boundaries.
    mesh_type A("A", size + 2, size + 2, size + 2);
    // Ai is the "inner" part of A, _not_ including the boundaries.
    //
    // A pair of indices in a particular dimension means the contiguous
    // zero-based index range in that dimension, including the first
    // entry of the pair but _not_ including the second entry.
    inner_mesh_type Ai = subview(A, pair<size_type, size_type>(1, size + 1),
                                 pair<size_type, size_type>(1, size + 1),
                                 pair<size_type, size_type>(1, size + 1));
    // A has six boundaries, one for each face of the cube.
    // Create a View of each of these boundaries.
    // ALL() means "select all indices in that dimension."
    xy_plane_type Zneg_halo = subview(A, ALL(), ALL(), 0);
    xy_plane_type Zpos_halo = subview(A, ALL(), ALL(), 101);
    xz_plane_type Yneg_halo = subview(A, ALL(), 0, ALL());
    xz_plane_type Ypos_halo = subview(A, ALL(), 101, ALL());
    yz_plane_type Xneg_halo = subview(A, 0, ALL(), ALL());
    yz_plane_type Xpos_halo = subview(A, 101, ALL(), ALL());

    // Set the boundaries to their initial conditions.
    parallel_for(Zneg_halo.extent(0),
                 set_boundary<xy_plane_type>(Zneg_halo, 1));
    parallel_for(Zpos_halo.extent(0),
                 set_boundary<xy_plane_type>(Zpos_halo, -1));
    parallel_for(Yneg_halo.extent(0),
                 set_boundary<xz_plane_type>(Yneg_halo, 2));
    parallel_for(Ypos_halo.extent(0),
                 set_boundary<xz_plane_type>(Ypos_halo, -2));
    parallel_for(Xneg_halo.extent(0),
                 set_boundary<yz_plane_type>(Xneg_halo, 3));
    parallel_for(Xpos_halo.extent(0),
                 set_boundary<yz_plane_type>(Xpos_halo, -3));

    // Set the interior of the mesh to its initial condition.
    parallel_for(Ai.extent(0), set_inner<inner_mesh_type>(Ai, 0));

    // Update the interior of the mesh.
    // This simulates one timestep with dt = 0.1.
    parallel_for(Ai.extent(0), update<mesh_type>(A, 0.1));

    Kokkos::fence();
    printf("Done\n");
  }
  Kokkos::finalize();
}
