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
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov)
// 
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOS_FEMESH_HPP
#define KOKKOS_FEMESH_HPP

#include <utility>
#include <limits>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include <Kokkos_Core.hpp>
#include <Kokkos_StaticCrsGraph.hpp>

#include <ParallelComm.hpp>
#include <ParallelDataMap.hpp>

namespace HybridFEM {

//----------------------------------------------------------------------------
/** \brief  Finite element mesh fixture for hybrid parallel performance tests.
 */
template< typename CoordScalarType , unsigned ElemNodeCount , class Device >
struct FEMesh {

  typedef typename Device::size_type size_type ;

  static const size_type element_node_count = ElemNodeCount ;

  typedef Kokkos::View< CoordScalarType*[3] , Device >       node_coords_type ;
  typedef Kokkos::View< size_type*[ElemNodeCount], Device >  elem_node_ids_type ;
  typedef Kokkos::StaticCrsGraph< size_type[2] ,  Device >   node_elem_ids_type ;

  node_coords_type         node_coords ;
  elem_node_ids_type       elem_node_ids ;
  node_elem_ids_type       node_elem_ids ;
  Kokkos::ParallelDataMap  parallel_data_map ;
};

//----------------------------------------------------------------------------

} /* namespace HybridFEM */

#endif /* #ifndef KOKKOS_FEMESH_HPP */

