/*
//@HEADER
// ************************************************************************
// 
//   Kokkos: Manycore Performance-Portable Multidimensional Arrays
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
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOS_CUDATYPES_HPP
#define KOKKOS_CUDATYPES_HPP

#include <Kokkos_Macros.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#if defined( __CUDACC__ )

namespace Kokkos {

typedef ::int2 int2 ;
typedef ::int3 int3 ;
typedef ::int4 int4 ;

typedef ::float2 float2 ;
typedef ::float3 float3 ;
typedef ::float4 float4 ;

typedef ::double2 double2 ;
typedef ::double3 double3 ;
typedef ::double4 double4 ;

} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#else /* NOT #if defined( __CUDACC__ ) */

namespace Kokkos {

struct int2 {
        int x;
        int y;
};

struct int3 {
        int x;
        int y;
        int z;
};

struct int4 {
        int x;
        int y;
        int z;
        int w;
};

struct float2 {
        float x;
        float y;
};

struct float3 {
        float x;
        float y;
        float z;
};

struct float4 {
        float x;
        float y;
        float z;
        float w;
};

struct double2 {
        double x;
        double y;
};

struct double3 {
        double x;
        double y;
        double z;
};

struct double4 {
        double x;
        double y;
        double z;
        double w;
};

} // namespace Kokkos

#endif

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #define KOKKOS_CUDATYPES_HPP */

