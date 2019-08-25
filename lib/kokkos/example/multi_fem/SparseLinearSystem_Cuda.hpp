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

#ifndef SPARSELINEARSYSTEM_CUDA_HPP
#define SPARSELINEARSYSTEM_CUDA_HPP

#if defined( BUILD_FROM_CU_FILE )

#include <cusparse_v2.h>
#include <Kokkos_Core.hpp>

namespace Kokkos {
namespace Impl {


struct CudaSparseSingleton {
  cusparseHandle_t   handle;
  cusparseMatDescr_t descra;

  CudaSparseSingleton()
  {
    cusparseCreate( & handle );
    cusparseCreateMatDescr( & descra );
    cusparseSetMatType(       descra , CUSPARSE_MATRIX_TYPE_GENERAL );
    cusparseSetMatIndexBase(  descra , CUSPARSE_INDEX_BASE_ZERO );
  }

  static CudaSparseSingleton & singleton();

};

CudaSparseSingleton & CudaSparseSingleton::singleton()
{ static CudaSparseSingleton s ; return s ; }


template<>
struct Multiply< CrsMatrix<double,Cuda> ,
                 View<double*,Cuda > ,
                 View<double*,Cuda > >
{
  typedef Cuda                                      execution_space ;
  typedef execution_space::size_type                    size_type ;
  typedef double                                    scalar_type ;
  typedef View< scalar_type* , execution_space >        vector_type ;
  typedef CrsMatrix< scalar_type , execution_space >    matrix_type ;

public:

  Multiply( const matrix_type & A ,
            const size_type nrow ,
            const size_type ncol ,
            const vector_type & x ,
            const vector_type & y )
  {
    CudaSparseSingleton & s = CudaSparseSingleton::singleton();
    const scalar_type alpha = 1 , beta = 0 ;

    cusparseStatus_t status =
      cusparseDcsrmv( s.handle ,
                      CUSPARSE_OPERATION_NON_TRANSPOSE ,
                      nrow , ncol , A.coefficients.dimension_0() ,
                      &alpha ,
                      s.descra ,
                      A.coefficients.ptr_on_device() ,
                      A.graph.row_map.ptr_on_device() ,
                      A.graph.entries.ptr_on_device() ,
                      x.ptr_on_device() ,
                      &beta ,
                      y.ptr_on_device() );

    if ( CUSPARSE_STATUS_SUCCESS != status ) {
      throw std::runtime_error( std::string("ERROR - cusparseDcsrmv " ) );
    }
  }
};


template<>
struct Multiply< CrsMatrix<float,Cuda> ,
                 View<float*,Cuda > ,
                 View<float*,Cuda > >
{
  typedef Cuda                                      execution_space ;
  typedef execution_space::size_type                    size_type ;
  typedef float                                     scalar_type ;
  typedef View< scalar_type* , execution_space >        vector_type ;
  typedef CrsMatrix< scalar_type , execution_space >    matrix_type ;

public:

  Multiply( const matrix_type & A ,
            const size_type nrow ,
            const size_type ncol ,
            const vector_type & x ,
            const vector_type & y )
  {
    CudaSparseSingleton & s = CudaSparseSingleton::singleton();
    const scalar_type alpha = 1 , beta = 0 ;

    cusparseStatus_t status =
      cusparseScsrmv( s.handle ,
                      CUSPARSE_OPERATION_NON_TRANSPOSE ,
                      nrow , ncol , A.coefficients.dimension_0() ,
                      &alpha ,
                      s.descra ,
                      A.coefficients.ptr_on_device() ,
                      A.graph.row_map.ptr_on_device() ,
                      A.graph.entries.ptr_on_device() ,
                      x.ptr_on_device() ,
                      &beta ,
                      y.ptr_on_device() );

    if ( CUSPARSE_STATUS_SUCCESS != status ) {
      throw std::runtime_error( std::string("ERROR - cusparseDcsrmv " ) );
    }
  }
};

} /* namespace Impl */
} /* namespace Kokkos */

#endif /* #if defined( __CUDACC__ ) */
#endif /* #ifndef SPARSELINEARSYSTEM_CUDA_HPP */

