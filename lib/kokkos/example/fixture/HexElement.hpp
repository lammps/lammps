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

#ifndef KOKKOS_HEXELEMENT_HPP
#define KOKKOS_HEXELEMENT_HPP

namespace Kokkos {
namespace Example {

template< unsigned NodeCount >
class HexElement_TensorData ;

template< unsigned NodeCount , class Device >
class HexElement_TensorEval ;

//----------------------------------------------------------------------------
/** \brief  Evaluate Hex element on interval [-1,1]^3 */
template<>
class HexElement_TensorData< 8 > {
public:

  static const unsigned element_node_count    = 8 ;
  static const unsigned spatial_dimension     = 3 ;
  static const unsigned integration_count_1d  = 2 ;
  static const unsigned function_count_1d     = 2 ;

  float values_1d [ function_count_1d ][ integration_count_1d ];
  float derivs_1d [ function_count_1d ][ integration_count_1d ];
  float weights_1d[ integration_count_1d ];

  unsigned char eval_map[ element_node_count ][4] ;

  static float eval_value_1d( const unsigned jf , const float x )
  {
    return 0 == jf ? 0.5 * ( 1.0 - x ) : (
           1 == jf ? 0.5 * ( 1.0 + x ) : 0 );
  }

  static float eval_deriv_1d( const unsigned jf , const float )
  {
    return 0 == jf ? -0.5 : (
           1 == jf ?  0.5 : 0 );
  }

  HexElement_TensorData()
  {
    const unsigned char tmp_map[ element_node_count ][ spatial_dimension ] =
      { { 0 , 0 , 0 },
        { 1 , 0 , 0 },
        { 1 , 1 , 0 },
        { 0 , 1 , 0 },
        { 0 , 0 , 1 },
        { 1 , 0 , 1 },
        { 1 , 1 , 1 },
        { 0 , 1 , 1 } };

    weights_1d[0] = 1 ;
    weights_1d[1] = 1 ;

    const float points_1d[ integration_count_1d ] =
      { -0.577350269 , 0.577350269 };

    for ( unsigned i = 0 ; i < element_node_count ; ++i ) {
      eval_map[i][0] = tmp_map[i][0];
      eval_map[i][1] = tmp_map[i][1];
      eval_map[i][2] = tmp_map[i][2];
    }

    for ( unsigned xp = 0 ; xp < integration_count_1d ; ++xp ) {
    for ( unsigned xf = 0 ; xf < function_count_1d ; ++xf ) {
      values_1d[xp][xf] = eval_value_1d( xf , points_1d[xp] );
      derivs_1d[xp][xf] = eval_deriv_1d( xf , points_1d[xp] );
    }}
  }
};

//----------------------------------------------------------------------------

template<>
class HexElement_TensorData< 27 > {
public:

  static const unsigned element_node_count    = 27 ;
  static const unsigned spatial_dimension     = 3 ;
  static const unsigned integration_count_1d  = 3 ;
  static const unsigned function_count_1d     = 3 ;

  float values_1d [ function_count_1d ][ integration_count_1d ];
  float derivs_1d [ function_count_1d ][ integration_count_1d ];
  float weights_1d[ integration_count_1d ];

  unsigned char eval_map[ element_node_count ][4] ;

  // sizeof(EvaluateElementHex) = 111 bytes =
  //   sizeof(float) * 9 +
  //   sizeof(float) * 9 +
  //   sizeof(float) * 3 +
  //   sizeof(char) * 27 

  static float eval_value_1d( const unsigned jf , const float p )
  {
    return 0 == jf ? 0.5 * p * ( p - 1 ) : (
           1 == jf ? 1.0 - p * p : (
           2 == jf ? 0.5 * p * ( p + 1 ) : 0 ));
  }

  static float eval_deriv_1d( const unsigned jf , const float p )
  {
    return 0 == jf ? p - 0.5 : (
           1 == jf ? -2.0 * p : (
           2 == jf ? p + 0.5 : 0 ));
  }

  HexElement_TensorData()
  {
    const unsigned char tmp_map[ element_node_count ][ spatial_dimension ] =
      { { 0 , 0 , 0 },
        { 2 , 0 , 0 },
        { 2 , 2 , 0 },
        { 0 , 2 , 0 },
        { 0 , 0 , 2 },
        { 2 , 0 , 2 },
        { 2 , 2 , 2 },
        { 0 , 2 , 2 },
        { 1 , 0 , 0 },
        { 2 , 1 , 0 },
        { 1 , 2 , 0 },
        { 0 , 1 , 0 },
        { 0 , 0 , 1 },
        { 2 , 0 , 1 },
        { 2 , 2 , 1 },
        { 0 , 2 , 1 },
        { 1 , 0 , 2 },
        { 2 , 1 , 2 },
        { 1 , 2 , 2 },
        { 0 , 1 , 2 },
        { 1 , 1 , 1 },
        { 1 , 1 , 0 },
        { 1 , 1 , 2 },
        { 0 , 1 , 1 },
        { 2 , 1 , 1 },
        { 1 , 0 , 1 },
        { 1 , 2 , 1 } };

    // Interval [-1,1]

    weights_1d[0] = 0.555555556 ;
    weights_1d[1] = 0.888888889 ;
    weights_1d[2] = 0.555555556 ;

    const float points_1d[3] = { -0.774596669 ,
                                  0.000000000 ,
                                  0.774596669 };

    for ( unsigned i = 0 ; i < element_node_count ; ++i ) {
      eval_map[i][0] = tmp_map[i][0];
      eval_map[i][1] = tmp_map[i][1];
      eval_map[i][2] = tmp_map[i][2];
    }

    for ( unsigned xp = 0 ; xp < integration_count_1d ; ++xp ) {
    for ( unsigned xf = 0 ; xf < function_count_1d ; ++xf ) {
      values_1d[xp][xf] = eval_value_1d( xf , points_1d[xp] );
      derivs_1d[xp][xf] = eval_deriv_1d( xf , points_1d[xp] );
    }}
  }
};

//----------------------------------------------------------------------------

template< unsigned NodeCount >
class HexElement_Data {
public:
  static const unsigned spatial_dimension   = 3 ;
  static const unsigned element_node_count  = NodeCount ;
  static const unsigned integration_count   = NodeCount ;
  static const unsigned function_count      = NodeCount ;

  float weights[   integration_count ] ;
  float values[    integration_count ][ function_count ];
  float gradients[ integration_count ][ spatial_dimension ][ function_count ];

  HexElement_Data()
  {
    HexElement_TensorData< NodeCount > tensor_data ;

    for ( unsigned ip = 0 ; ip < integration_count ; ++ip ) {

      const unsigned ipx = tensor_data.eval_map[ip][0] ;
      const unsigned ipy = tensor_data.eval_map[ip][1] ;
      const unsigned ipz = tensor_data.eval_map[ip][2] ;

      weights[ip] = tensor_data.weights_1d[ ipx ] *
                    tensor_data.weights_1d[ ipy ] *
                    tensor_data.weights_1d[ ipz ] ;

      for ( unsigned jf = 0 ; jf < function_count ; ++jf ) {

        const unsigned jfx = tensor_data.eval_map[jf][0] ;
        const unsigned jfy = tensor_data.eval_map[jf][1] ;
        const unsigned jfz = tensor_data.eval_map[jf][2] ;

        values[ip][jf] = tensor_data.values_1d[ ipx ][ jfx ] *
                         tensor_data.values_1d[ ipy ][ jfy ] *
                         tensor_data.values_1d[ ipz ][ jfz ] ;

        gradients[ip][0][jf] = tensor_data.derivs_1d[ ipx ][ jfx ] *
                               tensor_data.values_1d[ ipy ][ jfy ] *
                               tensor_data.values_1d[ ipz ][ jfz ] ;

        gradients[ip][1][jf] = tensor_data.values_1d[ ipx ][ jfx ] *
                               tensor_data.derivs_1d[ ipy ][ jfy ] *
                               tensor_data.values_1d[ ipz ][ jfz ] ;

        gradients[ip][2][jf] = tensor_data.values_1d[ ipx ][ jfx ] *
                               tensor_data.values_1d[ ipy ][ jfy ] *
                               tensor_data.derivs_1d[ ipz ][ jfz ] ;
      }
    }
  }
};

//----------------------------------------------------------------------------

} /* namespace Example */
} /* namespace Kokkos */

#endif /* #ifndef KOKKOS_HEXELEMENT_HPP */


