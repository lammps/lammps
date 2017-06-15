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

#include <Kokkos_Core.hpp>

namespace TestAtomicOperations {

//-----------------------------------------------
//--------------zero_functor---------------------
//-----------------------------------------------

template< class T, class DEVICE_TYPE >
struct ZeroFunctor {
  typedef DEVICE_TYPE execution_space;
  typedef typename Kokkos::View< T, execution_space > type;
  typedef typename Kokkos::View< T, execution_space >::HostMirror h_type;

  type data;

  KOKKOS_INLINE_FUNCTION
  void operator()( int ) const {
    data() = 0;
  }
};

//-----------------------------------------------
//--------------init_functor---------------------
//-----------------------------------------------

template< class T, class DEVICE_TYPE >
struct InitFunctor {
  typedef DEVICE_TYPE execution_space;
  typedef typename Kokkos::View< T, execution_space > type;
  typedef typename Kokkos::View< T, execution_space >::HostMirror h_type;

  type data;
  T init_value;

  KOKKOS_INLINE_FUNCTION
  void operator()( int ) const {
    data() = init_value;
  }

  InitFunctor( T _init_value ) : init_value( _init_value ) {}
};

//---------------------------------------------------
//--------------atomic_fetch_max---------------------
//---------------------------------------------------

template< class T, class DEVICE_TYPE >
struct MaxFunctor {
  typedef DEVICE_TYPE execution_space;
  typedef Kokkos::View< T, execution_space > type;

  type data;
  T i0;
  T i1;

  KOKKOS_INLINE_FUNCTION
  void operator()( int ) const {
    //Kokkos::atomic_fetch_max( &data(), (T) 1 );
    Kokkos::atomic_fetch_max( &data(), (T) i1 );
  }
  MaxFunctor( T _i0, T _i1 ) : i0( _i0 ), i1( _i1 ) {}
};

template< class T, class execution_space >
T MaxAtomic( T i0, T i1 ) {
  struct InitFunctor< T, execution_space > f_init( i0 );
  typename InitFunctor< T, execution_space >::type data( "Data" );
  typename InitFunctor< T, execution_space >::h_type h_data( "HData" );

  f_init.data = data;
  Kokkos::parallel_for( 1, f_init );
  execution_space::fence();

  struct MaxFunctor< T, execution_space > f( i0, i1 );

  f.data = data;
  Kokkos::parallel_for( 1, f );
  execution_space::fence();

  Kokkos::deep_copy( h_data, data );
  T val = h_data();

  return val;
}

template< class T >
T MaxAtomicCheck( T i0, T i1 ) {
  T* data = new T[1];
  data[0] = 0;

  *data = ( i0 > i1 ? i0 : i1 );

  T val = *data;
  delete [] data;

  return val;
}

template< class T, class DeviceType >
bool MaxAtomicTest( T i0, T i1 )
{
  T res       = MaxAtomic< T, DeviceType >( i0, i1 );
  T resSerial = MaxAtomicCheck<T>( i0, i1 );

  bool passed = true;

  if ( resSerial != res ) {
    passed = false;

    std::cout << "Loop<"
              << typeid( T ).name()
              << ">( test = MaxAtomicTest"
              << " FAILED : "
              << resSerial << " != " << res
              << std::endl;
  }

  return passed;
}

//---------------------------------------------------
//--------------atomic_fetch_min---------------------
//---------------------------------------------------

template< class T, class DEVICE_TYPE >
struct MinFunctor {
  typedef DEVICE_TYPE execution_space;
  typedef Kokkos::View< T, execution_space > type;

  type data;
  T i0;
  T i1;

  KOKKOS_INLINE_FUNCTION
  void operator()( int ) const {
    Kokkos::atomic_fetch_min( &data(), (T) i1 );
  }

  MinFunctor( T _i0, T _i1 ) : i0( _i0 ), i1( _i1 ) {}
};

template< class T, class execution_space >
T MinAtomic( T i0, T i1 ) {
  struct InitFunctor< T, execution_space > f_init( i0 );
  typename InitFunctor< T, execution_space >::type data( "Data" );
  typename InitFunctor< T, execution_space >::h_type h_data( "HData" );

  f_init.data = data;
  Kokkos::parallel_for( 1, f_init );
  execution_space::fence();

  struct MinFunctor< T, execution_space > f( i0, i1 );

  f.data = data;
  Kokkos::parallel_for( 1, f );
  execution_space::fence();

  Kokkos::deep_copy( h_data, data );
  T val = h_data();

  return val;
}

template< class T >
T MinAtomicCheck( T i0, T i1 ) {
  T* data = new T[1];
  data[0] = 0;

  *data = ( i0 < i1 ? i0 : i1 );

  T val = *data;
  delete [] data;

  return val;
}

template< class T, class DeviceType >
bool MinAtomicTest( T i0, T i1 )
{
  T res       = MinAtomic< T, DeviceType >( i0, i1 );
  T resSerial = MinAtomicCheck< T >( i0, i1 );

  bool passed = true;

  if ( resSerial != res ) {
    passed = false;

    std::cout << "Loop<"
              << typeid( T ).name()
              << ">( test = MinAtomicTest"
              << " FAILED : "
              << resSerial << " != " << res
              << std::endl;
  }

  return passed;
}

//---------------------------------------------------
//--------------atomic_increment---------------------
//---------------------------------------------------

template< class T, class DEVICE_TYPE >
struct IncFunctor {
  typedef DEVICE_TYPE execution_space;
  typedef Kokkos::View< T, execution_space > type;

  type data;
  T i0;

  KOKKOS_INLINE_FUNCTION
  void operator()( int ) const {
    Kokkos::atomic_increment( &data() );
  }

  IncFunctor( T _i0 ) : i0( _i0 ) {}
};

template< class T, class execution_space >
T IncAtomic( T i0 ) {
  struct InitFunctor< T, execution_space > f_init( i0 );
  typename InitFunctor< T, execution_space >::type data( "Data" );
  typename InitFunctor< T, execution_space >::h_type h_data( "HData" );

  f_init.data = data;
  Kokkos::parallel_for( 1, f_init );
  execution_space::fence();

  struct IncFunctor< T, execution_space > f( i0 );

  f.data = data;
  Kokkos::parallel_for( 1, f );
  execution_space::fence();

  Kokkos::deep_copy( h_data, data );
  T val = h_data();

  return val;
}

template< class T >
T IncAtomicCheck( T i0 ) {
  T* data = new T[1];
  data[0] = 0;

  *data = i0 + 1;

  T val = *data;
  delete [] data;

  return val;
}

template< class T, class DeviceType >
bool IncAtomicTest( T i0 )
{
  T res       = IncAtomic< T, DeviceType >( i0 );
  T resSerial = IncAtomicCheck< T >( i0 );

  bool passed = true;

  if ( resSerial != res ) {
    passed = false;

    std::cout << "Loop<"
              << typeid( T ).name()
              << ">( test = IncAtomicTest"
              << " FAILED : "
              << resSerial << " != " << res
              << std::endl;
  }

  return passed;
}

//---------------------------------------------------
//--------------atomic_decrement---------------------
//---------------------------------------------------

template< class T, class DEVICE_TYPE >
struct DecFunctor {
  typedef DEVICE_TYPE execution_space;
  typedef Kokkos::View< T, execution_space > type;

  type data;
  T i0;

  KOKKOS_INLINE_FUNCTION
  void operator()( int ) const {
    Kokkos::atomic_decrement( &data() );
  }

  DecFunctor( T _i0 ) : i0( _i0 ) {}
};

template< class T, class execution_space >
T DecAtomic( T i0 ) {
  struct InitFunctor< T, execution_space > f_init( i0 );
  typename InitFunctor< T, execution_space >::type data( "Data" );
  typename InitFunctor< T, execution_space >::h_type h_data( "HData" );

  f_init.data = data;
  Kokkos::parallel_for( 1, f_init );
  execution_space::fence();

  struct DecFunctor< T, execution_space > f( i0 );

  f.data = data;
  Kokkos::parallel_for( 1, f );
  execution_space::fence();

  Kokkos::deep_copy( h_data, data );
  T val = h_data();

  return val;
}

template< class T >
T DecAtomicCheck( T i0 ) {
  T* data = new T[1];
  data[0] = 0;

  *data = i0 - 1;

  T val = *data;
  delete [] data;

  return val;
}

template< class T, class DeviceType >
bool DecAtomicTest( T i0 )
{
  T res       = DecAtomic< T, DeviceType >( i0 );
  T resSerial = DecAtomicCheck< T >( i0 );

  bool passed = true;

  if ( resSerial != res ) {
    passed = false;

    std::cout << "Loop<"
              << typeid( T ).name()
              << ">( test = DecAtomicTest"
              << " FAILED : "
              << resSerial << " != " << res
              << std::endl;
  }

  return passed;
}

//---------------------------------------------------
//--------------atomic_fetch_mul---------------------
//---------------------------------------------------

template< class T, class DEVICE_TYPE >
struct MulFunctor {
  typedef DEVICE_TYPE execution_space;
  typedef Kokkos::View< T, execution_space > type;

  type data;
  T i0;
  T i1;

  KOKKOS_INLINE_FUNCTION
  void operator()( int ) const {
    Kokkos::atomic_fetch_mul( &data(), (T) i1 );
  }

  MulFunctor( T _i0, T _i1 ) : i0( _i0 ), i1( _i1 ) {}
};

template< class T, class execution_space >
T MulAtomic( T i0, T i1 ) {
  struct InitFunctor< T, execution_space > f_init( i0 );
  typename InitFunctor< T, execution_space >::type data( "Data" );
  typename InitFunctor< T, execution_space >::h_type h_data( "HData" );

  f_init.data = data;
  Kokkos::parallel_for( 1, f_init );
  execution_space::fence();

  struct MulFunctor< T, execution_space > f( i0, i1 );

  f.data = data;
  Kokkos::parallel_for( 1, f );
  execution_space::fence();

  Kokkos::deep_copy( h_data, data );
  T val = h_data();

  return val;
}

template< class T >
T MulAtomicCheck( T i0, T i1 ) {
  T* data = new T[1];
  data[0] = 0;

  *data = i0*i1;

  T val = *data;
  delete [] data;

  return val;
}

template< class T, class DeviceType >
bool MulAtomicTest( T i0, T i1 )
{
  T res       = MulAtomic< T, DeviceType >( i0, i1 );
  T resSerial = MulAtomicCheck< T >( i0, i1 );

  bool passed = true;

  if ( resSerial != res ) {
    passed = false;

    std::cout << "Loop<"
              << typeid( T ).name()
              << ">( test = MulAtomicTest"
              << " FAILED : "
              << resSerial << " != " << res
              << std::endl;
  }

  return passed;
}

//---------------------------------------------------
//--------------atomic_fetch_div---------------------
//---------------------------------------------------

template< class T, class DEVICE_TYPE >
struct DivFunctor {
  typedef DEVICE_TYPE execution_space;
  typedef Kokkos::View< T, execution_space > type;

  type data;
  T i0;
  T i1;

  KOKKOS_INLINE_FUNCTION
  void operator()( int ) const {
    Kokkos::atomic_fetch_div( &data(), (T) i1 );
  }

  DivFunctor( T _i0, T _i1 ) : i0( _i0 ), i1( _i1 ) {}
};

template< class T, class execution_space >
T DivAtomic( T i0, T i1 ) {
  struct InitFunctor< T, execution_space > f_init( i0 );
  typename InitFunctor< T, execution_space >::type data( "Data" );
  typename InitFunctor< T, execution_space >::h_type h_data( "HData" );

  f_init.data = data;
  Kokkos::parallel_for( 1, f_init );
  execution_space::fence();

  struct DivFunctor< T, execution_space > f( i0, i1 );

  f.data = data;
  Kokkos::parallel_for( 1, f );
  execution_space::fence();

  Kokkos::deep_copy( h_data, data );
  T val = h_data();

  return val;
}

template< class T >
T DivAtomicCheck( T i0, T i1 ) {
  T* data = new T[1];
  data[0] = 0;

  *data = i0 / i1;

  T val = *data;
  delete [] data;

  return val;
}

template< class T, class DeviceType >
bool DivAtomicTest( T i0, T i1 )
{
  T res       = DivAtomic< T, DeviceType >( i0, i1 );
  T resSerial = DivAtomicCheck< T >( i0, i1 );

  bool passed = true;

  if ( resSerial != res ) {
    passed = false;

    std::cout << "Loop<"
              << typeid( T ).name()
              << ">( test = DivAtomicTest"
              << " FAILED : "
              << resSerial << " != " << res
              << std::endl;
  }

  return passed;
}

//---------------------------------------------------
//--------------atomic_fetch_mod---------------------
//---------------------------------------------------

template< class T, class DEVICE_TYPE >
struct ModFunctor {
  typedef DEVICE_TYPE execution_space;
  typedef Kokkos::View< T, execution_space > type;

  type data;
  T i0;
  T i1;

  KOKKOS_INLINE_FUNCTION
  void operator()( int ) const {
    Kokkos::atomic_fetch_mod( &data(), (T) i1 );
  }

  ModFunctor( T _i0, T _i1 ) : i0( _i0 ), i1( _i1 ) {}
};

template< class T, class execution_space >
T ModAtomic( T i0, T i1 ) {
  struct InitFunctor< T, execution_space > f_init( i0 );
  typename InitFunctor< T, execution_space >::type data( "Data" );
  typename InitFunctor< T, execution_space >::h_type h_data( "HData" );

  f_init.data = data;
  Kokkos::parallel_for( 1, f_init );
  execution_space::fence();

  struct ModFunctor< T, execution_space > f( i0, i1 );

  f.data = data;
  Kokkos::parallel_for( 1, f );
  execution_space::fence();

  Kokkos::deep_copy( h_data, data );
  T val = h_data();

  return val;
}

template< class T >
T ModAtomicCheck( T i0, T i1 ) {
  T* data = new T[1];
  data[0] = 0;

  *data = i0 % i1;

  T val = *data;
  delete [] data;

  return val;
}

template< class T, class DeviceType >
bool ModAtomicTest( T i0, T i1 )
{
  T res       = ModAtomic< T, DeviceType >( i0, i1 );
  T resSerial = ModAtomicCheck< T >( i0, i1 );

  bool passed = true;

  if ( resSerial != res ) {
    passed = false;

    std::cout << "Loop<"
              << typeid( T ).name()
              << ">( test = ModAtomicTest"
              << " FAILED : "
              << resSerial << " != " << res
              << std::endl;
  }

  return passed;
}

//---------------------------------------------------
//--------------atomic_fetch_and---------------------
//---------------------------------------------------

template< class T, class DEVICE_TYPE >
struct AndFunctor {
  typedef DEVICE_TYPE execution_space;
  typedef Kokkos::View< T, execution_space > type;

  type data;
  T i0;
  T i1;

  KOKKOS_INLINE_FUNCTION
  void operator()( int ) const {
    Kokkos::atomic_fetch_and( &data(), (T) i1 );
  }

  AndFunctor( T _i0, T _i1 ) : i0( _i0 ), i1( _i1 ) {}
};

template< class T, class execution_space >
T AndAtomic( T i0, T i1 ) {
  struct InitFunctor< T, execution_space > f_init( i0 );
  typename InitFunctor< T, execution_space >::type data( "Data" );
  typename InitFunctor< T, execution_space >::h_type h_data( "HData" );

  f_init.data = data;
  Kokkos::parallel_for( 1, f_init );
  execution_space::fence();

  struct AndFunctor< T, execution_space > f( i0, i1 );

  f.data = data;
  Kokkos::parallel_for( 1, f );
  execution_space::fence();

  Kokkos::deep_copy( h_data, data );
  T val = h_data();

  return val;
}

template< class T >
T AndAtomicCheck( T i0, T i1 ) {
  T* data = new T[1];
  data[0] = 0;

  *data = i0 & i1;

  T val = *data;
  delete [] data;

  return val;
}

template< class T, class DeviceType >
bool AndAtomicTest( T i0, T i1 )
{
  T res       = AndAtomic< T, DeviceType >( i0, i1 );
  T resSerial = AndAtomicCheck< T >( i0, i1 );

  bool passed = true;

  if ( resSerial != res ) {
    passed = false;

    std::cout << "Loop<"
              << typeid( T ).name()
              << ">( test = AndAtomicTest"
              << " FAILED : "
              << resSerial << " != " << res
              << std::endl;
  }

  return passed;
}

//---------------------------------------------------
//--------------atomic_fetch_or----------------------
//---------------------------------------------------

template< class T, class DEVICE_TYPE >
struct OrFunctor {
  typedef DEVICE_TYPE execution_space;
  typedef Kokkos::View< T, execution_space > type;

  type data;
  T i0;
  T i1;

  KOKKOS_INLINE_FUNCTION
  void operator()( int ) const {
    Kokkos::atomic_fetch_or( &data(), (T) i1 );
  }

  OrFunctor( T _i0, T _i1 ) : i0( _i0 ), i1( _i1 ) {}
};

template< class T, class execution_space >
T OrAtomic( T i0, T i1 ) {
  struct InitFunctor< T, execution_space > f_init( i0 );
  typename InitFunctor< T, execution_space >::type data( "Data" );
  typename InitFunctor< T, execution_space >::h_type h_data( "HData" );

  f_init.data = data;
  Kokkos::parallel_for( 1, f_init );
  execution_space::fence();

  struct OrFunctor< T, execution_space > f( i0, i1 );

  f.data = data;
  Kokkos::parallel_for( 1, f );
  execution_space::fence();

  Kokkos::deep_copy( h_data, data );
  T val = h_data();

  return val;
}

template< class T >
T OrAtomicCheck( T i0, T i1 ) {
  T* data = new T[1];
  data[0] = 0;

  *data = i0 | i1;

  T val = *data;
  delete [] data;

  return val;
}

template< class T, class DeviceType >
bool OrAtomicTest( T i0, T i1 )
{
  T res       = OrAtomic< T, DeviceType >( i0, i1 );
  T resSerial = OrAtomicCheck< T >( i0, i1 );

  bool passed = true;

  if ( resSerial != res ) {
    passed = false;

    std::cout << "Loop<"
              << typeid( T ).name()
              << ">( test = OrAtomicTest"
              << " FAILED : "
              << resSerial << " != " << res
              << std::endl;
  }

  return passed;
}

//---------------------------------------------------
//--------------atomic_fetch_xor---------------------
//---------------------------------------------------

template< class T, class DEVICE_TYPE >
struct XorFunctor {
  typedef DEVICE_TYPE execution_space;
  typedef Kokkos::View< T, execution_space > type;

  type data;
  T i0;
  T i1;

  KOKKOS_INLINE_FUNCTION
  void operator()( int ) const {
    Kokkos::atomic_fetch_xor( &data(), (T) i1 );
  }

  XorFunctor( T _i0, T _i1 ) : i0( _i0 ), i1( _i1 ) {}
};

template< class T, class execution_space >
T XorAtomic( T i0, T i1 ) {
  struct InitFunctor< T, execution_space > f_init( i0 );
  typename InitFunctor< T, execution_space >::type data( "Data" );
  typename InitFunctor< T, execution_space >::h_type h_data( "HData" );

  f_init.data = data;
  Kokkos::parallel_for( 1, f_init );
  execution_space::fence();

  struct XorFunctor< T, execution_space > f( i0, i1 );

  f.data = data;
  Kokkos::parallel_for( 1, f );
  execution_space::fence();

  Kokkos::deep_copy( h_data, data );
  T val = h_data();

  return val;
}

template< class T >
T XorAtomicCheck( T i0, T i1 ) {
  T* data = new T[1];
  data[0] = 0;

  *data = i0 ^ i1;

  T val = *data;
  delete [] data;

  return val;
}

template< class T, class DeviceType >
bool XorAtomicTest( T i0, T i1 )
{
  T res       = XorAtomic< T, DeviceType >( i0, i1 );
  T resSerial = XorAtomicCheck< T >( i0, i1 );

  bool passed = true;

  if ( resSerial != res ) {
    passed = false;

    std::cout << "Loop<"
              << typeid( T ).name()
              << ">( test = XorAtomicTest"
              << " FAILED : "
              << resSerial << " != " << res
              << std::endl;
  }

  return passed;
}

//---------------------------------------------------
//--------------atomic_fetch_lshift---------------------
//---------------------------------------------------

template< class T, class DEVICE_TYPE >
struct LShiftFunctor {
  typedef DEVICE_TYPE execution_space;
  typedef Kokkos::View< T, execution_space > type;

  type data;
  T i0;
  T i1;

  KOKKOS_INLINE_FUNCTION
  void operator()( int ) const {
    Kokkos::atomic_fetch_lshift( &data(), (T) i1 );
  }

  LShiftFunctor( T _i0, T _i1 ) : i0( _i0 ), i1( _i1 ) {}
};

template< class T, class execution_space >
T LShiftAtomic( T i0, T i1 ) {
  struct InitFunctor< T, execution_space > f_init( i0 );
  typename InitFunctor< T, execution_space >::type data( "Data" );
  typename InitFunctor< T, execution_space >::h_type h_data( "HData" );

  f_init.data = data;
  Kokkos::parallel_for( 1, f_init );
  execution_space::fence();

  struct LShiftFunctor< T, execution_space > f( i0, i1 );

  f.data = data;
  Kokkos::parallel_for( 1, f );
  execution_space::fence();

  Kokkos::deep_copy( h_data, data );
  T val = h_data();

  return val;
}

template< class T >
T LShiftAtomicCheck( T i0, T i1 ) {
  T* data = new T[1];
  data[0] = 0;

  *data = i0 << i1;

  T val = *data;
  delete [] data;

  return val;
}

template< class T, class DeviceType >
bool LShiftAtomicTest( T i0, T i1 )
{
  T res       = LShiftAtomic< T, DeviceType >( i0, i1 );
  T resSerial = LShiftAtomicCheck< T >( i0, i1 );

  bool passed = true;

  if ( resSerial != res ) {
    passed = false;

    std::cout << "Loop<"
              << typeid( T ).name()
              << ">( test = LShiftAtomicTest"
              << " FAILED : "
              << resSerial << " != " << res
              << std::endl;
  }

  return passed;
}

//---------------------------------------------------
//--------------atomic_fetch_rshift---------------------
//---------------------------------------------------

template< class T, class DEVICE_TYPE >
struct RShiftFunctor {
  typedef DEVICE_TYPE execution_space;
  typedef Kokkos::View< T, execution_space > type;

  type data;
  T i0;
  T i1;

  KOKKOS_INLINE_FUNCTION
  void operator()( int ) const {
    Kokkos::atomic_fetch_rshift( &data(), (T) i1 );
  }

  RShiftFunctor( T _i0, T _i1 ) : i0( _i0 ), i1( _i1 ) {}
};

template< class T, class execution_space >
T RShiftAtomic( T i0, T i1 ) {
  struct InitFunctor< T, execution_space > f_init( i0 );
  typename InitFunctor< T, execution_space >::type data( "Data" );
  typename InitFunctor< T, execution_space >::h_type h_data( "HData" );

  f_init.data = data;
  Kokkos::parallel_for( 1, f_init );
  execution_space::fence();

  struct RShiftFunctor< T, execution_space > f( i0, i1 );

  f.data = data;
  Kokkos::parallel_for( 1, f );
  execution_space::fence();

  Kokkos::deep_copy( h_data, data );
  T val = h_data();

  return val;
}

template< class T >
T RShiftAtomicCheck( T i0, T i1 ) {
  T* data = new T[1];
  data[0] = 0;

  *data = i0 >> i1;

  T val = *data;
  delete [] data;

  return val;
}

template< class T, class DeviceType >
bool RShiftAtomicTest( T i0, T i1 )
{
  T res       = RShiftAtomic< T, DeviceType >( i0, i1 );
  T resSerial = RShiftAtomicCheck< T >( i0, i1 );

  bool passed = true;

  if ( resSerial != res ) {
    passed = false;

    std::cout << "Loop<"
              << typeid( T ).name()
              << ">( test = RShiftAtomicTest"
              << " FAILED : "
              << resSerial << " != " << res
              << std::endl;
  }

  return passed;
}

//---------------------------------------------------
//--------------atomic_test_control------------------
//---------------------------------------------------

template< class T, class DeviceType >
bool AtomicOperationsTestIntegralType( int i0, int i1, int test )
{
  switch ( test ) {
    case 1: return MaxAtomicTest< T, DeviceType >( (T) i0, (T) i1 );
    case 2: return MinAtomicTest< T, DeviceType >( (T) i0, (T) i1 );
    case 3: return MulAtomicTest< T, DeviceType >( (T) i0, (T) i1 );
    case 4: return DivAtomicTest< T, DeviceType >( (T) i0, (T) i1 );
    case 5: return ModAtomicTest< T, DeviceType >( (T) i0, (T) i1 );
    case 6: return AndAtomicTest< T, DeviceType >( (T) i0, (T) i1 );
    case 7: return OrAtomicTest< T, DeviceType >( (T) i0, (T) i1 );
    case 8: return XorAtomicTest< T, DeviceType >( (T) i0, (T) i1 );
    case 9: return LShiftAtomicTest< T, DeviceType >( (T) i0, (T) i1 );
    case 10: return RShiftAtomicTest< T, DeviceType >( (T) i0, (T) i1 );
    case 11: return IncAtomicTest< T, DeviceType >( (T) i0 );
    case 12: return DecAtomicTest< T, DeviceType >( (T) i0 );
  }

  return 0;
}

template< class T, class DeviceType >
bool AtomicOperationsTestNonIntegralType( int i0, int i1, int test )
{
  switch ( test ) {
    case 1: return MaxAtomicTest< T, DeviceType >( (T) i0, (T) i1 );
    case 2: return MinAtomicTest< T, DeviceType >( (T) i0, (T) i1 );
    case 3: return MulAtomicTest< T, DeviceType >( (T) i0, (T) i1 );
    case 4: return DivAtomicTest< T, DeviceType >( (T) i0, (T) i1 );
  }

  return 0;
}

}

namespace Test {

TEST_F( TEST_CATEGORY , atomic_operations )
{
  const int start = 1; // Avoid zero for division.
  const int end = 11;
  for ( int i = start; i < end; ++i )
  {
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< int, TEST_EXECSPACE >( start, end - i, 1 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< int, TEST_EXECSPACE >( start, end - i, 2 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< int, TEST_EXECSPACE >( start, end - i, 3 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< int, TEST_EXECSPACE >( start, end - i, 4 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< int, TEST_EXECSPACE >( start, end - i, 5 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< int, TEST_EXECSPACE >( start, end - i, 6 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< int, TEST_EXECSPACE >( start, end - i, 7 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< int, TEST_EXECSPACE >( start, end - i, 8 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< int, TEST_EXECSPACE >( start, end - i, 9 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< int, TEST_EXECSPACE >( start, end - i, 11 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< int, TEST_EXECSPACE >( start, end - i, 12 ) ) );

    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< unsigned int, TEST_EXECSPACE >( start, end - i, 1 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< unsigned int, TEST_EXECSPACE >( start, end - i, 2 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< unsigned int, TEST_EXECSPACE >( start, end - i, 3 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< unsigned int, TEST_EXECSPACE >( start, end - i, 4 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< unsigned int, TEST_EXECSPACE >( start, end - i, 5 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< unsigned int, TEST_EXECSPACE >( start, end - i, 6 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< unsigned int, TEST_EXECSPACE >( start, end - i, 7 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< unsigned int, TEST_EXECSPACE >( start, end - i, 8 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< unsigned int, TEST_EXECSPACE >( start, end - i, 9 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< unsigned int, TEST_EXECSPACE >( start, end - i, 11 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< unsigned int, TEST_EXECSPACE >( start, end - i, 12 ) ) );

    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< long int, TEST_EXECSPACE >( start, end - i, 1 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< long int, TEST_EXECSPACE >( start, end - i, 2 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< long int, TEST_EXECSPACE >( start, end - i, 3 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< long int, TEST_EXECSPACE >( start, end - i, 4 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< long int, TEST_EXECSPACE >( start, end - i, 5 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< long int, TEST_EXECSPACE >( start, end - i, 6 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< long int, TEST_EXECSPACE >( start, end - i, 7 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< long int, TEST_EXECSPACE >( start, end - i, 8 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< long int, TEST_EXECSPACE >( start, end - i, 9 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< long int, TEST_EXECSPACE >( start, end - i, 11 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< long int, TEST_EXECSPACE >( start, end - i, 12 ) ) );

    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< unsigned long int, TEST_EXECSPACE >( start, end - i, 1 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< unsigned long int, TEST_EXECSPACE >( start, end - i, 2 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< unsigned long int, TEST_EXECSPACE >( start, end - i, 3 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< unsigned long int, TEST_EXECSPACE >( start, end - i, 4 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< unsigned long int, TEST_EXECSPACE >( start, end - i, 5 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< unsigned long int, TEST_EXECSPACE >( start, end - i, 6 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< unsigned long int, TEST_EXECSPACE >( start, end - i, 7 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< unsigned long int, TEST_EXECSPACE >( start, end - i, 8 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< unsigned long int, TEST_EXECSPACE >( start, end - i, 9 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< unsigned long int, TEST_EXECSPACE >( start, end - i, 11 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< unsigned long int, TEST_EXECSPACE >( start, end - i, 12 ) ) );

    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< long long int, TEST_EXECSPACE >( start, end - i, 1 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< long long int, TEST_EXECSPACE >( start, end - i, 2 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< long long int, TEST_EXECSPACE >( start, end - i, 3 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< long long int, TEST_EXECSPACE >( start, end - i, 4 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< long long int, TEST_EXECSPACE >( start, end - i, 5 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< long long int, TEST_EXECSPACE >( start, end - i, 6 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< long long int, TEST_EXECSPACE >( start, end - i, 7 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< long long int, TEST_EXECSPACE >( start, end - i, 8 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< long long int, TEST_EXECSPACE >( start, end - i, 9 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< long long int, TEST_EXECSPACE >( start, end - i, 11 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< long long int, TEST_EXECSPACE >( start, end - i, 12 ) ) );

    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestNonIntegralType< double, TEST_EXECSPACE >( start, end - i, 1 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestNonIntegralType< double, TEST_EXECSPACE >( start, end - i, 2 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestNonIntegralType< double, TEST_EXECSPACE >( start, end - i, 3 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestNonIntegralType< double, TEST_EXECSPACE >( start, end - i, 4 ) ) );

    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestNonIntegralType< float, TEST_EXECSPACE >( start, end - i, 1 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestNonIntegralType< float, TEST_EXECSPACE >( start, end - i, 2 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestNonIntegralType< float, TEST_EXECSPACE >( start, end - i, 3 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestNonIntegralType< float, TEST_EXECSPACE >( start, end - i, 4 ) ) );
  }
}
} // namespace Test
