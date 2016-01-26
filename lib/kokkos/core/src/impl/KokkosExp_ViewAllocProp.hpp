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

#ifndef KOKKOS_EXPERIMENTAL_IMPL_VIEW_ALLOC_PROP_HPP
#define KOKKOS_EXPERIMENTAL_IMPL_VIEW_ALLOC_PROP_HPP

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#if defined( KOKKOS_USING_EXPERIMENTAL_VIEW )

namespace Kokkos {

/* For backward compatibility */

struct ViewAllocateWithoutInitializing {

  const std::string label ;

  ViewAllocateWithoutInitializing() : label() {}
  ViewAllocateWithoutInitializing( const std::string & arg_label ) : label( arg_label ) {}
  ViewAllocateWithoutInitializing( const char * const  arg_label ) : label( arg_label ) {}
};

} /* namespace Kokkos */

#endif

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Experimental {
namespace Impl {

struct WithoutInitializing_t {};
struct AllowPadding_t {};

template< class ... Parameters >
struct ViewAllocProp ;

template<>
struct ViewAllocProp<> {

  struct NullSpace {};

  typedef std::false_type  allow_padding_t ;
  typedef std::true_type   initialize_t ;
  typedef NullSpace        memory_space ;
  typedef NullSpace        execution_space ;

  const std::string      label ;
  const memory_space     memory ;
  const execution_space  execution ;
  const allow_padding_t  allow_padding ;
  const initialize_t     initialize ;

  ViewAllocProp()
    : label()
    , memory()
    , execution()
    , allow_padding()
    , initialize()
    {}

  ViewAllocProp( const std::string & arg_label )
    : label( arg_label )
    , memory()
    , execution()
    , allow_padding()
    , initialize()
    {}
};

template< class ... Parameters >
struct ViewAllocProp< const char * , Parameters ... >
{
  typedef ViewAllocProp< Parameters ... >  base_prop_type ;

  typedef typename base_prop_type::allow_padding_t  allow_padding_t ;
  typedef typename base_prop_type::initialize_t     initialize_t ;
  typedef typename base_prop_type::memory_space     memory_space ;
  typedef typename base_prop_type::execution_space  execution_space ;

  const std::string      label ;
  const memory_space     memory ;
  const execution_space  execution ;
  const allow_padding_t  allow_padding ;
  const initialize_t     initialize ;

  ViewAllocProp( const char * const arg_label , Parameters ... arg_param )
    : label( arg_label )
    , memory( base_prop_type( arg_param ... ).memory )
    , execution( base_prop_type( arg_param ... ).execution )
    , allow_padding()
    , initialize()
    {}
};

template< class ... Parameters >
struct ViewAllocProp< std::string , Parameters ... >
{
  typedef ViewAllocProp< Parameters ... >  base_prop_type ;

  typedef typename base_prop_type::allow_padding_t  allow_padding_t ;
  typedef typename base_prop_type::initialize_t     initialize_t ;
  typedef typename base_prop_type::memory_space     memory_space ;
  typedef typename base_prop_type::execution_space  execution_space ;

  const std::string      label ;
  const memory_space     memory ;
  const execution_space  execution ;
  const allow_padding_t  allow_padding ;
  const initialize_t     initialize ;

  ViewAllocProp( const std::string & arg_label , Parameters ... arg_param )
    : label( arg_label )
    , memory( base_prop_type( arg_param ... ).memory )
    , execution( base_prop_type( arg_param ... ).execution )
    , allow_padding()
    , initialize()
    {}
};

template< class ... Parameters >
struct ViewAllocProp< WithoutInitializing_t , Parameters ... >
{
  typedef ViewAllocProp< Parameters ... >  base_prop_type ;

  typedef typename base_prop_type::allow_padding_t  allow_padding_t ;
  typedef std::false_type                           initialize_t ;
  typedef typename base_prop_type::memory_space     memory_space ;
  typedef typename base_prop_type::execution_space  execution_space ;

  const std::string      label ;
  const memory_space     memory ;
  const execution_space  execution ;
  const allow_padding_t  allow_padding ;
  const initialize_t     initialize ;

  ViewAllocProp( const WithoutInitializing_t & , Parameters ... arg_param )
    : label( base_prop_type( arg_param ... ).label )
    , memory( base_prop_type( arg_param ... ).memory )
    , execution( base_prop_type( arg_param ... ).execution )
    , allow_padding()
    , initialize()
    {}
};

template< class ... Parameters >
struct ViewAllocProp< AllowPadding_t , Parameters ... >
{
  typedef ViewAllocProp< Parameters ... >  base_prop_type ;

  typedef std::true_type                            allow_padding_t ;
  typedef typename base_prop_type::initialize_t     initialize_t ;
  typedef typename base_prop_type::memory_space     memory_space ;
  typedef typename base_prop_type::execution_space  execution_space ;

  const std::string label ;
  const memory_space     memory ;
  const execution_space  execution ;
  const allow_padding_t  allow_padding ;
  const initialize_t     initialize ;

  ViewAllocProp( const AllowPadding_t & , Parameters ... arg_param )
    : label( base_prop_type( arg_param ... ).label )
    , memory( base_prop_type( arg_param ... ).memory )
    , execution( base_prop_type( arg_param ... ).execution )
    , allow_padding()
    , initialize()
    {}
};

template< class Space , class ... Parameters >
struct ViewAllocProp< Space , Parameters ... >
{
  enum { is_exec = Kokkos::Impl::is_execution_space< Space >::value };
  enum { is_mem  = Kokkos::Impl::is_memory_space< Space >::value };

  static_assert( is_exec || is_mem , "View allocation given unknown parameter" );

  typedef ViewAllocProp< Parameters ... >  base_prop_type ;

  typedef typename base_prop_type::allow_padding_t  allow_padding_t ;
  typedef typename base_prop_type::initialize_t     initialize_t ;
  typedef typename std::conditional< is_mem  , Space , typename base_prop_type::memory_space >::type     memory_space ;
  typedef typename std::conditional< is_exec , Space , typename base_prop_type::execution_space >::type  execution_space ;

  const std::string      label ;
  const memory_space     memory ;
  const execution_space  execution ;
  const allow_padding_t  allow_padding ;
  const initialize_t     initialize ;

  // Templated so that 'base_prop_type( args ... ).execution'
  // is not used unless arg_space == memory_space.
  template< class ... Args >
  ViewAllocProp( const memory_space & arg_space , Args ... args )
    : label( base_prop_type( args ... ).label )
    , memory( arg_space )
    , execution( base_prop_type( args ... ).execution )
    , allow_padding()
    , initialize()
    {}

  // Templated so that 'base_prop_type( args ... ).memory'
  // is not used unless arg_space == execution_space.
  template< class ... Args >
  ViewAllocProp( const execution_space & arg_space , Args ... args )
    : label( base_prop_type( args ... ).label )
    , memory( base_prop_type( args ... ).memory )
    , execution( arg_space )
    , allow_padding()
    , initialize()
    {}
};

template< class ExecSpace , class MemSpace >
struct ViewAllocProp< Kokkos::Device< ExecSpace , MemSpace > , std::string >
{
  typedef ViewAllocProp<>  base_prop_type ;

  typedef typename base_prop_type::allow_padding_t  allow_padding_t ;
  typedef typename base_prop_type::initialize_t     initialize_t ;
  typedef MemSpace   memory_space ;
  typedef ExecSpace  execution_space ;

  const std::string      label ;
  const memory_space     memory ;
  const execution_space  execution ;
  const allow_padding_t  allow_padding ;
  const initialize_t     initialize ;

  ViewAllocProp( const std::string & arg_label )
    : label( arg_label )
    , memory()
    , execution()
    , allow_padding()
    , initialize()
    {}
};

template< class ExecSpace , class MemSpace , unsigned N >
struct ViewAllocProp< Kokkos::Device< ExecSpace , MemSpace > , char[N] >
{
  typedef ViewAllocProp<>  base_prop_type ;

  typedef typename base_prop_type::allow_padding_t  allow_padding_t ;
  typedef typename base_prop_type::initialize_t     initialize_t ;
  typedef MemSpace   memory_space ;
  typedef ExecSpace  execution_space  ;

  const std::string      label ;
  const memory_space     memory ;
  const execution_space  execution ;
  const allow_padding_t  allow_padding ;
  const initialize_t     initialize ;

  ViewAllocProp( const char * const arg_label )
    : label( arg_label )
    , memory()
    , execution()
    , allow_padding()
    , initialize()
    {}
};


// Deprecate in favor of view_alloc( Kokkos::WithoutInitializing )
template< class ExecSpace , class MemSpace >
struct ViewAllocProp< Kokkos::Device< ExecSpace , MemSpace >
                    , Kokkos::ViewAllocateWithoutInitializing
                    >
{
  typedef ViewAllocProp<>  base_prop_type ;

  typedef typename base_prop_type::allow_padding_t  allow_padding_t ;
  typedef std::false_type                           initialize_t ;
  typedef MemSpace   memory_space ;
  typedef ExecSpace  execution_space  ;

  const std::string      label ;
  const memory_space     memory ;
  const execution_space  execution ;
  const allow_padding_t  allow_padding ;
  const initialize_t     initialize ;

  ViewAllocProp( const Kokkos::ViewAllocateWithoutInitializing & arg )
    : label( arg.label )
    , memory()
    , execution()
    , allow_padding()
    , initialize()
    {}
};

template< class ExecSpace , class MemSpace , class ... Parameters >
struct ViewAllocProp< Kokkos::Device< ExecSpace , MemSpace >
                    , ViewAllocProp< Parameters ... >
                    >
{
  typedef ViewAllocProp< Parameters ... >  base_prop_type ;

  typedef typename base_prop_type::allow_padding_t  allow_padding_t ;
  typedef typename base_prop_type::initialize_t     initialize_t ;
  typedef MemSpace  memory_space ;

  typedef
    typename std::conditional
      < Kokkos::Impl::is_execution_space< typename base_prop_type::execution_space >::value
      , typename base_prop_type::execution_space
      , ExecSpace
      >::type  execution_space ;

  static_assert( std::is_same< typename base_prop_type::memory_space , ViewAllocProp<>::NullSpace >::value ||
                 std::is_same< typename base_prop_type::memory_space , memory_space >::value
               , "View allocation given incompatible memory space" );

  static_assert( Kokkos::Impl::VerifyExecutionCanAccessMemorySpace< typename execution_space::memory_space
                                                                  , memory_space >::value
               , "View allocation given incompatible execution space" );

  const std::string      label ;
  const memory_space     memory ;
  const execution_space  execution ;
  const allow_padding_t  allow_padding ;
  const initialize_t     initialize ;

  // If the input properties have a memory or execution space then copy construct those spaces
  // otherwise default construct those spaces.

  template< class P >
  ViewAllocProp( const P & arg_prop
               , typename std::enable_if
                   < std::is_same< P , base_prop_type >::value &&
                     Kokkos::Impl::is_memory_space< typename P::memory_space >::value &&
                     Kokkos::Impl::is_execution_space< typename P::memory_space >::value
                   >::type * = 0 )
    : label( arg_prop.label )
    , memory( arg_prop.memory )
    , execution( arg_prop.execution )
    , allow_padding()
    , initialize()
    {}

  template< class P >
  ViewAllocProp( const P & arg_prop
               , typename std::enable_if
                   < std::is_same< P , base_prop_type >::value &&
                     Kokkos::Impl::is_memory_space< typename P::memory_space >::value &&
                     ! Kokkos::Impl::is_execution_space< typename P::execution_space >::value
                   >::type * = 0 )
    : label( arg_prop.label )
    , memory( arg_prop.memory )
    , execution()
    , allow_padding()
    , initialize()
    {}

  template< class P >
  ViewAllocProp( const P & arg_prop
               , typename std::enable_if
                   < std::is_same< P , base_prop_type >::value &&
                     ! Kokkos::Impl::is_memory_space< typename P::memory_space >::value &&
                     Kokkos::Impl::is_execution_space< typename P::execution_space >::value
                   >::type * = 0 )
    : label( arg_prop.label )
    , memory()
    , execution( arg_prop.execution )
    , allow_padding()
    , initialize()
    {}

  template< class P >
  ViewAllocProp( const P & arg_prop
               , typename std::enable_if
                   < std::is_same< P , base_prop_type >::value &&
                     ! Kokkos::Impl::is_memory_space< typename P::memory_space >::value &&
                     ! Kokkos::Impl::is_execution_space< typename P::execution_space >::value
                   >::type * = 0 )
    : label( arg_prop.label )
    , memory()
    , execution()
    , allow_padding()
    , initialize()
    {}
};

} /* namespace Impl */
} /* namespace Experimental */
} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif

