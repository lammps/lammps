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

namespace Kokkos {
namespace Experimental {
namespace Impl {

template< class MemorySpace = void , class DestroyFunctor = void >
class SharedAllocationRecord ;

class SharedAllocationHeader {
private:

  typedef SharedAllocationRecord<void,void>  Record ;

  static constexpr unsigned maximum_label_length = ( 1u << 7 /* 128 */ ) - sizeof(Record*);

  template< class , class > friend class SharedAllocationRecord ;

  Record * m_record ;
  char     m_label[ maximum_label_length ];

public:

  /* Given user memory get pointer to the header */
  KOKKOS_INLINE_FUNCTION static
  const SharedAllocationHeader * get_header( void * alloc_ptr )
    { return reinterpret_cast<SharedAllocationHeader*>( reinterpret_cast<char*>(alloc_ptr) - sizeof(SharedAllocationHeader) ); }
};

template<>
class SharedAllocationRecord< void , void > {
protected:

  static_assert( sizeof(SharedAllocationHeader) == ( 1u << 7 /* 128 */ ) , "sizeof(SharedAllocationHeader) != 128" );

  template< class , class > friend class SharedAllocationRecord ;

  typedef void (* function_type )( SharedAllocationRecord<void,void> * );

  SharedAllocationHeader * const m_alloc_ptr ;
  size_t                   const m_alloc_size ;
  function_type            const m_dealloc ;
  SharedAllocationRecord * const m_root ;
  SharedAllocationRecord *       m_prev ;
  SharedAllocationRecord *       m_next ;
  int                            m_count ;

  SharedAllocationRecord( const SharedAllocationRecord & ) = delete ;
  SharedAllocationRecord & operator = ( const SharedAllocationRecord & ) = delete ;

  /**\brief  Construct and insert into 'arg_root' tracking set.
   *         use_count is zero.
   */
  SharedAllocationRecord( SharedAllocationRecord * arg_root
                        , SharedAllocationHeader * arg_alloc_ptr
                        , size_t                   arg_alloc_size
                        , function_type            arg_dealloc
                        );

public:

  ~SharedAllocationRecord() = default ;

  constexpr SharedAllocationRecord()
    : m_alloc_ptr( 0 )
    , m_alloc_size( 0 )
    , m_dealloc( 0 )
    , m_root( this )
    , m_prev( this )
    , m_next( this )
    , m_count( 0 )
    {}

  static constexpr unsigned maximum_label_length = SharedAllocationHeader::maximum_label_length ;

  KOKKOS_INLINE_FUNCTION
  const SharedAllocationHeader * head() const { return m_alloc_ptr ; }

  /* User's memory begins at the end of the header */
  KOKKOS_INLINE_FUNCTION
  void * data() const { return reinterpret_cast<void*>( m_alloc_ptr + 1 ); }

  /* User's memory begins at the end of the header */
  constexpr size_t size() const { return m_alloc_size - sizeof(SharedAllocationHeader) ; }

  /* Cannot be 'constexpr' because 'm_count' is volatile */
  int use_count() const { return m_count ; }

  /* Increment use count */
  static void increment( SharedAllocationRecord * );

  /* Decrement use count. If 1->0 then remove from the tracking list and invoke m_dealloc */
  static SharedAllocationRecord * decrement( SharedAllocationRecord * );

  /* Given a root record and data pointer find the record */
  static SharedAllocationRecord * find( SharedAllocationRecord * const , void * const );

  /*  Sanity check for the whole set of records to which the input record belongs.
   *  Locks the set's insert/erase operations until the sanity check is complete.
   */
  static bool is_sane( SharedAllocationRecord * );

  /*  Print host-accessible records */
  static void print_host_accessible_records( std::ostream &
                                           , const char * const space_name
                                           , const SharedAllocationRecord * const root
                                           , const bool detail );
};

/*
 *  Memory space specialization of SharedAllocationRecord< Space , void > requires :
 *
 *  SharedAllocationRecord< Space , void > : public SharedAllocationRecord< void , void >
 *  {
 *    // delete allocated user memory via static_cast to this type.
 *    static void deallocate( const SharedAllocationRecord<void,void> * );
 *    Space m_space ;
 *  }
 */

template< class MemorySpace , class DestroyFunctor >
class SharedAllocationRecord : public SharedAllocationRecord< MemorySpace , void >
{
private:

  static void deallocate( SharedAllocationRecord<void,void> * record_ptr )
    { delete static_cast<SharedAllocationRecord<MemorySpace,DestroyFunctor>*>(record_ptr); }

  SharedAllocationRecord( const MemorySpace & arg_space
                        , const std::string & arg_label
                        , const size_t        arg_alloc
                        )
    /*  Allocate user memory as [ SharedAllocationHeader , user_memory ] */
    : SharedAllocationRecord< MemorySpace , void >( arg_space , arg_label , arg_alloc , & deallocate )
    , m_destroy()
    {}

  ~SharedAllocationRecord() { m_destroy.destroy_shared_allocation(); }

public:

  DestroyFunctor  m_destroy ;

  // Allocate with a zero use count.  Incrementing the use count from zero to one
  // inserts the record into the tracking list.  Decrementing the count from one to zero
  // removes from the trakcing list and deallocates.
  KOKKOS_INLINE_FUNCTION static
  SharedAllocationRecord * allocate( const MemorySpace & arg_space 
                                   , const std::string & arg_label
                                   , const size_t        arg_alloc
                                   )
    {
#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
      return new SharedAllocationRecord( arg_space , arg_label , arg_alloc );
#else
      return (SharedAllocationRecord *) 0 ;
#endif
    }
};

union SharedAllocationTracker {
private:

  typedef SharedAllocationRecord<void,void>  Record ;

  enum : unsigned long {
    DO_NOT_DEREF_FLAG = 0x01ul
  };

  // The allocation record resides in Host memory space
  Record * m_record ;
  unsigned long m_record_bits;

  KOKKOS_INLINE_FUNCTION
  static Record * disable( Record * rec )
    { return reinterpret_cast<Record*>( reinterpret_cast<unsigned long>( rec ) & DO_NOT_DEREF_FLAG ); }

  KOKKOS_INLINE_FUNCTION
  void increment() const
    {
#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
      if ( ! ( m_record_bits & DO_NOT_DEREF_FLAG ) ) Record::increment( m_record );
#endif
    }

  KOKKOS_INLINE_FUNCTION
  void decrement() const
    {
#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
       if ( ! ( m_record_bits & DO_NOT_DEREF_FLAG ) ) Record::decrement( m_record );
#endif
    }

public:

  KOKKOS_INLINE_FUNCTION
  constexpr SharedAllocationTracker() : m_record_bits( DO_NOT_DEREF_FLAG ) {}

  template< class MemorySpace >
  constexpr
  SharedAllocationRecord< MemorySpace , void > & get_record() const
    { return * static_cast< SharedAllocationRecord< MemorySpace , void > * >( m_record ); }

  template< class MemorySpace >
  std::string get_label() const
    { return static_cast< SharedAllocationRecord< MemorySpace , void > * >( m_record )->get_label(); }

  KOKKOS_INLINE_FUNCTION
  SharedAllocationTracker( Record * arg_record )
    : m_record( arg_record ) { increment(); }

  KOKKOS_INLINE_FUNCTION
  ~SharedAllocationTracker() { decrement(); }

  KOKKOS_INLINE_FUNCTION
  SharedAllocationTracker( const SharedAllocationTracker & rhs )
    : m_record( rhs.m_record ) { increment(); }

  KOKKOS_INLINE_FUNCTION
  SharedAllocationTracker( SharedAllocationTracker && rhs )
    : m_record( rhs.m_record ) { rhs.m_record_bits = DO_NOT_DEREF_FLAG ; }

  KOKKOS_INLINE_FUNCTION
  SharedAllocationTracker & operator = ( const SharedAllocationTracker & rhs )
    {
      decrement();
      m_record = rhs.m_record ;
      increment();
      return *this ;
    }

  KOKKOS_INLINE_FUNCTION
  SharedAllocationTracker & operator = ( SharedAllocationTracker && rhs )
    {
      m_record = rhs.m_record ;
      rhs.m_record_bits = DO_NOT_DEREF_FLAG ;
      return *this ;
    }
};


} /* namespace Impl */
} /* namespace Experimental */
} /* namespace Kokkos */


