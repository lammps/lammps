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

#include <Kokkos_Core_fwd.hpp>

#if ! defined( KOKKOS_USING_EXPERIMENTAL_VIEW )

#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )

#include <Kokkos_Atomic.hpp>

#include <impl/Kokkos_Singleton.hpp>
#include <impl/Kokkos_AllocationTracker.hpp>
#include <impl/Kokkos_Error.hpp>


#include <string>
#include <vector>
#include <sstream>
#include <algorithm>
#include <utility>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <iomanip>

/* Enable clean up of memory leaks */
#define CLEAN_UP_MEMORY_LEAKS 0

namespace Kokkos { namespace Impl {

namespace {


//-----------------------------------------------------------------------------
// AllocationRecord
//-----------------------------------------------------------------------------
//
// Used to track details about an allocation and provide a ref count
// sizeof(AllocationRecord) == 128
struct AllocationRecord
{
  enum {
     OFFSET = sizeof(AllocatorBase*)          // allocator
            + sizeof(void*)                   // alloc_ptr
            + sizeof(uint64_t)                // alloc_size
            + sizeof(AllocatorAttributeBase*) // attribute
            + sizeof(uint32_t)                // node_index
            + sizeof(uint32_t)                // ref_count
   , LABEL_LENGTH = 128 - OFFSET
  };

  AllocatorBase * const          allocator;
  void * const                   alloc_ptr;
  const uint64_t                 alloc_size;
  AllocatorAttributeBase * const attribute;
  const int32_t                  node_index;
  volatile uint32_t              ref_count;
  const char                     label[LABEL_LENGTH];


  AllocationRecord(  AllocatorBase * const arg_allocator
                   , void *   arg_alloc_ptr
                   , uint64_t arg_alloc_size
                   , int32_t  arg_node_index
                   , const std::string & arg_label
                  )
    : allocator(arg_allocator)
    , alloc_ptr(arg_alloc_ptr)
    , alloc_size(arg_alloc_size)
    , attribute(NULL)
    , node_index(arg_node_index)
    , ref_count(1)
    , label() // zero fill
  {
    const size_t length = static_cast<size_t>(LABEL_LENGTH-1u) < arg_label.size() ? static_cast<size_t>(LABEL_LENGTH-1u) : arg_label.size();
    strncpy( const_cast<char *>(label), arg_label.c_str(), length );
  }

  ~AllocationRecord()
  {
    if (attribute) {
      delete attribute;
    }
  }

  uint32_t increment_ref_count()
  {
    uint32_t old_value = atomic_fetch_add( &ref_count, static_cast<uint32_t>(1) );
    return old_value + 1u;
  }

  uint32_t decrement_ref_count()
  {
    uint32_t old_value = atomic_fetch_sub( &ref_count, static_cast<uint32_t>(1) );
    return old_value - 1u;
  }

  void print( std::ostream & oss ) const
  {
    oss << "{ " << allocator->name()
        << " } : \"" << label
        << "\" ref_count(" << ref_count
        << ") memory[ " << alloc_ptr
        << " + " << alloc_size
        << " ]" ;
  }

  bool set_attribute( AllocatorAttributeBase * attr )
  {
    bool result = false;
    if (attribute == NULL) {
      result = NULL == atomic_compare_exchange(  const_cast<AllocatorAttributeBase **>(&attribute)
                                               , reinterpret_cast<AllocatorAttributeBase *>(NULL)
                                               , attr );
    }

    return result;
  }

  // disallow copy and assignment
  AllocationRecord( const AllocationRecord & );
  AllocationRecord & operator=(const AllocationRecord &);
};

template <int NumBlocks>
struct Bitset
{
  enum { blocks = NumBlocks };
  enum { size = blocks * 64 };
  enum { block_mask = 63u };
  enum { block_shift = 6 };

  // used to find free bits in a bitset
  static int count_trailing_zeros(uint64_t x)
  {
    #if defined( KOKKOS_COMPILER_GNU ) || defined( KOKKOS_COMPILER_CLANG ) || defined( KOKKOS_COMPILER_APPLECC )
      return x ? __builtin_ctzll(x) : 64;
    #elif defined( KOKKOS_COMPILER_INTEL )
      enum { shift = 32 };
      enum { mask = (static_cast<uint64_t>(1) << shift) - 1u };
      return (x & mask) ? _bit_scan_forward(static_cast<int>(x & mask)) :
             (x >> shift) ? shift + _bit_scan_forward(static_cast<int>(x >> shift)) :
             64 ;
    #elif defined( KOKKOS_COMPILER_IBM )
      return x ? __cnttz8(x) : 64;
    #else
      int i = 0;
      for (; ((x & (static_cast<uint64_t>(1) << i)) == 0u) && i < 64; ++i ) {}
      return i;
    #endif
  }

  Bitset()
    : m_bits()
  {
    for (int i=0; i < blocks; ++i) {
      m_bits[i] = 0u;
    }
  }

  bool set( int i )
  {
    const uint64_t bit = static_cast<uint64_t>(1) << ( i & block_mask );
    return !( atomic_fetch_or( m_bits + (i >> block_shift), bit ) & bit );
  }

  bool reset( int i )
  {
    const uint64_t bit = static_cast<uint64_t>(1) << ( i & block_mask );
    return atomic_fetch_and( m_bits + (i >> block_shift), ~bit ) & bit;
  }

  bool test( int i )
  {
    const uint64_t block = m_bits[ i >> block_shift ];
    const uint64_t bit = static_cast<uint64_t>(1) << ( i & block_mask );
    return block & bit;
  }

  int find_first_unset() const
  {
    for (int i=0; i < blocks; ++i) {
      const uint64_t block = m_bits[i];
      int b = count_trailing_zeros( ~block );

      if ( b < 64 ) {
        return (i << block_shift) + b;
      }
    }
    return size;
  }

  volatile uint64_t m_bits[blocks];
};

//-----------------------------------------------------------------------------
// AllocationRecordPool -- singleton class
//
// global_alloc_rec_pool is the ONLY instance of this class
//
//-----------------------------------------------------------------------------
// Record AllocationRecords in a lock-free circular list.
// Each node in the list has a buffer with space for 959 ((15*64)-1) records
// managed by a bitset.  Atomics are used to set and reset bits in the bit set.
// The head of the list is atomically updated to the last node found with
// unused space.
//
// Cost time to create an allocation record: amortized O(1), worst case O(num nodes)
// Cost to destroy an allocation recored: O(1)
//
// Singleton allocations are pushed onto a lock-free stack that is destroyed
// after the circular list of allocation records.
struct AllocationRecordPool
{
  enum { BITSET_BLOCKS = 15 };

  typedef Bitset<BITSET_BLOCKS> bitset_type;

  enum { BUFFER_SIZE = (bitset_type::size - 1) * sizeof(AllocationRecord) };

  struct AllocationNode
  {
    AllocationNode()
      : next()
      , bitset()
      , buffer()
    {
      // set the first bit to used
      bitset.set(0);
    }

    void * get_buffer( int32_t node_index )
    {
      return buffer + (node_index-1) * sizeof(AllocationRecord);
    }

    // return 0 if no space is available in the node
    int32_t get_node_index()
    {
      int32_t node_index = 0;
      do {
        node_index = bitset.find_first_unset();

        // successfully claimed a bit
        if ( node_index != bitset.size && bitset.set(node_index) )
        {
          return node_index;
        }
      } while ( node_index != bitset.size );
      return 0;
    }

    void clear_node_index( int32_t node_index )
    {
      bitset.reset(node_index);
    }

    AllocationNode * next;
    bitset_type      bitset;
    char             buffer[BUFFER_SIZE];
  };

  struct SingletonNode
  {
    void * buffer;
    SingletonNode * next;
    Impl::singleton_destroy_function_type destroy;

    SingletonNode( size_t size, Impl::singleton_create_function_type create_func, Impl::singleton_destroy_function_type destroy_func  )
      : buffer(NULL)
      , next(NULL)
      , destroy(destroy_func)
    {
      if (size) {
        buffer = malloc(size);
        create_func(buffer);
      }
    }

    ~SingletonNode()
    {
      if (buffer) {
        try {
          destroy(buffer);
        } catch(...) {}
        free(buffer);
      }
    }
  };

  AllocationRecordPool()
    : head( new AllocationNode() )
    , singleton_head(NULL)
  {
    // setup ring
    head->next = head;
  }

  ~AllocationRecordPool()
  {
    // delete allocation records
    {
      AllocationNode * start = head;

      AllocationNode * curr = start;

      std::vector< std::string > string_vec;

      do {
        AllocationNode * next = curr->next;

        #if defined( KOKKOS_DEBUG_PRINT_ALLOCATION_BITSET )
        // print node bitset
        for (int i=0; i < bitset_type::blocks; ++i ) {
          std::cout << std::hex << std::showbase << curr->bitset.m_bits[i] << "   ";
        }
        std::cout << std::endl;
        #endif

        // bit zero does not map to an AllocationRecord
        for ( int32_t i=1; i < bitset_type::size; ++i )
        {
          if (curr->bitset.test(i)) {
            AllocationRecord * alloc_rec = reinterpret_cast<AllocationRecord *>( curr->get_buffer(i) );

            std::ostringstream oss;
            alloc_rec->print( oss );
            string_vec.push_back( oss.str() );

#if CLEAN_UP_MEMORY_LEAKS
/* Cleaning up memory leaks prevents memory error detection tools
 * from reporting the original source of allocation, which can
 * impede debugging with such tools.
 */
            try {
              destroy(alloc_rec);
            }
            catch(...) {}
#endif
          }
        }

        curr->next = NULL;

        delete curr;

        curr = next;
      } while ( curr != start );

      //if ( !string_vec.empty() ) {
      //  std::sort( string_vec.begin(), string_vec.end() );
      //
      //  std::ostringstream oss;
      //  oss << "Error: Allocation pool destroyed with the following memory leak(s):\n";
      //  for (size_t i=0; i< string_vec.size(); ++i)
      //  {
      //    oss << "   " << string_vec[i] << std::endl;
      //  }
      //
      //  std::cerr << oss.str() << std::endl;
      //}
    }

    // delete singletons
    {
      SingletonNode * curr = singleton_head;

      while (curr) {
        SingletonNode * next = curr->next;
        delete curr;
        curr = next;
      }
    }
  }

  AllocationRecord * create(  AllocatorBase * arg_allocator
                            , void * arg_alloc_ptr
                            , size_t arg_alloc_size
                            , const std::string & arg_label
                           )
  {
    AllocationNode * start = volatile_load(&head);

    AllocationNode * curr = start;


    int32_t node_index = curr->get_node_index();

    if (node_index == 0) {
      curr = volatile_load(&curr->next);
    }

    while (node_index == 0 && curr != start)
    {
      node_index = curr->get_node_index();
      if (node_index == 0) {
        curr = volatile_load(&curr->next);
      }
    }

    // Need to allocate and insert a new node
    if (node_index == 0 && curr == start)
    {
      AllocationNode * new_node = new AllocationNode();

      node_index = new_node->get_node_index();

      AllocationNode * next = NULL;
      do {
        next = volatile_load(&curr->next);
        new_node->next = next;
        memory_fence();
      } while ( next != atomic_compare_exchange( &(curr->next), next, new_node ) );

      curr = new_node;
    }

    void * buffer = curr->get_buffer(node_index);

    // try to set head to curr
    if ( start != curr )
    {
      atomic_compare_exchange( & head, start, curr );
    }

    return new (buffer) AllocationRecord(  arg_allocator
                                         , arg_alloc_ptr
                                         , arg_alloc_size
                                         , node_index
                                         , arg_label
                                        );
  }

  void destroy( AllocationRecord * alloc_rec )
  {
    if (alloc_rec) {
      const int32_t node_index = alloc_rec->node_index;
      AllocationNode * node = get_node( alloc_rec );

      // deallocate memory
      alloc_rec->allocator->deallocate( alloc_rec->alloc_ptr, alloc_rec->alloc_size );

      // call destructor
      alloc_rec->~AllocationRecord();

      // wait for writes to complete
      memory_fence();

      // clear node index
      node->clear_node_index( node_index );
    }
  }

  void * create_singleton( size_t size, Impl::singleton_create_function_type create_func, Impl::singleton_destroy_function_type destroy_func )
  {
    SingletonNode * node = new SingletonNode( size, create_func, destroy_func );
    SingletonNode * next;

    // insert new node at the head of the list
    do {
      next = volatile_load(&singleton_head);
      node->next = next;
    } while ( next != atomic_compare_exchange( &singleton_head, next, node ) );

    return node->buffer;
  }

  void print_memory( std::ostream & out ) const
  {
    AllocationNode * start = head;

    AllocationNode * curr = start;

    std::vector< std::string > string_vec;

    do {
      AllocationNode * next = curr->next;

      // bit zero does not map to an AllocationRecord
      for ( int32_t i=1; i < bitset_type::size; ++i )
      {
        if (curr->bitset.test(i)) {
          AllocationRecord * alloc_rec = reinterpret_cast<AllocationRecord *>( curr->get_buffer(i) );

          std::ostringstream oss;
          alloc_rec->print( oss );
          string_vec.push_back( oss.str() );
        }
      }
      curr = next;
    } while ( curr != start );

    if ( !string_vec.empty() ) {
      std::sort( string_vec.begin(), string_vec.end() );

      std::ostringstream oss;
      oss << "Tracked Memory:" << std::endl;
      for (size_t i=0; i< string_vec.size(); ++i)
      {
        oss << "   " << string_vec[i] << std::endl;
      }
      out << oss.str() << std::endl;
    }
    else {
      out << "No Tracked Memory" << std::endl;
    }
  }

  // find an AllocationRecord such that
  // alloc_ptr <= ptr < alloc_ptr + alloc_size
  // otherwise return NULL
  AllocationRecord * find( void const * ptr, AllocatorBase const * allocator ) const
  {
    AllocationNode * start = head;

    AllocationNode * curr = start;

    char const * const char_ptr = reinterpret_cast<const char *>(ptr);

    do {
      AllocationNode * next = curr->next;

      // bit zero does not map to an AllocationRecord
      for ( int32_t i=1; i < bitset_type::size; ++i )
      {
        if (curr->bitset.test(i)) {
          AllocationRecord * alloc_rec = reinterpret_cast<AllocationRecord *>( curr->get_buffer(i) );

          char const * const alloc_ptr = reinterpret_cast<char const *>(alloc_rec->alloc_ptr);

          if (   (allocator == alloc_rec->allocator)
              && (alloc_ptr <= char_ptr)
              && (char_ptr < (alloc_ptr + alloc_rec->alloc_size)) )
          {
            return alloc_rec;
          }
        }
      }
      curr = next;
    } while ( curr != start );

    return NULL;
  }

private:

  AllocationNode * get_node( AllocationRecord * alloc_rec )
  {
    return reinterpret_cast<AllocationNode *>( alloc_rec - alloc_rec->node_index);
  }

  AllocationNode * head;
  SingletonNode * singleton_head;
};

// create the global pool for allocation records
AllocationRecordPool global_alloc_rec_pool;



// convert a uintptr_t to an AllocationRecord pointer
inline
AllocationRecord * to_alloc_rec( uintptr_t alloc_rec )
{
  return reinterpret_cast<AllocationRecord *>( alloc_rec & ~static_cast<uintptr_t>(1) );
}

} // unnamed namespace

//-----------------------------------------------------------------------------
// Allocation Tracker methods
//-----------------------------------------------------------------------------

// Create a reference counted AllocationTracker
void AllocationTracker::initalize(  AllocatorBase * arg_allocator
                                  , void * arg_alloc_ptr
                                  , size_t arg_alloc_size
                                  , const std::string & arg_label
                                 )
{
  if ( arg_allocator && arg_alloc_ptr && arg_alloc_size) {
    // create record
    AllocationRecord * alloc_rec = global_alloc_rec_pool.create(  arg_allocator
                                                                , arg_alloc_ptr
                                                                , arg_alloc_size
                                                                , arg_label
                                                               );

    m_alloc_rec = reinterpret_cast<uintptr_t>(alloc_rec) | REF_COUNT_BIT;
  }
}

void AllocationTracker::reallocate( size_t size ) const
{
  AllocationRecord * rec = to_alloc_rec( m_alloc_rec );

  void * the_alloc_ptr = rec->allocator->reallocate( rec->alloc_ptr, rec->alloc_size, size );

  if ( NULL != the_alloc_ptr )
  {
    *const_cast<void **>(&rec->alloc_ptr) = the_alloc_ptr;
    *const_cast<uint64_t *>(&rec->alloc_size) = size;
  }
  else {
    Impl::throw_runtime_exception( "Error: unable to reallocate allocation tracker");
  }
}


void AllocationTracker::increment_ref_count() const
{
  to_alloc_rec( m_alloc_rec )->increment_ref_count();
}


void AllocationTracker::decrement_ref_count() const
{
  AllocationRecord * alloc_rec = to_alloc_rec( m_alloc_rec );
  uint32_t the_ref_count = alloc_rec->decrement_ref_count();
  if (the_ref_count == 0u) {
    try {
      global_alloc_rec_pool.destroy( alloc_rec );
    }
    catch(...) {}
  }
}

namespace {

struct NullAllocator { static const char * name() { return "Null Allocator"; } };

}

AllocatorBase * AllocationTracker::allocator() const
{
  if (m_alloc_rec & REF_COUNT_MASK) {
    return to_alloc_rec(m_alloc_rec)->allocator;
  }
  return Allocator<NullAllocator>::singleton();
}

void * AllocationTracker::alloc_ptr()  const
{
  if (m_alloc_rec & REF_COUNT_MASK) {
    return to_alloc_rec(m_alloc_rec)->alloc_ptr;
  }
  return NULL;
}

size_t AllocationTracker::alloc_size() const
{
  if (m_alloc_rec & REF_COUNT_MASK) {
    return to_alloc_rec(m_alloc_rec)->alloc_size;
  }
  return 0u;
}

size_t AllocationTracker::ref_count()  const
{
  if (m_alloc_rec & REF_COUNT_MASK) {
    return to_alloc_rec(m_alloc_rec)->ref_count;
  }
  return 0u;
}

char const * AllocationTracker::label() const
{
  if (m_alloc_rec & REF_COUNT_MASK) {
    return to_alloc_rec(m_alloc_rec)->label;
  }
  return "[Empty Allocation Tracker]";
}

void AllocationTracker::print( std::ostream & oss) const
{
  if (m_alloc_rec & REF_COUNT_MASK) {
    to_alloc_rec(m_alloc_rec)->print(oss);
  }
  else {
    oss << label();
  }
}

bool AllocationTracker::set_attribute( AllocatorAttributeBase * attr ) const
{
  bool result = false;
  if (m_alloc_rec & REF_COUNT_MASK) {
    result = to_alloc_rec(m_alloc_rec)->set_attribute(attr);
  }
  return result;
}

AllocatorAttributeBase * AllocationTracker::attribute() const
{
  if (m_alloc_rec & REF_COUNT_MASK) {
    return to_alloc_rec(m_alloc_rec)->attribute;
  }
  return NULL;
}

void AllocationTracker::print_tracked_memory( std::ostream & out )
{
  global_alloc_rec_pool.print_memory( out );
}


AllocationTracker AllocationTracker::find( void const * ptr, AllocatorBase const * arg_allocator )
{
  AllocationRecord * alloc_rec = global_alloc_rec_pool.find(ptr, arg_allocator);

  AllocationTracker tracker;

  if ( alloc_rec != NULL )
  {
    if ( tracking_enabled() ) {
      alloc_rec->increment_ref_count();
      tracker.m_alloc_rec = reinterpret_cast<uintptr_t>(alloc_rec) | REF_COUNT_BIT;
    }
    else {
      tracker.m_alloc_rec = reinterpret_cast<uintptr_t>(alloc_rec);
    }
  }

  return tracker ;
}



//-----------------------------------------------------------------------------
// static AllocationTracker
//-----------------------------------------------------------------------------
#if defined( KOKKOS_USE_DECENTRALIZED_HOST )
namespace {

  // TODO : Detect compiler support for thread local variables
  #if defined( KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_OPENMP )
    bool g_thread_local_tracking_enabled = true;
    #pragma omp threadprivate(g_thread_local_tracking_enabled)
  #elif defined( KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_THREADS )
    __thread bool g_thread_local_tracking_enabled = true;
  #elif defined( KOKKOS_HAVE_OPENMP )
    bool g_thread_local_tracking_enabled = true;
    #pragma omp threadprivate(g_thread_local_tracking_enabled)
  #elif defined( KOKKOS_HAVE_PTHREAD )
    __thread bool g_thread_local_tracking_enabled = true;
  #elif defined( KOKKOS_HAVE_SERIAL )
      bool g_thread_local_tracking_enabled = true;
  #endif
} // unnamed namespace

void AllocationTracker::disable_tracking()
{
  g_thread_local_tracking_enabled = false;
}

void AllocationTracker::enable_tracking()
{
  g_thread_local_tracking_enabled = true;
}

bool AllocationTracker::tracking_enabled()
{
  return g_thread_local_tracking_enabled;
}
#else
namespace {
enum TrackingEnum { TRACKING_ENABLED, TRACKING_DISABLED };
volatile TrackingEnum g_tracking_enabled = TRACKING_ENABLED;
}

void AllocationTracker::disable_tracking()
{
  if ( TRACKING_ENABLED != atomic_compare_exchange( &g_tracking_enabled, TRACKING_ENABLED, TRACKING_DISABLED ) ) {
    Impl::throw_runtime_exception("Error: Tracking already disabled");
  }
}

void AllocationTracker::enable_tracking()
{
  if ( TRACKING_DISABLED != atomic_compare_exchange( &g_tracking_enabled, TRACKING_DISABLED, TRACKING_ENABLED ) ) {
    Impl::throw_runtime_exception("Error: Tracking already enabled");
  }
}

bool AllocationTracker::tracking_enabled()
{
  return g_tracking_enabled == TRACKING_ENABLED;
}
#endif


//-----------------------------------------------------------------------------
// create singleton free function
//-----------------------------------------------------------------------------
void * create_singleton(  size_t size
                        , Impl::singleton_create_function_type create_func
                        , Impl::singleton_destroy_function_type destroy_func )
{
  return global_alloc_rec_pool.create_singleton( size, create_func, destroy_func );
}

}} // namespace Kokkos::Impl

#endif /* #if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST ) */

#endif /* #if ! defined( KOKKOS_USING_EXPERIMENTAL_VIEW ) */

