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

#ifndef KOKKOS_ALLOCATION_TRACKER_HPP
#define KOKKOS_ALLOCATION_TRACKER_HPP

#include <Kokkos_Macros.hpp>

#if ! defined( KOKKOS_USING_EXPERIMENTAL_VIEW )

#include <impl/Kokkos_Traits.hpp>
#include <impl/Kokkos_Error.hpp>

#include <stdint.h>
#include <cstdlib>
#include <string>
#include <iosfwd>

namespace Kokkos { namespace Impl {

//-----------------------------------------------------------------------------
// Create Singleton objects
//-----------------------------------------------------------------------------

typedef void * (*singleton_create_function_type)(void * buffer);
typedef void (*singleton_destroy_function_type)(void *);

void * create_singleton(  size_t size
                        , singleton_create_function_type create_func
                        , singleton_destroy_function_type destroy_func
                       );



/// class Singleton
///
/// Default construct a singleton type.  This method is used to circumvent
/// order of construction issues.  Singleton objects are destroyed after all
/// other allocations in the reverse order of their creation.
template <typename Type>
class Singleton
{
public:
  /// Get a pointer to the Singleton. Default construct the singleton if it does not already exist
  static Type * get()
  {
    static Type * singleton = NULL;
    if (singleton == NULL) {
      Impl::singleton_create_function_type  create_func = &create;
      Impl::singleton_destroy_function_type destroy_func = &destroy;
      singleton = reinterpret_cast<Type*>( Impl::create_singleton( sizeof(Type), create_func, destroy_func ) );
    }
    return singleton;
  }

private:

  /// Call the Type constructor
  static void destroy(void * ptr)
  {
    reinterpret_cast<Type*>(ptr)->~Type();
  }

  /// placement new the Type in buffer
  static void * create(void * buffer)
  {
    return new (buffer) Type();
  }
};


//-----------------------------------------------------------------------------
// AllocatorBase
//-----------------------------------------------------------------------------

/// class AllocatorBase
///
/// Abstract base class for all Allocators.
/// Allocators should be singleton objects, use Singleton<Allocator>::get to create
/// to avoid order of destruction issues
class AllocatorBase
{
public:
  /// name of the allocator
  /// used to report memory leaks
  virtual const char * name() const = 0;

  /// Allocate a buffer of size number of bytes
  virtual void* allocate(size_t size) const = 0;

  /// Deallocate a buffer with size number of bytes
  /// The pointer must have been allocated with a call to corresponding allocate
  virtual void deallocate(void * ptr, size_t size) const = 0;

  /// Changes the size of the memory block pointed to by ptr.
  /// Ptr must have been allocated with the corresponding allocate call
  /// The function may move the memory block to a new location
  /// (whose address is returned by the function).
  ///
  /// The content of the memory block is preserved up to the lesser of the new and
  /// old sizes, even if the block is moved to a new location. If the new size is larger,
  /// the value of the newly allocated portion is indeterminate.
  ///
  /// In case that ptr is a null pointer, the function behaves like allocate, assigning a
  /// new block of size bytes and returning a pointer to its beginning.
  virtual void * reallocate(void * old_ptr, size_t old_size, size_t new_size) const = 0;

  /// can a texture object be bound to the allocated memory
  virtual bool support_texture_binding() const = 0;

  /// virtual destructor
  virtual ~AllocatorBase() {}
};

/// class AllocatorAttributeBase
class AllocatorAttributeBase
{
public:
  virtual ~AllocatorAttributeBase() {}
};

//-----------------------------------------------------------------------------
// Allocator< StaticAllocator > : public AllocatorBase
//-----------------------------------------------------------------------------

// HasStaticName
template<typename T>
class HasStaticName
{
  typedef const char * (*static_method)();
  template<typename U, static_method> struct SFINAE {};
  template<typename U> static char Test(SFINAE<U, &U::name>*);
  template<typename U> static int Test(...);
public:
  enum { value = sizeof(Test<T>(0)) == sizeof(char) };
};


template <typename T>
inline
typename enable_if<HasStaticName<T>::value, const char *>::type
allocator_name()
{
  return T::name();
}

template <typename T>
inline
typename enable_if<!HasStaticName<T>::value, const char *>::type
allocator_name()
{
  return "Unnamed Allocator";
}


// HasStaticAllocate
template<typename T>
class HasStaticAllocate
{
  typedef void * (*static_method)(size_t);
  template<typename U, static_method> struct SFINAE {};
  template<typename U> static char Test(SFINAE<U, &U::allocate>*);
  template<typename U> static int Test(...);
public:
  enum { value = sizeof(Test<T>(0)) == sizeof(char) };
};

template <typename T>
inline
typename enable_if<HasStaticAllocate<T>::value, void *>::type
allocator_allocate(size_t size)
{
  return T::allocate(size);
}

template <typename T>
inline
typename enable_if<!HasStaticAllocate<T>::value, void *>::type
allocator_allocate(size_t)
{
  throw_runtime_exception(  std::string("Error: ")
                          + std::string(allocator_name<T>())
                          + std::string(" cannot allocate memory!") );
  return NULL;
}

// HasStaticDeallocate
template<typename T>
class HasStaticDeallocate
{
  typedef void (*static_method)(void *, size_t);
  template<typename U, static_method> struct SFINAE {};
  template<typename U> static char Test(SFINAE<U, &U::deallocate>*);
  template<typename U> static int Test(...);
public:
  enum { value = sizeof(Test<T>(0)) == sizeof(char) };
};

template <typename T>
inline
typename enable_if<HasStaticDeallocate<T>::value, void>::type
allocator_deallocate(void * ptr, size_t size)
{
  T::deallocate(ptr,size);
}

template <typename T>
inline
typename enable_if<!HasStaticDeallocate<T>::value, void>::type
allocator_deallocate(void *, size_t)
{
  throw_runtime_exception(  std::string("Error: ")
                          + std::string(allocator_name<T>())
                          + std::string(" cannot deallocate memory!") );
}

// HasStaticReallocate
template<typename T>
class HasStaticReallocate
{
  typedef void * (*static_method)(void *, size_t, size_t);
  template<typename U, static_method> struct SFINAE {};
  template<typename U> static char Test(SFINAE<U, &U::reallocate>*);
  template<typename U> static int Test(...);
public:
  enum { value = sizeof(Test<T>(0)) == sizeof(char) };
};

template <typename T>
inline
typename enable_if<HasStaticReallocate<T>::value, void *>::type
allocator_reallocate(void * old_ptr, size_t old_size, size_t new_size)
{
  return T::reallocate(old_ptr, old_size, new_size);
}

template <typename T>
inline
typename enable_if<!HasStaticReallocate<T>::value, void *>::type
allocator_reallocate(void *, size_t, size_t)
{
  throw_runtime_exception(  std::string("Error: ")
                          + std::string(allocator_name<T>())
                          + std::string(" cannot reallocate memory!") );
  return NULL;
}

// HasStaticReallocate
template<typename T>
class HasStaticSupportTextureBinding
{
  typedef bool (*static_method)();
  template<typename U, static_method> struct SFINAE {};
  template<typename U> static char Test(SFINAE<U, &U::support_texture_binding>*);
  template<typename U> static int Test(...);
public:
  enum { value = sizeof(Test<T>(0)) == sizeof(char) };
};

template <typename T>
inline
typename enable_if<HasStaticSupportTextureBinding<T>::value, bool>::type
allocator_support_texture_binding()
{
  return T::support_texture_binding();
}

template <typename T>
inline
typename enable_if<!HasStaticSupportTextureBinding<T>::value, bool>::type
allocator_support_texture_binding()
{
  return false;
}

template <typename T>
class Allocator : public AllocatorBase
{
public:
  virtual const char * name() const
  {
    return allocator_name<T>();
  }

  virtual void* allocate(size_t size) const
  {
    return allocator_allocate<T>(size);
  }

  virtual void deallocate(void * ptr, size_t size) const
  {
    allocator_deallocate<T>(ptr,size);
  }

  virtual void * reallocate(void * old_ptr, size_t old_size, size_t new_size) const
  {
    return allocator_reallocate<T>(old_ptr, old_size, new_size);
  }

  virtual bool support_texture_binding() const
  {
    return allocator_support_texture_binding<T>();
  }

  static AllocatorBase * singleton()
  {
    return Singleton< Allocator<T> >::get();
  }
};

//-----------------------------------------------------------------------------
// AllocationTracker
//-----------------------------------------------------------------------------

// forward declaration for friend classes
struct MallocHelper;

/// class AllocationTracker
/// Will call deallocate from the AllocatorBase when the reference count reaches 0.
/// Reference counting is disabled when the host is in parallel.
class AllocationTracker
{
  // use the least significant bit of the AllocationRecord pointer to indicate if the
  // AllocationTracker should reference count
  enum {
     REF_COUNT_BIT = static_cast<uintptr_t>(1)
   , REF_COUNT_MASK = ~static_cast<uintptr_t>(1)
  };

public:

  /// Find an AllocationTracker such that
  /// alloc_ptr <= ptr < alloc_ptr + alloc_size
  /// O(n) where n is the number of tracked allocations.
  template <typename StaticAllocator>
  static AllocationTracker find( void const * ptr )
  {
    return find( ptr, Allocator<StaticAllocator>::singleton() );
  }


  /// Pretty print all the currently tracked memory
  static void print_tracked_memory( std::ostream & out );

  /// Default constructor
  KOKKOS_INLINE_FUNCTION
  AllocationTracker()
    : m_alloc_rec(0)
  {}

  /// Create a AllocationTracker
  ///
  /// Start reference counting the alloc_ptr.
  /// When the reference count reachs 0 the allocator deallocate method
  /// will be call with the given size.  The alloc_ptr should have been
  /// allocated with the allocator's allocate method.
  ///
  /// If arg_allocator == NULL OR arg_alloc_ptr == NULL OR size == 0
  /// do nothing
  template <typename StaticAllocator>
  AllocationTracker(  StaticAllocator const &
                    , void * arg_alloc_ptr
                    , size_t arg_alloc_size
                    , const std::string & arg_label = std::string("") )
    : m_alloc_rec(0)
  {
    AllocatorBase * arg_allocator = Allocator<StaticAllocator>::singleton();
    initalize( arg_allocator, arg_alloc_ptr, arg_alloc_size, arg_label);
  }

  /// Create a AllocationTracker
  ///
  /// Start reference counting the alloc_ptr.
  /// When the reference count reachs 0 the allocator deallocate method
  /// will be call with the given size.  The alloc_ptr should have been
  /// allocated with the allocator's allocate method.
  ///
  /// If arg_allocator == NULL OR arg_alloc_ptr == NULL OR size == 0
  /// do nothing
  template <typename StaticAllocator>
  AllocationTracker(  StaticAllocator const &
                    , size_t arg_alloc_size
                    , const std::string & arg_label = std::string("")
                   )
    : m_alloc_rec(0)
  {
    AllocatorBase * arg_allocator = Allocator<StaticAllocator>::singleton();
    void * arg_alloc_ptr = arg_allocator->allocate( arg_alloc_size );

    initalize( arg_allocator, arg_alloc_ptr, arg_alloc_size, arg_label);
  }

  /// Copy an AllocatorTracker
  KOKKOS_INLINE_FUNCTION
  AllocationTracker( const AllocationTracker & rhs )
    : m_alloc_rec( rhs.m_alloc_rec)
  {
#if !defined( __CUDA_ARCH__ )
    if ( rhs.ref_counting() && tracking_enabled() ) {
      increment_ref_count();
    }
    else {
      m_alloc_rec = m_alloc_rec & REF_COUNT_MASK;
    }
#else
    m_alloc_rec = m_alloc_rec & REF_COUNT_MASK;
#endif
  }

  /// Copy an AllocatorTracker
  /// Decrement the reference count of the current tracker if necessary
  KOKKOS_INLINE_FUNCTION
  AllocationTracker & operator=( const AllocationTracker & rhs )
  {
    if (this != &rhs) {
#if !defined( __CUDA_ARCH__ )
      if ( ref_counting() ) {
        decrement_ref_count();
      }

      m_alloc_rec = rhs.m_alloc_rec;

      if ( rhs.ref_counting() && tracking_enabled() ) {
        increment_ref_count();
      }
      else {
        m_alloc_rec = m_alloc_rec & REF_COUNT_MASK;
      }
#else
      m_alloc_rec = rhs.m_alloc_rec & REF_COUNT_MASK;
#endif
    }

    return * this;
  }

  /// Destructor
  /// Decrement the reference count if necessary
  KOKKOS_INLINE_FUNCTION
  ~AllocationTracker()
  {
#if !defined( __CUDA_ARCH__ )
    if ( ref_counting() ) {
      decrement_ref_count();
    }
#endif
  }

  /// Is the tracker valid?
  KOKKOS_INLINE_FUNCTION
  bool is_valid() const
  {
    return (m_alloc_rec & REF_COUNT_MASK);
  }



  /// clear the tracker
  KOKKOS_INLINE_FUNCTION
  void clear()
  {
#if !defined( __CUDA_ARCH__ )
    if ( ref_counting() ) {
      decrement_ref_count();
    }
#endif
    m_alloc_rec = 0;
  }

  /// is this tracker currently counting allocations?
  KOKKOS_INLINE_FUNCTION
  bool ref_counting() const
  {
    return (m_alloc_rec & REF_COUNT_BIT);
  }

  AllocatorBase * allocator() const;

  /// pointer to the allocated memory
  void * alloc_ptr()  const;

  /// size in bytes of the allocated memory
  size_t alloc_size() const;

  /// the current reference count
  size_t ref_count()  const;

  /// the label given to the allocation
  char const * label() const;

  /// pretty print all the tracker's information to the std::ostream
  void print( std::ostream & oss) const;


  /// set an attribute ptr on the allocation record
  /// the arg_attribute pointer will be deleted when the record is destroyed
  /// the attribute ptr can only be set once
  bool set_attribute( AllocatorAttributeBase * arg_attribute) const;

  /// get the attribute ptr from the allocation record
  AllocatorAttributeBase * attribute() const;


  /// reallocate the memory tracked by this allocation
  /// NOT thread-safe
  void reallocate( size_t size ) const;

  static void disable_tracking();
  static void enable_tracking();
  static bool tracking_enabled();

private:

  static AllocationTracker find( void const * ptr, AllocatorBase const * arg_allocator );

  void initalize(  AllocatorBase * arg_allocator
                 , void * arg_alloc_ptr
                 , size_t arg_alloc_size
                 , std::string const & label );

  void increment_ref_count() const;
  void decrement_ref_count() const;

  friend struct Impl::MallocHelper;

  uintptr_t m_alloc_rec;
};

}} // namespace Kokkos::Impl

#endif /* #if ! defined( KOKKOS_USING_EXPERIMENTAL_VIEW ) */

#endif //KOKKOS_ALLOCATION_TRACKER_HPP

