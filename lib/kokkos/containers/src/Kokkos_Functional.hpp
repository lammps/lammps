#ifndef KOKKOS_FUNCTIONAL_HPP
#define KOKKOS_FUNCTIONAL_HPP

#include <Kokkos_Macros.hpp>
#include <impl/Kokkos_Functional_impl.hpp>

namespace Kokkos {

// These should work for most types

template <typename T>
struct pod_hash
{
  typedef T argument_type;
  typedef T first_argument_type;
  typedef uint32_t second_argument_type;
  typedef uint32_t result_type;

  KOKKOS_FORCEINLINE_FUNCTION
  uint32_t operator()(T const & t) const
  { return Impl::MurmurHash3_x86_32( &t, sizeof(T), 0); }

  KOKKOS_FORCEINLINE_FUNCTION
  uint32_t operator()(T const & t, uint32_t seed) const
  { return Impl::MurmurHash3_x86_32( &t, sizeof(T), seed); }
};

template <typename T>
struct pod_equal_to
{
  typedef T first_argument_type;
  typedef T second_argument_type;
  typedef bool result_type;

  KOKKOS_FORCEINLINE_FUNCTION
  bool operator()(T const & a, T const & b) const
  { return Impl::bitwise_equal(&a,&b); }
};

template <typename T>
struct pod_not_equal_to
{
  typedef T first_argument_type;
  typedef T second_argument_type;
  typedef bool result_type;

  KOKKOS_FORCEINLINE_FUNCTION
  bool operator()(T const & a, T const & b) const
  { return !Impl::bitwise_equal(&a,&b); }
};

template <typename T>
struct equal_to
{
  typedef T first_argument_type;
  typedef T second_argument_type;
  typedef bool result_type;

  KOKKOS_FORCEINLINE_FUNCTION
  bool operator()(T const & a, T const & b) const
  { return a == b; }
};

template <typename T>
struct not_equal_to
{
  typedef T first_argument_type;
  typedef T second_argument_type;
  typedef bool result_type;

  KOKKOS_FORCEINLINE_FUNCTION
  bool operator()(T const & a, T const & b) const
  { return a != b; }
};


template <typename T>
struct greater
{
  typedef T first_argument_type;
  typedef T second_argument_type;
  typedef bool result_type;

  KOKKOS_FORCEINLINE_FUNCTION
  bool operator()(T const & a, T const & b) const
  { return a > b; }
};


template <typename T>
struct less
{
  typedef T first_argument_type;
  typedef T second_argument_type;
  typedef bool result_type;

  KOKKOS_FORCEINLINE_FUNCTION
  bool operator()(T const & a, T const & b) const
  { return a < b; }
};

template <typename T>
struct greater_equal
{
  typedef T first_argument_type;
  typedef T second_argument_type;
  typedef bool result_type;

  KOKKOS_FORCEINLINE_FUNCTION
  bool operator()(T const & a, T const & b) const
  { return a >= b; }
};


template <typename T>
struct less_equal
{
  typedef T first_argument_type;
  typedef T second_argument_type;
  typedef bool result_type;

  KOKKOS_FORCEINLINE_FUNCTION
  bool operator()(T const & a, T const & b) const
  { return a <= b; }
};

} // namespace Kokkos


#endif //KOKKOS_FUNCTIONAL_HPP


