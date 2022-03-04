#include <Kokkos_NumericTraits.hpp>

// NOTE These out-of class definitions are only required with C++14.  Since
// C++17, a static data member declared constrexpr is impllictly inline.

#if !defined(KOKKOS_ENABLE_CXX17)
namespace Kokkos {
namespace Experimental {
namespace Impl {
#define OUT_OF_CLASS_DEFINTION_FLOATING_POINT(TRAIT) \
  constexpr float TRAIT##_helper<float>::value;      \
  constexpr double TRAIT##_helper<double>::value;    \
  constexpr long double TRAIT##_helper<long double>::value

#define OUT_OF_CLASS_DEFINTION_INTEGRAL(TRAIT)                          \
  constexpr bool TRAIT##_helper<bool>::value;                           \
  constexpr char TRAIT##_helper<char>::value;                           \
  constexpr signed char TRAIT##_helper<signed char>::value;             \
  constexpr unsigned char TRAIT##_helper<unsigned char>::value;         \
  constexpr short TRAIT##_helper<short>::value;                         \
  constexpr unsigned short TRAIT##_helper<unsigned short>::value;       \
  constexpr int TRAIT##_helper<int>::value;                             \
  constexpr unsigned int TRAIT##_helper<unsigned int>::value;           \
  constexpr long int TRAIT##_helper<long int>::value;                   \
  constexpr unsigned long int TRAIT##_helper<unsigned long int>::value; \
  constexpr long long int TRAIT##_helper<long long int>::value;         \
  constexpr unsigned long long int TRAIT##_helper<unsigned long long int>::value

#define OUT_OF_CLASS_DEFINTION_FLOATING_POINT_2(TRAIT) \
  constexpr int TRAIT##_helper<float>::value;          \
  constexpr int TRAIT##_helper<double>::value;         \
  constexpr int TRAIT##_helper<long double>::value

#define OUT_OF_CLASS_DEFINTION_INTEGRAL_2(TRAIT)          \
  constexpr int TRAIT##_helper<bool>::value;              \
  constexpr int TRAIT##_helper<char>::value;              \
  constexpr int TRAIT##_helper<signed char>::value;       \
  constexpr int TRAIT##_helper<unsigned char>::value;     \
  constexpr int TRAIT##_helper<short>::value;             \
  constexpr int TRAIT##_helper<unsigned short>::value;    \
  constexpr int TRAIT##_helper<int>::value;               \
  constexpr int TRAIT##_helper<unsigned int>::value;      \
  constexpr int TRAIT##_helper<long int>::value;          \
  constexpr int TRAIT##_helper<unsigned long int>::value; \
  constexpr int TRAIT##_helper<long long int>::value;     \
  constexpr int TRAIT##_helper<unsigned long long int>::value

OUT_OF_CLASS_DEFINTION_FLOATING_POINT(infinity);
OUT_OF_CLASS_DEFINTION_FLOATING_POINT(epsilon);
OUT_OF_CLASS_DEFINTION_FLOATING_POINT(round_error);
OUT_OF_CLASS_DEFINTION_FLOATING_POINT(norm_min);

OUT_OF_CLASS_DEFINTION_INTEGRAL(finite_min);
OUT_OF_CLASS_DEFINTION_FLOATING_POINT(finite_min);
OUT_OF_CLASS_DEFINTION_INTEGRAL(finite_max);
OUT_OF_CLASS_DEFINTION_FLOATING_POINT(finite_max);

OUT_OF_CLASS_DEFINTION_INTEGRAL_2(digits);
OUT_OF_CLASS_DEFINTION_FLOATING_POINT_2(digits);
OUT_OF_CLASS_DEFINTION_INTEGRAL_2(digits10);
OUT_OF_CLASS_DEFINTION_FLOATING_POINT_2(digits10);
OUT_OF_CLASS_DEFINTION_FLOATING_POINT_2(max_digits10);
OUT_OF_CLASS_DEFINTION_INTEGRAL_2(radix);
OUT_OF_CLASS_DEFINTION_FLOATING_POINT_2(radix);

OUT_OF_CLASS_DEFINTION_FLOATING_POINT_2(min_exponent);
OUT_OF_CLASS_DEFINTION_FLOATING_POINT_2(min_exponent10);
OUT_OF_CLASS_DEFINTION_FLOATING_POINT_2(max_exponent);
OUT_OF_CLASS_DEFINTION_FLOATING_POINT_2(max_exponent10);
}  // namespace Impl
}  // namespace Experimental
}  // namespace Kokkos
#endif
