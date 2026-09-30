// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#include <Kokkos_Macros.hpp>
static_assert(false,
              "Including non-public Kokkos header files is not allowed.");
#endif
#ifndef KOKKOS_BUILTIN_REDUCERS_HPP
#define KOKKOS_BUILTIN_REDUCERS_HPP

#include <Kokkos_ReductionIdentity.hpp>
#include <Kokkos_View.hpp>
#include <type_traits>

namespace Kokkos {

// \brief Class offering functionalities common to all reducers
//
// In order to be a valid reducer, a class must implement the functions and
// define the types offered in this class.
// To facilitate implementation, a new reducer class can simply inherit from
// BaseReducer.
namespace Impl {
template <class Scalar, class Space>
struct BaseReducer {
 public:
  // Following types need to be available for the reducer to be valid
  using value_type       = std::remove_cv_t<Scalar>;
  using result_view_type = Kokkos::View<value_type, Space>;

  static_assert(!std::is_pointer_v<value_type> && !std::is_array_v<value_type>);

 protected:
  // Contains the value of the reduction
  result_view_type value;
  // Whether the reducer returns its value through a Kokkos::View or a scalar
  bool references_scalar_v;

 public:
  // Construct from a scalar value
  KOKKOS_INLINE_FUNCTION
  BaseReducer(value_type& value_) : value(&value_), references_scalar_v(true) {}

  // Construct from a View
  KOKKOS_INLINE_FUNCTION
  BaseReducer(const result_view_type& value_)
      : value(value_), references_scalar_v(false) {}

  // Reducers also need to implement the two following functions:
  // KOKKOS_INLINE_FUNCTION
  // void join(value_type& dest, const value_type& src) const {
  //    // Do the reduction operation here
  // }

  // KOKKOS_INLINE_FUNCTION
  // void init(value_type& val) const {
  //   // Return the neutral value for the reduction operation, for instance
  //   // FLOAT_MIN if searching for the max, 0 for a sum or 1 for a product).
  // }
  //

  // Needed accessors
  KOKKOS_INLINE_FUNCTION
  value_type& reference() const { return *value.data(); }

  KOKKOS_INLINE_FUNCTION
  result_view_type view() const { return value; }

  KOKKOS_INLINE_FUNCTION
  bool references_scalar() const { return references_scalar_v; }
};
}  // namespace Impl

template <class Scalar, class Space>
struct Sum : Impl::BaseReducer<Scalar, Space> {
 private:
  using parent_type = Impl::BaseReducer<Scalar, Space>;

 public:
  using reducer    = Sum<Scalar, Space>;
  using value_type = typename parent_type::value_type;

  // Inherit constructors
  using parent_type::parent_type;

  KOKKOS_INLINE_FUNCTION
  void join(value_type& dest, const value_type& src) const { dest += src; }

  KOKKOS_INLINE_FUNCTION
  void init(value_type& val) const {
    val = reduction_identity<value_type>::sum();
  }
};

template <typename Scalar, typename... Properties>
KOKKOS_DEDUCTION_GUIDE Sum(View<Scalar, Properties...> const&)
    -> Sum<Scalar, typename View<Scalar, Properties...>::memory_space>;

template <class Scalar, class Space>
struct Prod : Impl::BaseReducer<Scalar, Space> {
 private:
  using parent_type = Impl::BaseReducer<Scalar, Space>;

 public:
  using reducer    = Prod<Scalar, Space>;
  using value_type = typename parent_type::value_type;

  // Inherit constructors
  using parent_type::parent_type;

  KOKKOS_INLINE_FUNCTION
  void join(value_type& dest, const value_type& src) const { dest *= src; }

  KOKKOS_INLINE_FUNCTION
  void init(value_type& val) const {
    val = reduction_identity<value_type>::prod();
  }
};

template <typename Scalar, typename... Properties>
KOKKOS_DEDUCTION_GUIDE Prod(View<Scalar, Properties...> const&)
    -> Prod<Scalar, typename View<Scalar, Properties...>::memory_space>;

template <class Scalar, class Space>
struct Min : Impl::BaseReducer<Scalar, Space> {
 private:
  using parent_type = Impl::BaseReducer<Scalar, Space>;

 public:
  using reducer    = Min<Scalar, Space>;
  using value_type = typename parent_type::value_type;

  // Inherit constructors
  using parent_type::parent_type;

  KOKKOS_INLINE_FUNCTION
  void join(value_type& dest, const value_type& src) const {
    if (src < dest) dest = src;
  }

  KOKKOS_INLINE_FUNCTION
  void init(value_type& val) const {
    val = reduction_identity<value_type>::min();
  }
};

template <typename Scalar, typename... Properties>
KOKKOS_DEDUCTION_GUIDE Min(View<Scalar, Properties...> const&)
    -> Min<Scalar, typename View<Scalar, Properties...>::memory_space>;

template <class Scalar, class Space>
struct Max : Impl::BaseReducer<Scalar, Space> {
 private:
  using parent_type = Impl::BaseReducer<Scalar, Space>;

 public:
  using reducer    = Max<Scalar, Space>;
  using value_type = typename parent_type::value_type;

  // Inherit constructors
  using parent_type::parent_type;

  KOKKOS_INLINE_FUNCTION
  void join(value_type& dest, const value_type& src) const {
    if (src > dest) dest = src;
  }

  KOKKOS_INLINE_FUNCTION
  void init(value_type& val) const {
    val = reduction_identity<value_type>::max();
  }
};

template <typename Scalar, typename... Properties>
KOKKOS_DEDUCTION_GUIDE Max(View<Scalar, Properties...> const&)
    -> Max<Scalar, typename View<Scalar, Properties...>::memory_space>;

template <class Scalar, class Space>
struct LAnd : Impl::BaseReducer<Scalar, Space> {
 private:
  using parent_type = Impl::BaseReducer<Scalar, Space>;

 public:
  using reducer    = LAnd<Scalar, Space>;
  using value_type = typename parent_type::value_type;

  // Inherit constructors
  using parent_type::parent_type;

  KOKKOS_INLINE_FUNCTION
  void join(value_type& dest, const value_type& src) const {
    dest = dest && src;
  }

  KOKKOS_INLINE_FUNCTION
  void init(value_type& val) const {
    val = reduction_identity<value_type>::land();
  }
};

template <typename Scalar, typename... Properties>
KOKKOS_DEDUCTION_GUIDE LAnd(View<Scalar, Properties...> const&)
    -> LAnd<Scalar, typename View<Scalar, Properties...>::memory_space>;

template <class Scalar, class Space>
struct LOr : Impl::BaseReducer<Scalar, Space> {
 private:
  using parent_type = Impl::BaseReducer<Scalar, Space>;

 public:
  using reducer    = LOr<Scalar, Space>;
  using value_type = typename parent_type::value_type;

  // Inherit constructors
  using parent_type::parent_type;

  KOKKOS_INLINE_FUNCTION
  void join(value_type& dest, const value_type& src) const {
    dest = dest || src;
  }

  KOKKOS_INLINE_FUNCTION
  void init(value_type& val) const {
    val = reduction_identity<value_type>::lor();
  }
};

template <typename Scalar, typename... Properties>
KOKKOS_DEDUCTION_GUIDE LOr(View<Scalar, Properties...> const&)
    -> LOr<Scalar, typename View<Scalar, Properties...>::memory_space>;

template <class Scalar, class Space>
struct BAnd : Impl::BaseReducer<Scalar, Space> {
 private:
  using parent_type = Impl::BaseReducer<Scalar, Space>;

 public:
  using reducer    = BAnd<Scalar, Space>;
  using value_type = typename parent_type::value_type;

  // Inherit constructors
  using parent_type::parent_type;

  KOKKOS_INLINE_FUNCTION
  void join(value_type& dest, const value_type& src) const {
    dest = dest & src;
  }

  KOKKOS_INLINE_FUNCTION
  void init(value_type& val) const {
    val = reduction_identity<value_type>::band();
  }
};

template <typename Scalar, typename... Properties>
KOKKOS_DEDUCTION_GUIDE BAnd(View<Scalar, Properties...> const&)
    -> BAnd<Scalar, typename View<Scalar, Properties...>::memory_space>;

template <class Scalar, class Space>
struct BOr : Impl::BaseReducer<Scalar, Space> {
 private:
  using parent_type = Impl::BaseReducer<Scalar, Space>;

 public:
  using reducer    = BOr<Scalar, Space>;
  using value_type = typename parent_type::value_type;

  // Inherit constructors
  using parent_type::parent_type;

  KOKKOS_INLINE_FUNCTION
  void join(value_type& dest, const value_type& src) const {
    dest = dest | src;
  }

  KOKKOS_INLINE_FUNCTION
  void init(value_type& val) const {
    val = reduction_identity<value_type>::bor();
  }
};

template <typename Scalar, typename... Properties>
KOKKOS_DEDUCTION_GUIDE BOr(View<Scalar, Properties...> const&)
    -> BOr<Scalar, typename View<Scalar, Properties...>::memory_space>;

template <class Scalar, class Index>
struct ValLocScalar {
  Scalar val;
  Index loc;
};

template <class Scalar, class Index, class Space>
struct MinLoc
    : Impl::BaseReducer<
          ValLocScalar<std::remove_cv_t<Scalar>, std::remove_cv_t<Index>>,
          Space> {
 private:
  using scalar_type = std::remove_cv_t<Scalar>;
  using index_type  = std::remove_cv_t<Index>;
  using parent_type =
      Impl::BaseReducer<ValLocScalar<scalar_type, index_type>, Space>;

  static_assert(!std::is_pointer_v<scalar_type> &&
                !std::is_array_v<scalar_type>);

 public:
  using reducer    = MinLoc<Scalar, Index, Space>;
  using value_type = typename parent_type::value_type;

  // Inherit constructors
  using parent_type::parent_type;

  KOKKOS_INLINE_FUNCTION
  void join(value_type& dest, const value_type& src) const {
    if (src.val < dest.val) {
      dest = src;
    } else if (src.val == dest.val &&
               dest.loc == reduction_identity<index_type>::min()) {
      dest.loc = src.loc;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void init(value_type& val) const {
    val.val = reduction_identity<scalar_type>::min();
    val.loc = reduction_identity<index_type>::min();
  }
};

template <typename Scalar, typename Index, typename... Properties>
KOKKOS_DEDUCTION_GUIDE
MinLoc(View<ValLocScalar<Scalar, Index>, Properties...> const&) -> MinLoc<
    Scalar, Index,
    typename View<ValLocScalar<Scalar, Index>, Properties...>::memory_space>;

template <class Scalar, class Index, class Space>
struct MaxLoc
    : Impl::BaseReducer<
          ValLocScalar<std::remove_cv_t<Scalar>, std::remove_cv_t<Index>>,
          Space> {
 private:
  using scalar_type = std::remove_cv_t<Scalar>;
  using index_type  = std::remove_cv_t<Index>;
  using parent_type =
      Impl::BaseReducer<ValLocScalar<scalar_type, index_type>, Space>;

  static_assert(!std::is_pointer_v<scalar_type> &&
                !std::is_array_v<scalar_type>);

 public:
  using value_type = typename parent_type::value_type;

  using reducer = MaxLoc<Scalar, Index, Space>;

  // Inherit constructors
  using parent_type::parent_type;

  KOKKOS_INLINE_FUNCTION
  void join(value_type& dest, const value_type& src) const {
    if (src.val > dest.val) {
      dest = src;
    } else if (src.val == dest.val &&
               dest.loc == reduction_identity<index_type>::min()) {
      dest.loc = src.loc;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void init(value_type& val) const {
    val.val = reduction_identity<scalar_type>::max();
    val.loc = reduction_identity<index_type>::min();
  }
};

template <typename Scalar, typename Index, typename... Properties>
KOKKOS_DEDUCTION_GUIDE
MaxLoc(View<ValLocScalar<Scalar, Index>, Properties...> const&) -> MaxLoc<
    Scalar, Index,
    typename View<ValLocScalar<Scalar, Index>, Properties...>::memory_space>;

template <class Scalar>
struct MinMaxScalar {
  Scalar min_val, max_val;
};

template <class Scalar, class Space>
struct MinMax
    : Impl::BaseReducer<MinMaxScalar<std::remove_cv_t<Scalar>>, Space> {
 private:
  using scalar_type = std::remove_cv_t<Scalar>;
  using parent_type =
      Impl::BaseReducer<MinMaxScalar<std::remove_cv_t<Scalar>>, Space>;

  static_assert(!std::is_pointer_v<scalar_type> &&
                !std::is_array_v<scalar_type>);

 public:
  using value_type = typename parent_type::value_type;

  using reducer = MinMax<Scalar, Space>;

  // Inherit constructors
  using parent_type::parent_type;

  KOKKOS_INLINE_FUNCTION
  void join(value_type& dest, const value_type& src) const {
    if (src.min_val < dest.min_val) {
      dest.min_val = src.min_val;
    }
    if (src.max_val > dest.max_val) {
      dest.max_val = src.max_val;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void init(value_type& val) const {
    val.max_val = reduction_identity<scalar_type>::max();
    val.min_val = reduction_identity<scalar_type>::min();
  }
};

template <typename Scalar, typename... Properties>
KOKKOS_DEDUCTION_GUIDE MinMax(View<MinMaxScalar<Scalar>, Properties...> const&)
    -> MinMax<Scalar,
              typename View<MinMaxScalar<Scalar>, Properties...>::memory_space>;

template <class Scalar, class Index>
struct MinMaxLocScalar {
  Scalar min_val, max_val;
  Index min_loc, max_loc;
};

template <class Scalar, class Index, class Space>
struct MinMaxLoc
    : Impl::BaseReducer<
          MinMaxLocScalar<std::remove_cv_t<Scalar>, std::remove_cv_t<Index>>,
          Space> {
 private:
  using scalar_type = std::remove_cv_t<Scalar>;
  using index_type  = std::remove_cv_t<Index>;
  using parent_type =
      Impl::BaseReducer<MinMaxLocScalar<scalar_type, index_type>, Space>;

  static_assert(!std::is_pointer_v<scalar_type> &&
                !std::is_array_v<scalar_type>);

 public:
  using reducer    = MinMaxLoc<Scalar, Index, Space>;
  using value_type = typename parent_type::value_type;

  // Inherit constructors
  using parent_type::parent_type;

  KOKKOS_INLINE_FUNCTION
  void join(value_type& dest, const value_type& src) const {
    if (src.min_val < dest.min_val) {
      dest.min_val = src.min_val;
      dest.min_loc = src.min_loc;
    } else if (dest.min_val == src.min_val &&
               dest.min_loc == reduction_identity<index_type>::min()) {
      dest.min_loc = src.min_loc;
    }
    if (src.max_val > dest.max_val) {
      dest.max_val = src.max_val;
      dest.max_loc = src.max_loc;
    } else if (dest.max_val == src.max_val &&
               dest.max_loc == reduction_identity<index_type>::min()) {
      dest.max_loc = src.max_loc;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void init(value_type& val) const {
    val.max_val = reduction_identity<scalar_type>::max();
    val.min_val = reduction_identity<scalar_type>::min();
    val.max_loc = reduction_identity<index_type>::min();
    val.min_loc = reduction_identity<index_type>::min();
  }
};

template <typename Scalar, typename Index, typename... Properties>
KOKKOS_DEDUCTION_GUIDE MinMaxLoc(
    View<MinMaxLocScalar<Scalar, Index>, Properties...> const&)
    -> MinMaxLoc<Scalar, Index,
                 typename View<MinMaxLocScalar<Scalar, Index>,
                               Properties...>::memory_space>;

// --------------------------------------------------
// reducers added to support std algorithms
// --------------------------------------------------

//
// MaxFirstLoc
//
template <class Scalar, class Index, class Space>
struct MaxFirstLoc
    : Impl::BaseReducer<
          ValLocScalar<std::remove_cv_t<Scalar>, std::remove_cv_t<Index>>,
          Space> {
 private:
  using scalar_type = std::remove_cv_t<Scalar>;
  using index_type  = std::remove_cv_t<Index>;
  using parent_type =
      Impl::BaseReducer<ValLocScalar<scalar_type, index_type>, Space>;

  static_assert(!std::is_pointer_v<scalar_type> &&
                !std::is_array_v<scalar_type>);
  static_assert(std::is_integral_v<index_type>);

 public:
  using reducer    = MaxFirstLoc<Scalar, Index, Space>;
  using value_type = typename parent_type::value_type;

  // Inherit constructors
  using parent_type::parent_type;

  KOKKOS_INLINE_FUNCTION
  void join(value_type& dest, const value_type& src) const {
    if (dest.val < src.val) {
      dest = src;
    } else if (!(src.val < dest.val)) {
      dest.loc = (src.loc < dest.loc) ? src.loc : dest.loc;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void init(value_type& val) const {
    val.val = reduction_identity<scalar_type>::max();
    val.loc = reduction_identity<index_type>::min();
  }
};

template <typename Scalar, typename Index, typename... Properties>
KOKKOS_DEDUCTION_GUIDE MaxFirstLoc(
    View<ValLocScalar<Scalar, Index>, Properties...> const&)
    -> MaxFirstLoc<Scalar, Index,
                   typename View<ValLocScalar<Scalar, Index>,
                                 Properties...>::memory_space>;

//
// MaxFirstLocCustomComparator
// recall that comp(a,b) returns true is a < b
//
template <class Scalar, class Index, class ComparatorType, class Space>
struct MaxFirstLocCustomComparator
    : Impl::BaseReducer<
          ValLocScalar<std::remove_cv_t<Scalar>, std::remove_cv_t<Index>>,
          Space> {
 private:
  using scalar_type = std::remove_cv_t<Scalar>;
  using index_type  = std::remove_cv_t<Index>;
  using parent_type =
      Impl::BaseReducer<ValLocScalar<scalar_type, index_type>, Space>;

  static_assert(!std::is_pointer_v<scalar_type> &&
                !std::is_array_v<scalar_type>);
  static_assert(std::is_integral_v<index_type>);

 public:
  using reducer =
      MaxFirstLocCustomComparator<Scalar, Index, ComparatorType, Space>;
  using value_type       = typename parent_type::value_type;
  using result_view_type = typename parent_type::result_view_type;

 private:
  ComparatorType m_comp;

 public:
  KOKKOS_INLINE_FUNCTION
  MaxFirstLocCustomComparator(value_type& value_, ComparatorType comp_)
      : parent_type(value_), m_comp(comp_) {}

  KOKKOS_INLINE_FUNCTION
  MaxFirstLocCustomComparator(const result_view_type& value_,
                              ComparatorType comp_)
      : parent_type(value_), m_comp(comp_) {}

  KOKKOS_INLINE_FUNCTION
  void join(value_type& dest, const value_type& src) const {
    if (m_comp(dest.val, src.val)) {
      dest = src;
    } else if (!m_comp(src.val, dest.val)) {
      dest.loc = (src.loc < dest.loc) ? src.loc : dest.loc;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void init(value_type& val) const {
    val.val = reduction_identity<scalar_type>::max();
    val.loc = reduction_identity<index_type>::min();
  }
};

template <typename Scalar, typename Index, typename ComparatorType,
          typename... Properties>
KOKKOS_DEDUCTION_GUIDE MaxFirstLocCustomComparator(
    View<ValLocScalar<Scalar, Index>, Properties...> const&, ComparatorType)
    -> MaxFirstLocCustomComparator<Scalar, Index, ComparatorType,
                                   typename View<ValLocScalar<Scalar, Index>,
                                                 Properties...>::memory_space>;

//
// MinFirstLoc
//
template <class Scalar, class Index, class Space>
struct MinFirstLoc
    : Impl::BaseReducer<
          ValLocScalar<std::remove_cv_t<Scalar>, std::remove_cv_t<Index>>,
          Space> {
 private:
  using scalar_type = std::remove_cv_t<Scalar>;
  using index_type  = std::remove_cv_t<Index>;
  using parent_type =
      Impl::BaseReducer<ValLocScalar<scalar_type, index_type>, Space>;

  static_assert(!std::is_pointer_v<scalar_type> &&
                !std::is_array_v<scalar_type>);
  static_assert(std::is_integral_v<index_type>);

 public:
  using reducer    = MinFirstLoc<Scalar, Index, Space>;
  using value_type = typename parent_type::value_type;

  // Inherit constructors
  using parent_type::parent_type;

  KOKKOS_INLINE_FUNCTION
  void join(value_type& dest, const value_type& src) const {
    if (src.val < dest.val) {
      dest = src;
    } else if (!(dest.val < src.val)) {
      dest.loc = (src.loc < dest.loc) ? src.loc : dest.loc;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void init(value_type& val) const {
    val.val = reduction_identity<scalar_type>::min();
    val.loc = reduction_identity<index_type>::min();
  }
};

template <typename Scalar, typename Index, typename... Properties>
KOKKOS_DEDUCTION_GUIDE MinFirstLoc(
    View<ValLocScalar<Scalar, Index>, Properties...> const&)
    -> MinFirstLoc<Scalar, Index,
                   typename View<ValLocScalar<Scalar, Index>,
                                 Properties...>::memory_space>;

//
// MinFirstLocCustomComparator
// recall that comp(a,b) returns true is a < b
//
template <class Scalar, class Index, class ComparatorType, class Space>
struct MinFirstLocCustomComparator
    : Impl::BaseReducer<
          ValLocScalar<std::remove_cv_t<Scalar>, std::remove_cv_t<Index>>,
          Space> {
 private:
  using scalar_type = std::remove_cv_t<Scalar>;
  using index_type  = std::remove_cv_t<Index>;
  using parent_type =
      Impl::BaseReducer<ValLocScalar<scalar_type, index_type>, Space>;

  static_assert(!std::is_pointer_v<scalar_type> &&
                !std::is_array_v<scalar_type>);
  static_assert(std::is_integral_v<index_type>);

 public:
  using reducer =
      MinFirstLocCustomComparator<Scalar, Index, ComparatorType, Space>;
  using value_type       = typename parent_type::value_type;
  using result_view_type = typename parent_type::result_view_type;

 private:
  ComparatorType m_comp;

 public:
  KOKKOS_INLINE_FUNCTION
  MinFirstLocCustomComparator(value_type& value_, ComparatorType comp_)
      : parent_type(value_), m_comp(comp_) {}

  KOKKOS_INLINE_FUNCTION
  MinFirstLocCustomComparator(const result_view_type& value_,
                              ComparatorType comp_)
      : parent_type(value_), m_comp(comp_) {}

  KOKKOS_INLINE_FUNCTION
  void join(value_type& dest, const value_type& src) const {
    if (m_comp(src.val, dest.val)) {
      dest = src;
    } else if (!m_comp(dest.val, src.val)) {
      dest.loc = (src.loc < dest.loc) ? src.loc : dest.loc;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void init(value_type& val) const {
    val.val = reduction_identity<scalar_type>::min();
    val.loc = reduction_identity<index_type>::min();
  }
};

template <typename Scalar, typename Index, typename ComparatorType,
          typename... Properties>
KOKKOS_DEDUCTION_GUIDE MinFirstLocCustomComparator(
    View<ValLocScalar<Scalar, Index>, Properties...> const&, ComparatorType)
    -> MinFirstLocCustomComparator<Scalar, Index, ComparatorType,
                                   typename View<ValLocScalar<Scalar, Index>,
                                                 Properties...>::memory_space>;

//
// MinMaxFirstLastLoc
//
template <class Scalar, class Index, class Space>
struct MinMaxFirstLastLoc
    : Impl::BaseReducer<
          MinMaxLocScalar<std::remove_cv_t<Scalar>, std::remove_cv_t<Index>>,
          Space> {
 private:
  using scalar_type = std::remove_cv_t<Scalar>;
  using index_type  = std::remove_cv_t<Index>;
  using parent_type =
      Impl::BaseReducer<MinMaxLocScalar<scalar_type, index_type>, Space>;

  static_assert(!std::is_pointer_v<scalar_type> &&
                !std::is_array_v<scalar_type>);

 public:
  using reducer    = MinMaxFirstLastLoc<Scalar, Index, Space>;
  using value_type = typename parent_type::value_type;

  // Inherit constructors
  using parent_type::parent_type;

  KOKKOS_INLINE_FUNCTION
  void join(value_type& dest, const value_type& src) const {
    if (src.min_val < dest.min_val) {
      dest.min_val = src.min_val;
      dest.min_loc = src.min_loc;
    } else if (!(dest.min_val < src.min_val)) {
      dest.min_loc = (src.min_loc < dest.min_loc) ? src.min_loc : dest.min_loc;
    }

    if (dest.max_val < src.max_val) {
      dest.max_val = src.max_val;
      dest.max_loc = src.max_loc;
    } else if (!(src.max_val < dest.max_val)) {
      dest.max_loc = (src.max_loc > dest.max_loc) ? src.max_loc : dest.max_loc;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void init(value_type& val) const {
    val.max_val = ::Kokkos::reduction_identity<scalar_type>::max();
    val.min_val = ::Kokkos::reduction_identity<scalar_type>::min();
    val.max_loc = ::Kokkos::reduction_identity<index_type>::max();
    val.min_loc = ::Kokkos::reduction_identity<index_type>::min();
  }
};

template <typename Scalar, typename Index, typename... Properties>
KOKKOS_DEDUCTION_GUIDE MinMaxFirstLastLoc(
    View<MinMaxLocScalar<Scalar, Index>, Properties...> const&)
    -> MinMaxFirstLastLoc<Scalar, Index,
                          typename View<MinMaxLocScalar<Scalar, Index>,
                                        Properties...>::memory_space>;

//
// MinMaxFirstLastLocCustomComparator
// recall that comp(a,b) returns true is a < b
//
template <class Scalar, class Index, class ComparatorType, class Space>
struct MinMaxFirstLastLocCustomComparator
    : Impl::BaseReducer<
          MinMaxLocScalar<std::remove_cv_t<Scalar>, std::remove_cv_t<Index>>,
          Space> {
 private:
  using scalar_type = std::remove_cv_t<Scalar>;
  using index_type  = std::remove_cv_t<Index>;
  using parent_type =
      Impl::BaseReducer<MinMaxLocScalar<scalar_type, index_type>, Space>;

  static_assert(!std::is_pointer_v<scalar_type> &&
                !std::is_array_v<scalar_type>);

 public:
  using reducer =
      MinMaxFirstLastLocCustomComparator<Scalar, Index, ComparatorType, Space>;
  using value_type       = typename parent_type::value_type;
  using result_view_type = typename parent_type::result_view_type;

 private:
  ComparatorType m_comp;

 public:
  KOKKOS_INLINE_FUNCTION
  MinMaxFirstLastLocCustomComparator(value_type& value_, ComparatorType comp_)
      : parent_type(value_), m_comp(comp_) {}

  KOKKOS_INLINE_FUNCTION
  MinMaxFirstLastLocCustomComparator(const result_view_type& value_,
                                     ComparatorType comp_)
      : parent_type(value_), m_comp(comp_) {}

  KOKKOS_INLINE_FUNCTION
  void join(value_type& dest, const value_type& src) const {
    if (m_comp(src.min_val, dest.min_val)) {
      dest.min_val = src.min_val;
      dest.min_loc = src.min_loc;
    } else if (!m_comp(dest.min_val, src.min_val)) {
      dest.min_loc = (src.min_loc < dest.min_loc) ? src.min_loc : dest.min_loc;
    }

    if (m_comp(dest.max_val, src.max_val)) {
      dest.max_val = src.max_val;
      dest.max_loc = src.max_loc;
    } else if (!m_comp(src.max_val, dest.max_val)) {
      dest.max_loc = (src.max_loc > dest.max_loc) ? src.max_loc : dest.max_loc;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void init(value_type& val) const {
    val.max_val = ::Kokkos::reduction_identity<scalar_type>::max();
    val.min_val = ::Kokkos::reduction_identity<scalar_type>::min();
    val.max_loc = ::Kokkos::reduction_identity<index_type>::max();
    val.min_loc = ::Kokkos::reduction_identity<index_type>::min();
  }
};

template <typename Scalar, typename Index, typename ComparatorType,
          typename... Properties>
KOKKOS_DEDUCTION_GUIDE MinMaxFirstLastLocCustomComparator(
    View<MinMaxLocScalar<Scalar, Index>, Properties...> const&, ComparatorType)
    -> MinMaxFirstLastLocCustomComparator<
        Scalar, Index, ComparatorType,
        typename View<MinMaxLocScalar<Scalar, Index>,
                      Properties...>::memory_space>;

//
// FirstLoc
//
template <class Index>
struct FirstLocScalar {
  Index min_loc_true;
};

template <class Index, class Space>
struct FirstLoc
    : Impl::BaseReducer<FirstLocScalar<std::remove_cv_t<Index>>, Space> {
 private:
  using index_type = std::remove_cv_t<Index>;
  static_assert(std::is_integral_v<index_type>);
  static_assert(!std::is_pointer_v<index_type> && !std::is_array_v<index_type>);

  using parent_type =
      Impl::BaseReducer<FirstLocScalar<std::remove_cv_t<Index>>, Space>;

 public:
  using reducer    = FirstLoc<Index, Space>;
  using value_type = typename parent_type::value_type;

  // Inherit constructors
  using parent_type::parent_type;

  KOKKOS_INLINE_FUNCTION
  void join(value_type& dest, const value_type& src) const {
    dest.min_loc_true = (src.min_loc_true < dest.min_loc_true)
                            ? src.min_loc_true
                            : dest.min_loc_true;
  }

  KOKKOS_INLINE_FUNCTION
  void init(value_type& val) const {
    val.min_loc_true = ::Kokkos::reduction_identity<index_type>::min();
  }
};

template <typename Index, typename... Properties>
KOKKOS_DEDUCTION_GUIDE
FirstLoc(View<FirstLocScalar<Index>, Properties...> const&) -> FirstLoc<
    Index, typename View<FirstLocScalar<Index>, Properties...>::memory_space>;

//
// LastLoc
//
template <class Index>
struct LastLocScalar {
  Index max_loc_true;
};

template <class Index, class Space>
struct LastLoc
    : Impl::BaseReducer<LastLocScalar<std::remove_cv_t<Index>>, Space> {
 private:
  using index_type = std::remove_cv_t<Index>;
  static_assert(std::is_integral_v<index_type>);
  static_assert(!std::is_pointer_v<index_type> && !std::is_array_v<index_type>);

  using parent_type = Impl::BaseReducer<LastLocScalar<index_type>, Space>;

 public:
  using reducer    = LastLoc<Index, Space>;
  using value_type = typename parent_type::value_type;

  // Inherit constructors
  using parent_type::parent_type;

  KOKKOS_INLINE_FUNCTION
  void join(value_type& dest, const value_type& src) const {
    dest.max_loc_true = (src.max_loc_true > dest.max_loc_true)
                            ? src.max_loc_true
                            : dest.max_loc_true;
  }

  KOKKOS_INLINE_FUNCTION
  void init(value_type& val) const {
    val.max_loc_true = ::Kokkos::reduction_identity<index_type>::max();
  }
};

template <typename Index, typename... Properties>
KOKKOS_DEDUCTION_GUIDE LastLoc(View<LastLocScalar<Index>, Properties...> const&)
    -> LastLoc<Index, typename View<LastLocScalar<Index>,
                                    Properties...>::memory_space>;

template <class Index>
struct StdIsPartScalar {
  Index max_loc_true, min_loc_false;
};

//
// StdIsPartitioned
//
template <class Index, class Space>
struct StdIsPartitioned
    : Impl::BaseReducer<StdIsPartScalar<std::remove_cv_t<Index>>, Space> {
 private:
  using index_type = std::remove_cv_t<Index>;
  static_assert(std::is_integral_v<index_type>);
  static_assert(!std::is_pointer_v<index_type> && !std::is_array_v<index_type>);

  using parent_type = Impl::BaseReducer<StdIsPartScalar<index_type>, Space>;

 public:
  using reducer    = StdIsPartitioned<Index, Space>;
  using value_type = typename parent_type::value_type;

  // Inherit constructors
  using parent_type::parent_type;

  KOKKOS_INLINE_FUNCTION
  void join(value_type& dest, const value_type& src) const {
    dest.max_loc_true = (dest.max_loc_true < src.max_loc_true)
                            ? src.max_loc_true
                            : dest.max_loc_true;

    dest.min_loc_false = (dest.min_loc_false < src.min_loc_false)
                             ? dest.min_loc_false
                             : src.min_loc_false;
  }

  KOKKOS_INLINE_FUNCTION
  void init(value_type& val) const {
    val.max_loc_true  = ::Kokkos::reduction_identity<index_type>::max();
    val.min_loc_false = ::Kokkos::reduction_identity<index_type>::min();
  }
};

template <typename Index, typename... Properties>
KOKKOS_DEDUCTION_GUIDE StdIsPartitioned(
    View<StdIsPartScalar<Index>, Properties...> const&)
    -> StdIsPartitioned<Index, typename View<StdIsPartScalar<Index>,
                                             Properties...>::memory_space>;

template <class Index>
struct StdPartPointScalar {
  Index min_loc_false;
};

//
// StdPartitionPoint
//
template <class Index, class Space>
struct StdPartitionPoint
    : Impl::BaseReducer<StdPartPointScalar<std::remove_cv_t<Index>>, Space> {
 private:
  using index_type = std::remove_cv_t<Index>;
  static_assert(std::is_integral_v<index_type>);
  static_assert(!std::is_pointer_v<index_type> && !std::is_array_v<index_type>);

  using parent_type = Impl::BaseReducer<StdPartPointScalar<index_type>, Space>;

 public:
  using reducer    = StdPartitionPoint<Index, Space>;
  using value_type = typename parent_type::value_type;

  // Inherit constructors
  using parent_type::parent_type;

  KOKKOS_INLINE_FUNCTION
  void join(value_type& dest, const value_type& src) const {
    dest.min_loc_false = (dest.min_loc_false < src.min_loc_false)
                             ? dest.min_loc_false
                             : src.min_loc_false;
  }

  KOKKOS_INLINE_FUNCTION
  void init(value_type& val) const {
    val.min_loc_false = ::Kokkos::reduction_identity<index_type>::min();
  }
};

template <typename Index, typename... Properties>
KOKKOS_DEDUCTION_GUIDE StdPartitionPoint(
    View<StdPartPointScalar<Index>, Properties...> const&)
    -> StdPartitionPoint<Index, typename View<StdPartPointScalar<Index>,
                                              Properties...>::memory_space>;

}  // namespace Kokkos

#endif
