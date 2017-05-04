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


template<class T, class Enable = void>
struct is_reducer_type {
  enum { value = 0 };
};


template<class T>
struct is_reducer_type<T,typename std::enable_if<
                       std::is_same<typename std::remove_cv<T>::type,
                                    typename std::remove_cv<typename T::reducer_type>::type>::value
                      >::type> {
  enum { value = 1 };
};

namespace Experimental {


template<class Scalar,class Space = HostSpace>
struct Sum {
public:
  //Required
  typedef Sum reducer_type;
  typedef Scalar value_type;

  typedef Kokkos::View<value_type, Space, Kokkos::MemoryTraits<Kokkos::Unmanaged> > result_view_type;

  value_type init_value;

private:
  result_view_type result;

  template<class ValueType, bool is_arithmetic = std::is_arithmetic<ValueType>::value >
  struct InitWrapper;

  template<class ValueType >
  struct InitWrapper<ValueType,true> {
    static ValueType value() {
      return static_cast<value_type>(0);
    }
  };

  template<class ValueType >
  struct InitWrapper<ValueType,false> {
    static ValueType value() {
      return value_type();
    }
  };

public:

  Sum(value_type& result_):
    init_value(InitWrapper<value_type>::value()),result(&result_) {}
  Sum(const result_view_type& result_):
    init_value(InitWrapper<value_type>::value()),result(result_) {}
  Sum(value_type& result_, const value_type& init_value_):
    init_value(init_value_),result(&result_) {}
  Sum(const result_view_type& result_, const value_type& init_value_):
    init_value(init_value_),result(result_) {}

  //Required
  KOKKOS_INLINE_FUNCTION
  void join(value_type& dest, const value_type& src)  const {
    dest += src;
  }

  KOKKOS_INLINE_FUNCTION
  void join(volatile value_type& dest, const volatile value_type& src) const {
    dest += src;
  }

  //Optional
  KOKKOS_INLINE_FUNCTION
  void init( value_type& val)  const {
    val = init_value;
  }

  result_view_type result_view() const {
    return result;
  }
};

template<class Scalar,class Space = HostSpace>
struct Prod {
public:
  //Required
  typedef Prod reducer_type;
  typedef Scalar value_type;

  typedef Kokkos::View<value_type, Space, Kokkos::MemoryTraits<Kokkos::Unmanaged> > result_view_type;

  value_type init_value;

private:
  result_view_type result;

  template<class ValueType, bool is_arithmetic = std::is_arithmetic<ValueType>::value >
  struct InitWrapper;

  template<class ValueType >
  struct InitWrapper<ValueType,true> {
    static ValueType value() {
      return static_cast<value_type>(1);
    }
  };

  template<class ValueType >
  struct InitWrapper<ValueType,false> {
    static ValueType value() {
      return value_type();
    }
  };

public:

  Prod(value_type& result_):
    init_value(InitWrapper<value_type>::value()),result(&result_) {}
  Prod(const result_view_type& result_):
    init_value(InitWrapper<value_type>::value()),result(result_) {}
  Prod(value_type& result_, const value_type& init_value_):
    init_value(init_value_),result(&result_) {}
  Prod(const result_view_type& result_, const value_type& init_value_):
    init_value(init_value_),result(result_) {}

  //Required
  KOKKOS_INLINE_FUNCTION
  void join(value_type& dest, const value_type& src)  const {
    dest *= src;
  }

  KOKKOS_INLINE_FUNCTION
  void join(volatile value_type& dest, const volatile value_type& src) const {
    dest *= src;
  }

  //Optional
  KOKKOS_INLINE_FUNCTION
  void init( value_type& val)  const {
    val = init_value;
  }

  result_view_type result_view() const {
    return result;
  }
};

template<class Scalar, class Space = HostSpace>
struct Min {
public:
  //Required
  typedef Min reducer_type;
  typedef typename std::remove_cv<Scalar>::type value_type;

  typedef Kokkos::View<value_type, Space, Kokkos::MemoryTraits<Kokkos::Unmanaged> > result_view_type;

  value_type init_value;

private:
  result_view_type result;

  template<class ValueType, bool is_arithmetic = std::is_arithmetic<ValueType>::value >
  struct InitWrapper;

  template<class ValueType >
  struct InitWrapper<ValueType,true> {
    static ValueType value() {
      return std::numeric_limits<value_type>::max();
    }
  };

  template<class ValueType >
  struct InitWrapper<ValueType,false> {
    static ValueType value() {
      return value_type();
    }
  };

public:

  Min(value_type& result_):
    init_value(InitWrapper<value_type>::value()),result(&result_) {}
  Min(const result_view_type& result_):
    init_value(InitWrapper<value_type>::value()),result(result_) {}
  Min(value_type& result_, const value_type& init_value_):
    init_value(init_value_),result(&result_) {}
  Min(const result_view_type& result_, const value_type& init_value_):
    init_value(init_value_),result(result_) {}

  //Required
  KOKKOS_INLINE_FUNCTION
  void join(value_type& dest, const value_type& src)  const {
    if ( src < dest )
      dest = src;
  }

  KOKKOS_INLINE_FUNCTION
  void join(volatile value_type& dest, const volatile value_type& src) const {
    if ( src < dest )
      dest = src;
  }

  //Optional
  KOKKOS_INLINE_FUNCTION
  void init( value_type& val)  const {
    val = init_value;
  }

  result_view_type result_view() const {
    return result;
  }
};

template<class Scalar, class Space = HostSpace>
struct Max {
public:
  //Required
  typedef Max reducer_type;
  typedef typename std::remove_cv<Scalar>::type value_type;

  typedef Kokkos::View<value_type, Space, Kokkos::MemoryTraits<Kokkos::Unmanaged> > result_view_type;

  value_type init_value;

private:
  result_view_type result;

  template<class ValueType, bool is_arithmetic = std::is_arithmetic<ValueType>::value >
  struct InitWrapper;

  template<class ValueType >
  struct InitWrapper<ValueType,true> {
    static ValueType value() {
      return std::numeric_limits<value_type>::min();
    }
  };

  template<class ValueType >
  struct InitWrapper<ValueType,false> {
    static ValueType value() {
      return value_type();
    }
  };

public:

  Max(value_type& result_):
    init_value(InitWrapper<value_type>::value()),result(&result_) {}
  Max(const result_view_type& result_):
    init_value(InitWrapper<value_type>::value()),result(result_) {}
  Max(value_type& result_, const value_type& init_value_):
    init_value(init_value_),result(&result_) {}
  Max(const result_view_type& result_, const value_type& init_value_):
    init_value(init_value_),result(result_) {}

  //Required
  KOKKOS_INLINE_FUNCTION
  void join(value_type& dest, const value_type& src)  const {
    if ( src > dest )
      dest = src;
  }

  KOKKOS_INLINE_FUNCTION
  void join(volatile value_type& dest, const volatile value_type& src) const {
    if ( src > dest )
      dest = src;
  }

  //Optional
  KOKKOS_INLINE_FUNCTION
  void init( value_type& val)  const {
    val = init_value;
  }

  result_view_type result_view() const {
    return result;
  }
};

template<class Scalar, class Space = HostSpace>
struct LAnd {
public:
  //Required
  typedef LAnd reducer_type;
  typedef Scalar value_type;

  typedef Kokkos::View<value_type, Space, Kokkos::MemoryTraits<Kokkos::Unmanaged> > result_view_type;

private:
  result_view_type result;

public:

  LAnd(value_type& result_):result(&result_) {}
  LAnd(const result_view_type& result_):result(result_) {}

  //Required
  KOKKOS_INLINE_FUNCTION
  void join(value_type& dest, const value_type& src)  const {
    dest = dest && src;
  }

  KOKKOS_INLINE_FUNCTION
  void join(volatile value_type& dest, const volatile value_type& src) const {
    dest = dest && src;
  }

  //Optional
  KOKKOS_INLINE_FUNCTION
  void init( value_type& val)  const {
    val = 1;
  }

  result_view_type result_view() const {
    return result;
  }
};

template<class Scalar, class Space = HostSpace>
struct LOr {
public:
  //Required
  typedef LOr reducer_type;
  typedef Scalar value_type;

  typedef Kokkos::View<value_type, Space, Kokkos::MemoryTraits<Kokkos::Unmanaged> > result_view_type;

private:
  result_view_type result;

public:

  LOr(value_type& result_):result(&result_) {}
  LOr(const result_view_type& result_):result(result_) {}

  //Required
  KOKKOS_INLINE_FUNCTION
  void join(value_type& dest, const value_type& src)  const {
    dest = dest || src;
  }

  KOKKOS_INLINE_FUNCTION
  void join(volatile value_type& dest, const volatile value_type& src) const {
    dest = dest || src;
  }

  //Optional
  KOKKOS_INLINE_FUNCTION
  void init( value_type& val)  const {
    val = 0;
  }

  result_view_type result_view() const {
    return result;
  }
};

template<class Scalar, class Space = HostSpace>
struct LXor {
public:
  //Required
  typedef LXor reducer_type;
  typedef Scalar value_type;

  typedef Kokkos::View<value_type, Space, Kokkos::MemoryTraits<Kokkos::Unmanaged> > result_view_type;

private:
  result_view_type result;

public:

  LXor(value_type& result_):result(&result_) {}
  LXor(const result_view_type& result_):result(result_) {}

  //Required
  KOKKOS_INLINE_FUNCTION
  void join(value_type& dest, const value_type& src)  const {
    dest = dest? (!src) : src;
  }

  KOKKOS_INLINE_FUNCTION
  void join(volatile value_type& dest, const volatile value_type& src) const {
    dest = dest? (!src) : src;
  }

  //Optional
  KOKKOS_INLINE_FUNCTION
  void init( value_type& val)  const {
    val = 0;
  }

  result_view_type result_view() const {
    return result;
  }
};

template<class Scalar, class Space = HostSpace>
struct BAnd {
public:
  //Required
  typedef BAnd reducer_type;
  typedef typename std::remove_cv<Scalar>::type value_type;

  typedef Kokkos::View<value_type, Space, Kokkos::MemoryTraits<Kokkos::Unmanaged> > result_view_type;

  value_type init_value;

private:
  result_view_type result;

public:

  BAnd(value_type& result_):
    init_value(value_type() | (~value_type())),result(&result_) {}
  BAnd(const result_view_type& result_):
    init_value(value_type() | (~value_type())),result(result_) {}

  //Required
  KOKKOS_INLINE_FUNCTION
  void join(value_type& dest, const value_type& src)  const {
      dest = dest & src;
  }

  KOKKOS_INLINE_FUNCTION
  void join(volatile value_type& dest, const volatile value_type& src) const {
    dest = dest & src;
  }

  //Optional
  KOKKOS_INLINE_FUNCTION
  void init( value_type& val)  const {
    val = init_value;
  }

  result_view_type result_view() const {
    return result;
  }
};

template<class Scalar, class Space = HostSpace>
struct BOr {
public:
  //Required
  typedef BOr reducer_type;
  typedef typename std::remove_cv<Scalar>::type value_type;

  typedef Kokkos::View<value_type, Space, Kokkos::MemoryTraits<Kokkos::Unmanaged> > result_view_type;

  value_type init_value;

private:
  result_view_type result;

public:

  BOr(value_type& result_):
    init_value(value_type() & (~value_type())),result(&result_) {}
  BOr(const result_view_type& result_):
    init_value(value_type() & (~value_type())),result(result_) {}

  //Required
  KOKKOS_INLINE_FUNCTION
  void join(value_type& dest, const value_type& src)  const {
      dest = dest | src;
  }

  KOKKOS_INLINE_FUNCTION
  void join(volatile value_type& dest, const volatile value_type& src) const {
    dest = dest | src;
  }

  //Optional
  KOKKOS_INLINE_FUNCTION
  void init( value_type& val)  const {
    val = init_value;
  }

  result_view_type result_view() const {
    return result;
  }
};

template<class Scalar, class Space = HostSpace>
struct BXor {
public:
  //Required
  typedef BXor reducer_type;
  typedef typename std::remove_cv<Scalar>::type value_type;

  typedef Kokkos::View<value_type, Space, Kokkos::MemoryTraits<Kokkos::Unmanaged> > result_view_type;

  value_type init_value;

private:
  result_view_type result;

public:

  BXor(value_type& result_):
    init_value(value_type() & (~value_type())),result(&result_) {}
  BXor(const result_view_type& result_):
    init_value(value_type() & (~value_type())),result(result_) {}

  //Required
  KOKKOS_INLINE_FUNCTION
  void join(value_type& dest, const value_type& src)  const {
      dest = dest ^ src;
  }

  KOKKOS_INLINE_FUNCTION
  void join(volatile value_type& dest, const volatile value_type& src) const {
    dest = dest ^ src;
  }

  //Optional
  KOKKOS_INLINE_FUNCTION
  void init( value_type& val)  const {
    val = init_value;
  }

  result_view_type result_view() const {
    return result;
  }
};

template<class Scalar, class Index>
struct ValLocScalar {
  Scalar val;
  Index loc;

  KOKKOS_INLINE_FUNCTION
  void operator = (const ValLocScalar& rhs) {
    val = rhs.val;
    loc = rhs.loc;
  }

  KOKKOS_INLINE_FUNCTION
  void operator = (const volatile ValLocScalar& rhs) volatile {
    val = rhs.val;
    loc = rhs.loc;
  }
};

template<class Scalar, class Index, class Space = HostSpace>
struct MinLoc {
private:
  typedef typename std::remove_cv<Scalar>::type scalar_type;
  typedef typename std::remove_cv<Index>::type index_type;

public:
  //Required
  typedef MinLoc reducer_type;
  typedef ValLocScalar<scalar_type,index_type> value_type;

  typedef Kokkos::View<value_type, Space, Kokkos::MemoryTraits<Kokkos::Unmanaged> > result_view_type;

  scalar_type init_value;

private:
  result_view_type result;

  template<class ValueType, bool is_arithmetic = std::is_arithmetic<ValueType>::value >
  struct InitWrapper;

  template<class ValueType >
  struct InitWrapper<ValueType,true> {
    static ValueType value() {
      return std::numeric_limits<scalar_type>::max();
    }
  };

  template<class ValueType >
  struct InitWrapper<ValueType,false> {
    static ValueType value() {
      return scalar_type();
    }
  };

public:

  MinLoc(value_type& result_):
    init_value(InitWrapper<scalar_type>::value()),result(&result_) {}
  MinLoc(const result_view_type& result_):
    init_value(InitWrapper<scalar_type>::value()),result(result_) {}
  MinLoc(value_type& result_, const scalar_type& init_value_):
    init_value(init_value_),result(&result_) {}
  MinLoc(const result_view_type& result_, const scalar_type& init_value_):
    init_value(init_value_),result(result_) {}


  //Required
  KOKKOS_INLINE_FUNCTION
  void join(value_type& dest, const value_type& src)  const {
    if ( src.val < dest.val )
      dest = src;
  }

  KOKKOS_INLINE_FUNCTION
  void join(volatile value_type& dest, const volatile value_type& src) const {
    if ( src.val < dest.val )
      dest = src;
  }

  //Optional
  KOKKOS_INLINE_FUNCTION
  void init( value_type& val)  const {
    val.val = init_value;
  }

  result_view_type result_view() const {
    return result;
  }
};

template<class Scalar, class Index, class Space = HostSpace>
struct MaxLoc {
private:
  typedef typename std::remove_cv<Scalar>::type scalar_type;
  typedef typename std::remove_cv<Index>::type index_type;

public:
  //Required
  typedef MaxLoc reducer_type;
  typedef ValLocScalar<scalar_type,index_type> value_type;

  typedef Kokkos::View<value_type, Space, Kokkos::MemoryTraits<Kokkos::Unmanaged> > result_view_type;

  scalar_type init_value;

private:
  result_view_type result;

  template<class ValueType, bool is_arithmetic = std::is_arithmetic<ValueType>::value >
  struct InitWrapper;

  template<class ValueType >
  struct InitWrapper<ValueType,true> {
    static ValueType value() {
      return std::numeric_limits<scalar_type>::min();
    }
  };

  template<class ValueType >
  struct InitWrapper<ValueType,false> {
    static ValueType value() {
      return scalar_type();
    }
  };

public:

  MaxLoc(value_type& result_):
    init_value(InitWrapper<scalar_type>::value()),result(&result_) {}
  MaxLoc(const result_view_type& result_):
    init_value(InitWrapper<scalar_type>::value()),result(result_) {}
  MaxLoc(value_type& result_, const scalar_type& init_value_):
    init_value(init_value_),result(&result_) {}
  MaxLoc(const result_view_type& result_, const scalar_type& init_value_):
    init_value(init_value_),result(result_) {}

  //Required
  KOKKOS_INLINE_FUNCTION
  void join(value_type& dest, const value_type& src)  const {
    if ( src.val > dest.val )
      dest = src;
  }

  KOKKOS_INLINE_FUNCTION
  void join(volatile value_type& dest, const volatile value_type& src) const {
    if ( src.val > dest.val )
      dest = src;
  }

  //Optional
  KOKKOS_INLINE_FUNCTION
  void init( value_type& val)  const {
    val.val = init_value;
  }

  result_view_type result_view() const {
    return result;
  }
};

template<class Scalar>
struct MinMaxScalar {
  Scalar min_val,max_val;

  KOKKOS_INLINE_FUNCTION
  void operator = (const MinMaxScalar& rhs) {
    min_val = rhs.min_val;
    max_val = rhs.max_val;
  }

  KOKKOS_INLINE_FUNCTION
  void operator = (const volatile MinMaxScalar& rhs) volatile {
    min_val = rhs.min_val;
    max_val = rhs.max_val;
  }
};

template<class Scalar, class Space = HostSpace>
struct MinMax {
private:
  typedef typename std::remove_cv<Scalar>::type scalar_type;

public:
  //Required
  typedef MinMax reducer_type;
  typedef MinMaxScalar<scalar_type> value_type;

  typedef Kokkos::View<value_type, Space, Kokkos::MemoryTraits<Kokkos::Unmanaged> > result_view_type;

  scalar_type min_init_value;
  scalar_type max_init_value;

private:
  result_view_type result;

  template<class ValueType, bool is_arithmetic = std::is_arithmetic<ValueType>::value >
  struct MinInitWrapper;

  template<class ValueType >
  struct MinInitWrapper<ValueType,true> {
    static ValueType value() {
      return std::numeric_limits<scalar_type>::max();
    }
  };

  template<class ValueType >
  struct MinInitWrapper<ValueType,false> {
    static ValueType value() {
      return scalar_type();
    }
  };

  template<class ValueType, bool is_arithmetic = std::is_arithmetic<ValueType>::value >
  struct MaxInitWrapper;

  template<class ValueType >
  struct MaxInitWrapper<ValueType,true> {
    static ValueType value() {
      return std::numeric_limits<scalar_type>::min();
    }
  };

  template<class ValueType >
  struct MaxInitWrapper<ValueType,false> {
    static ValueType value() {
      return scalar_type();
    }
  };

public:

  MinMax(value_type& result_):
    min_init_value(MinInitWrapper<scalar_type>::value()),max_init_value(MaxInitWrapper<scalar_type>::value()),result(&result_) {}
  MinMax(const result_view_type& result_):
    min_init_value(MinInitWrapper<scalar_type>::value()),max_init_value(MaxInitWrapper<scalar_type>::value()),result(result_) {}
  MinMax(value_type& result_, const scalar_type& min_init_value_, const scalar_type& max_init_value_):
    min_init_value(min_init_value_),max_init_value(max_init_value_),result(&result_) {}
  MinMax(const result_view_type& result_, const scalar_type& min_init_value_, const scalar_type& max_init_value_):
    min_init_value(min_init_value_),max_init_value(max_init_value_),result(result_) {}

  //Required
  KOKKOS_INLINE_FUNCTION
  void join(value_type& dest, const value_type& src)  const {
    if ( src.min_val < dest.min_val ) {
      dest.min_val = src.min_val;
    }
    if ( src.max_val > dest.max_val ) {
      dest.max_val = src.max_val;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void join(volatile value_type& dest, const volatile value_type& src) const {
    if ( src.min_val < dest.min_val ) {
      dest.min_val = src.min_val;
    }
    if ( src.max_val > dest.max_val ) {
      dest.max_val = src.max_val;
    }
  }

  //Optional
  KOKKOS_INLINE_FUNCTION
  void init( value_type& val)  const {
    val.min_val = min_init_value;
    val.max_val = max_init_value;
  }

  result_view_type result_view() const {
    return result;
  }
};

template<class Scalar, class Index>
struct MinMaxLocScalar {
  Scalar min_val,max_val;
  Index min_loc,max_loc;

  KOKKOS_INLINE_FUNCTION
  void operator = (const MinMaxLocScalar& rhs) {
    min_val = rhs.min_val;
    min_loc = rhs.min_loc;
    max_val = rhs.max_val;
    max_loc = rhs.max_loc;
  }

  KOKKOS_INLINE_FUNCTION
  void operator = (const volatile MinMaxLocScalar& rhs) volatile {
    min_val = rhs.min_val;
    min_loc = rhs.min_loc;
    max_val = rhs.max_val;
    max_loc = rhs.max_loc;
  }
};

template<class Scalar, class Index, class Space = HostSpace>
struct MinMaxLoc {
private:
  typedef typename std::remove_cv<Scalar>::type scalar_type;
  typedef typename std::remove_cv<Index>::type index_type;

public:
  //Required
  typedef MinMaxLoc reducer_type;
  typedef MinMaxLocScalar<scalar_type,index_type> value_type;

  typedef Kokkos::View<value_type, Space, Kokkos::MemoryTraits<Kokkos::Unmanaged> > result_view_type;

  scalar_type min_init_value;
  scalar_type max_init_value;

private:
  result_view_type result;

  template<class ValueType, bool is_arithmetic = std::is_arithmetic<ValueType>::value >
  struct MinInitWrapper;

  template<class ValueType >
  struct MinInitWrapper<ValueType,true> {
    static ValueType value() {
      return std::numeric_limits<scalar_type>::max();
    }
  };

  template<class ValueType >
  struct MinInitWrapper<ValueType,false> {
    static ValueType value() {
      return scalar_type();
    }
  };

  template<class ValueType, bool is_arithmetic = std::is_arithmetic<ValueType>::value >
  struct MaxInitWrapper;

  template<class ValueType >
  struct MaxInitWrapper<ValueType,true> {
    static ValueType value() {
      return std::numeric_limits<scalar_type>::min();
    }
  };

  template<class ValueType >
  struct MaxInitWrapper<ValueType,false> {
    static ValueType value() {
      return scalar_type();
    }
  };

public:

  MinMaxLoc(value_type& result_):
    min_init_value(MinInitWrapper<scalar_type>::value()),max_init_value(MaxInitWrapper<scalar_type>::value()),result(&result_) {}
  MinMaxLoc(const result_view_type& result_):
    min_init_value(MinInitWrapper<scalar_type>::value()),max_init_value(MaxInitWrapper<scalar_type>::value()),result(result_) {}
  MinMaxLoc(value_type& result_, const scalar_type& min_init_value_, const scalar_type& max_init_value_):
    min_init_value(min_init_value_),max_init_value(max_init_value_),result(&result_) {}
  MinMaxLoc(const result_view_type& result_, const scalar_type& min_init_value_, const scalar_type& max_init_value_):
    min_init_value(min_init_value_),max_init_value(max_init_value_),result(result_) {}

  //Required
  KOKKOS_INLINE_FUNCTION
  void join(value_type& dest, const value_type& src)  const {
    if ( src.min_val < dest.min_val ) {
      dest.min_val = src.min_val;
      dest.min_loc = src.min_loc;
    }
    if ( src.max_val > dest.max_val ) {
      dest.max_val = src.max_val;
      dest.max_loc = src.max_loc;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void join(volatile value_type& dest, const volatile value_type& src) const {
    if ( src.min_val < dest.min_val ) {
      dest.min_val = src.min_val;
      dest.min_loc = src.min_loc;
    }
    if ( src.max_val > dest.max_val ) {
      dest.max_val = src.max_val;
      dest.max_loc = src.max_loc;
    }
  }

  //Optional
  KOKKOS_INLINE_FUNCTION
  void init( value_type& val)  const {
    val.min_val = min_init_value;
    val.max_val = max_init_value;
  }

  result_view_type result_view() const {
    return result;
  }
};
}
}


namespace Kokkos {
namespace Impl {

template< class T, class ReturnType , class ValueTraits>
struct ParallelReduceReturnValue;

template< class ReturnType , class FunctorType >
struct ParallelReduceReturnValue<typename std::enable_if<Kokkos::is_view<ReturnType>::value>::type, ReturnType, FunctorType> {
  typedef ReturnType return_type;
  typedef InvalidType reducer_type;

  typedef typename return_type::value_type value_type_scalar;
  typedef typename return_type::value_type* const value_type_array;

  typedef typename if_c<return_type::rank==0,value_type_scalar,value_type_array>::type value_type;

  static return_type& return_value(ReturnType& return_val, const FunctorType&) {
    return return_val;
  }
};

template< class ReturnType , class FunctorType>
struct ParallelReduceReturnValue<typename std::enable_if<
                                   !Kokkos::is_view<ReturnType>::value &&
                                  (!std::is_array<ReturnType>::value && !std::is_pointer<ReturnType>::value) &&
                                   !Kokkos::is_reducer_type<ReturnType>::value
                                 >::type, ReturnType, FunctorType> {
  typedef Kokkos::View<  ReturnType
                       , Kokkos::HostSpace
                       , Kokkos::MemoryUnmanaged
      > return_type;

  typedef InvalidType reducer_type;

  typedef typename return_type::value_type value_type;

  static return_type return_value(ReturnType& return_val, const FunctorType&) {
    return return_type(&return_val);
  }
};

template< class ReturnType , class FunctorType>
struct ParallelReduceReturnValue<typename std::enable_if<
                                  (is_array<ReturnType>::value || std::is_pointer<ReturnType>::value)
                                >::type, ReturnType, FunctorType> {
  typedef Kokkos::View<  typename std::remove_const<ReturnType>::type
                       , Kokkos::HostSpace
                       , Kokkos::MemoryUnmanaged
      > return_type;

  typedef InvalidType reducer_type;

  typedef typename return_type::value_type value_type[];

  static return_type return_value(ReturnType& return_val,
                                  const FunctorType& functor) {
    return return_type(return_val,functor.value_count);
  }
};

template< class ReturnType , class FunctorType>
struct ParallelReduceReturnValue<typename std::enable_if<
                                   Kokkos::is_reducer_type<ReturnType>::value
                                >::type, ReturnType, FunctorType> {
  typedef ReturnType return_type;
  typedef ReturnType reducer_type;
  typedef typename return_type::value_type value_type;

  static return_type return_value(ReturnType& return_val,
                                  const FunctorType& functor) {
    return return_val;
  }
};
}

namespace Impl {
template< class T, class ReturnType , class FunctorType>
struct ParallelReducePolicyType;

template< class PolicyType , class FunctorType >
struct ParallelReducePolicyType<typename std::enable_if<Kokkos::Impl::is_execution_policy<PolicyType>::value>::type, PolicyType,FunctorType> {

  typedef PolicyType policy_type;
  static PolicyType policy(const PolicyType& policy_) {
    return policy_;
  }
};

template< class PolicyType , class FunctorType >
struct ParallelReducePolicyType<typename std::enable_if<std::is_integral<PolicyType>::value>::type, PolicyType,FunctorType> {
  typedef typename
    Impl::FunctorPolicyExecutionSpace< FunctorType , void >::execution_space
      execution_space ;

  typedef Kokkos::RangePolicy<execution_space> policy_type;

  static policy_type policy(const PolicyType& policy_) {
    return policy_type(0,policy_);
  }
};

}

namespace Impl {
  template< class FunctorType, class ExecPolicy, class ValueType, class ExecutionSpace>
  struct ParallelReduceFunctorType {
    typedef FunctorType functor_type;
    static const functor_type& functor(const functor_type& functor) {
      return functor;
    }
  };
}

namespace Impl {

  template< class PolicyType, class FunctorType, class ReturnType >
  struct ParallelReduceAdaptor {
    typedef Impl::ParallelReduceReturnValue<void,ReturnType,FunctorType> return_value_adapter;
    #ifdef KOKKOS_IMPL_NEED_FUNCTOR_WRAPPER
    typedef Impl::ParallelReduceFunctorType<FunctorType,PolicyType,
                                            typename return_value_adapter::value_type,
                                            typename PolicyType::execution_space> functor_adaptor;
    #endif
    static inline
    void execute(const std::string& label,
        const PolicyType& policy,
        const FunctorType& functor,
        ReturnType& return_value) {
          #if defined(KOKKOS_ENABLE_PROFILING)
            uint64_t kpID = 0;
            if(Kokkos::Profiling::profileLibraryLoaded()) {
              Kokkos::Profiling::beginParallelReduce("" == label ? typeid(FunctorType).name() : label, 0, &kpID);
            }
          #endif

          Kokkos::Impl::shared_allocation_tracking_claim_and_disable();
          #ifdef KOKKOS_IMPL_NEED_FUNCTOR_WRAPPER
          Impl::ParallelReduce<typename functor_adaptor::functor_type, PolicyType, typename return_value_adapter::reducer_type >
             closure(functor_adaptor::functor(functor),
                     policy,
                     return_value_adapter::return_value(return_value,functor));
          #else
          Impl::ParallelReduce<FunctorType, PolicyType, typename return_value_adapter::reducer_type >
             closure(functor,
                     policy,
                     return_value_adapter::return_value(return_value,functor));
          #endif
          Kokkos::Impl::shared_allocation_tracking_release_and_enable();
          closure.execute();

          #if defined(KOKKOS_ENABLE_PROFILING)
            if(Kokkos::Profiling::profileLibraryLoaded()) {
              Kokkos::Profiling::endParallelReduce(kpID);
            }
          #endif
        }

  };
}
/*! \fn void parallel_reduce(label,policy,functor,return_argument)
    \brief Perform a parallel reduction.
    \param label An optional Label giving the call name. Must be able to construct a std::string from the argument.
    \param policy A Kokkos Execution Policy, such as an integer, a RangePolicy or a TeamPolicy.
    \param functor A functor with a reduction operator, and optional init, join and final functions.
    \param return_argument A return argument which can be a scalar, a View, or a ReducerStruct. This argument can be left out if the functor has a final function.
*/

/** \brief  Parallel reduction
 *
 * parallel_reduce performs parallel reductions with arbitrary functions - i.e.
 * it is not solely data based. The call expects up to 4 arguments:
 *
 *
 * Example of a parallel_reduce functor for a POD (plain old data) value type:
 * \code
 *  class FunctorType { // For POD value type
 *  public:
 *    typedef    ...     execution_space ;
 *    typedef <podType>  value_type ;
 *    void operator()( <intType> iwork , <podType> & update ) const ;
 *    void init( <podType> & update ) const ;
 *    void join( volatile       <podType> & update ,
 *               volatile const <podType> & input ) const ;
 *
 *    typedef true_type has_final ;
 *    void final( <podType> & update ) const ;
 *  };
 * \endcode
 *
 * Example of a parallel_reduce functor for an array of POD (plain old data) values:
 * \code
 *  class FunctorType { // For array of POD value
 *  public:
 *    typedef    ...     execution_space ;
 *    typedef <podType>  value_type[] ;
 *    void operator()( <intType> , <podType> update[] ) const ;
 *    void init( <podType> update[] ) const ;
 *    void join( volatile       <podType> update[] ,
 *               volatile const <podType> input[] ) const ;
 *
 *    typedef true_type has_final ;
 *    void final( <podType> update[] ) const ;
 *  };
 * \endcode
 */

// ReturnValue is scalar or array: take by reference

template< class PolicyType, class FunctorType, class ReturnType >
inline
void parallel_reduce(const std::string& label,
                     const PolicyType& policy,
                     const FunctorType& functor,
                     ReturnType& return_value,
                     typename Impl::enable_if<
                       Kokkos::Impl::is_execution_policy<PolicyType>::value
                     >::type * = 0) {
  Impl::ParallelReduceAdaptor<PolicyType,FunctorType,ReturnType>::execute(label,policy,functor,return_value);
}

template< class PolicyType, class FunctorType, class ReturnType >
inline
void parallel_reduce(const PolicyType& policy,
                     const FunctorType& functor,
                     ReturnType& return_value,
                     typename Impl::enable_if<
                       Kokkos::Impl::is_execution_policy<PolicyType>::value
                     >::type * = 0) {
  Impl::ParallelReduceAdaptor<PolicyType,FunctorType,ReturnType>::execute("",policy,functor,return_value);
}

template< class FunctorType, class ReturnType >
inline
void parallel_reduce(const size_t& policy,
                     const FunctorType& functor,
                     ReturnType& return_value) {
  typedef typename Impl::ParallelReducePolicyType<void,size_t,FunctorType>::policy_type policy_type;
  Impl::ParallelReduceAdaptor<policy_type,FunctorType,ReturnType>::execute("",policy_type(0,policy),functor,return_value);
}

template< class FunctorType, class ReturnType >
inline
void parallel_reduce(const std::string& label,
                     const size_t& policy,
                     const FunctorType& functor,
                     ReturnType& return_value) {
  typedef typename Impl::ParallelReducePolicyType<void,size_t,FunctorType>::policy_type policy_type;
  Impl::ParallelReduceAdaptor<policy_type,FunctorType,ReturnType>::execute(label,policy_type(0,policy),functor,return_value);
}

// ReturnValue as View or Reducer: take by copy to allow for inline construction

template< class PolicyType, class FunctorType, class ReturnType >
inline
void parallel_reduce(const std::string& label,
                     const PolicyType& policy,
                     const FunctorType& functor,
                     const ReturnType& return_value,
                     typename Impl::enable_if<
                       Kokkos::Impl::is_execution_policy<PolicyType>::value
                     >::type * = 0) {
  Impl::ParallelReduceAdaptor<PolicyType,FunctorType,const ReturnType>::execute(label,policy,functor,return_value);
}

template< class PolicyType, class FunctorType, class ReturnType >
inline
void parallel_reduce(const PolicyType& policy,
                     const FunctorType& functor,
                     const ReturnType& return_value,
                     typename Impl::enable_if<
                       Kokkos::Impl::is_execution_policy<PolicyType>::value
                     >::type * = 0) {
  ReturnType return_value_impl = return_value;
  Impl::ParallelReduceAdaptor<PolicyType,FunctorType,ReturnType>::execute("",policy,functor,return_value_impl);
}

template< class FunctorType, class ReturnType >
inline
void parallel_reduce(const size_t& policy,
                     const FunctorType& functor,
                     const ReturnType& return_value) {
  typedef typename Impl::ParallelReducePolicyType<void,size_t,FunctorType>::policy_type policy_type;
  ReturnType return_value_impl = return_value;
  Impl::ParallelReduceAdaptor<policy_type,FunctorType,ReturnType>::execute("",policy_type(0,policy),functor,return_value_impl);
}

template< class FunctorType, class ReturnType >
inline
void parallel_reduce(const std::string& label,
                     const size_t& policy,
                     const FunctorType& functor,
                     const ReturnType& return_value) {
  typedef typename Impl::ParallelReducePolicyType<void,size_t,FunctorType>::policy_type policy_type;
  ReturnType return_value_impl = return_value;
  Impl::ParallelReduceAdaptor<policy_type,FunctorType,ReturnType>::execute(label,policy_type(0,policy),functor,return_value_impl);
}

// No Return Argument

template< class PolicyType, class FunctorType>
inline
void parallel_reduce(const std::string& label,
                     const PolicyType& policy,
                     const FunctorType& functor,
                     typename Impl::enable_if<
                       Kokkos::Impl::is_execution_policy<PolicyType>::value
                     >::type * = 0) {
  typedef Kokkos::Impl::FunctorValueTraits< FunctorType , void >  ValueTraits ;
  typedef typename Kokkos::Impl::if_c< (ValueTraits::StaticValueSize != 0)
                                     , typename ValueTraits::value_type
                                     , typename ValueTraits::pointer_type
                                     >::type value_type ;

  typedef Kokkos::View< value_type
              , Kokkos::HostSpace
              , Kokkos::MemoryUnmanaged
              > result_view_type;
  result_view_type result_view ;

  Impl::ParallelReduceAdaptor<PolicyType,FunctorType,result_view_type>::execute(label,policy,functor,result_view);
}

template< class PolicyType, class FunctorType >
inline
void parallel_reduce(const PolicyType& policy,
                     const FunctorType& functor,
                     typename Impl::enable_if<
                       Kokkos::Impl::is_execution_policy<PolicyType>::value
                     >::type * = 0) {
  typedef Kokkos::Impl::FunctorValueTraits< FunctorType , void >  ValueTraits ;
  typedef typename Kokkos::Impl::if_c< (ValueTraits::StaticValueSize != 0)
                                     , typename ValueTraits::value_type
                                     , typename ValueTraits::pointer_type
                                     >::type value_type ;

  typedef Kokkos::View< value_type
              , Kokkos::HostSpace
              , Kokkos::MemoryUnmanaged
              > result_view_type;
  result_view_type result_view ;

  Impl::ParallelReduceAdaptor<PolicyType,FunctorType,result_view_type>::execute("",policy,functor,result_view);
}

template< class FunctorType >
inline
void parallel_reduce(const size_t& policy,
                     const FunctorType& functor) {
  typedef typename Impl::ParallelReducePolicyType<void,size_t,FunctorType>::policy_type policy_type;
  typedef Kokkos::Impl::FunctorValueTraits< FunctorType , void >  ValueTraits ;
  typedef typename Kokkos::Impl::if_c< (ValueTraits::StaticValueSize != 0)
                                     , typename ValueTraits::value_type
                                     , typename ValueTraits::pointer_type
                                     >::type value_type ;

  typedef Kokkos::View< value_type
              , Kokkos::HostSpace
              , Kokkos::MemoryUnmanaged
              > result_view_type;
  result_view_type result_view ;

  Impl::ParallelReduceAdaptor<policy_type,FunctorType,result_view_type>::execute("",policy_type(0,policy),functor,result_view);
}

template< class FunctorType>
inline
void parallel_reduce(const std::string& label,
                     const size_t& policy,
                     const FunctorType& functor) {
  typedef typename Impl::ParallelReducePolicyType<void,size_t,FunctorType>::policy_type policy_type;
  typedef Kokkos::Impl::FunctorValueTraits< FunctorType , void >  ValueTraits ;
  typedef typename Kokkos::Impl::if_c< (ValueTraits::StaticValueSize != 0)
                                     , typename ValueTraits::value_type
                                     , typename ValueTraits::pointer_type
                                     >::type value_type ;

  typedef Kokkos::View< value_type
              , Kokkos::HostSpace
              , Kokkos::MemoryUnmanaged
              > result_view_type;
  result_view_type result_view ;

  Impl::ParallelReduceAdaptor<policy_type,FunctorType,result_view_type>::execute(label,policy_type(0,policy),functor,result_view);
}



} //namespace Kokkos
