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

#include <type_traits>
#include <Kokkos_Macros.hpp>

#if !defined( KOKKOS_ROCM_INVOKE_H )
#define KOKKOS_ROCM_INVOKE_H

namespace Kokkos {
namespace Impl {

template<class Tag, class F, class... Ts, typename std::enable_if<(!std::is_void<Tag>()), int>::type = 0>
KOKKOS_INLINE_FUNCTION void rocm_invoke(F&& f, Ts&&... xs)
{
  f(Tag(), static_cast<Ts&&>(xs)...);
}

template<class Tag, class F, class... Ts, typename std::enable_if<(std::is_void<Tag>()), int>::type = 0>
KOKKOS_INLINE_FUNCTION void rocm_invoke(F&& f, Ts&&... xs)
{
  f(static_cast<Ts&&>(xs)...);
}


template<class F, class Tag=void>
struct rocm_invoke_fn
{
    F* f;
    rocm_invoke_fn(F& f_) : f(&f_)
    {}

    template<class... Ts>
    KOKKOS_INLINE_FUNCTION void operator()(Ts&&... xs) const
    {
        rocm_invoke<Tag>(*f, static_cast<Ts&&>(xs)...);
    }
};

template<class Tag, class F>
KOKKOS_INLINE_FUNCTION rocm_invoke_fn<F, Tag> make_rocm_invoke_fn(F& f)
{
    return {f};
}

template<class T>
KOKKOS_INLINE_FUNCTION T& rocm_unwrap(T& x)
{
    return x;
}

template<class T>
KOKKOS_INLINE_FUNCTION T& rocm_unwrap(std::reference_wrapper<T> x)
{
    return x;
}

template<class F, class T>
struct rocm_capture_fn
{
    F f;
    T data;

    KOKKOS_INLINE_FUNCTION rocm_capture_fn(F f_, T x) 
    : f(f_), data(x)
    {}

    template<class... Ts>
    KOKKOS_INLINE_FUNCTION void operator()(Ts&&... xs) const
    {
        f(rocm_unwrap(data), static_cast<Ts&&>(xs)...);
    }
};

template<class F, class T>
KOKKOS_INLINE_FUNCTION rocm_capture_fn<F, T> rocm_capture(F f, T x)
{
    return {f, x};
}

template<class F, class T, class U, class... Ts>
KOKKOS_INLINE_FUNCTION auto rocm_capture(F f, T x, U y, Ts... xs) -> decltype(rocm_capture(rocm_capture(f, x), y, xs...))
{
    return rocm_capture(rocm_capture(f, x), y, xs...);
}

struct rocm_apply_op
{
    template<class F, class... Ts>
    KOKKOS_INLINE_FUNCTION void operator()(F&& f, Ts&&... xs) const
    {
        f(static_cast<Ts&&>(xs)...);
    }
};

}}

#endif
