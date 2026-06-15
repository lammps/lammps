// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <Kokkos_TypeInfo.hpp>

#include <type_traits>

namespace {

using Kokkos::Impl::TypeInfo;

struct Foo {};
using FooAlias = Foo;
enum Bar { BAR_0, BAR_1, BAR_2 };
union Baz {
  int i;
  float f;
};

[[maybe_unused]] auto func = [](int) {};  // < line 20
//                           ^  column 30
using Lambda = decltype(func);

// clang-format off
#if defined(__NVCC__) && !defined(__CUDA_ARCH__)
// can't do much
// it looks like that there is 1st an EDG pass and then a host pass and they cannot both agree on what the type info is
#elif defined(__EDG__) || (defined(__NVCC__) && defined(__CUDA_ARCH__))
static_assert(TypeInfo<Foo>::name()      == "<unnamed>::Foo");
static_assert(TypeInfo<FooAlias>::name() == "<unnamed>::Foo");
static_assert(TypeInfo<Bar>::name()      == "<unnamed>::Bar");
static_assert(TypeInfo<Baz>::name()      == "<unnamed>::Baz");
static_assert(TypeInfo<Lambda>::name()   == "lambda [](int)->void");
#elif defined(__clang__)
static_assert(TypeInfo<Foo>::name()      == "(anonymous namespace)::Foo");
static_assert(TypeInfo<FooAlias>::name() == "(anonymous namespace)::Foo");
static_assert(TypeInfo<Bar>::name()      == "(anonymous namespace)::Bar");
static_assert(TypeInfo<Baz>::name()      == "(anonymous namespace)::Baz");
static_assert(TypeInfo<Lambda>::name()   == "(anonymous namespace)::(lambda at "  __FILE__  ":20:30)");
#elif defined(__GNUC__)
static_assert(TypeInfo<Foo>::name()      == "{anonymous}::Foo");
static_assert(TypeInfo<FooAlias>::name() == "{anonymous}::Foo");
static_assert(TypeInfo<Bar>::name()      == "{anonymous}::Bar");
static_assert(TypeInfo<Baz>::name()      == "{anonymous}::Baz");
static_assert(TypeInfo<Lambda>::name()   == "{anonymous}::<lambda(int)>");
#elif defined(_MSC_VER)
static_assert(TypeInfo<Foo>::name()      == "struct `anonymous-namespace'::Foo");
static_assert(TypeInfo<FooAlias>::name() == "struct `anonymous-namespace'::Foo");
static_assert(TypeInfo<Bar>::name()      == "enum `anonymous-namespace'::Bar");
static_assert(TypeInfo<Baz>::name()      == "union `anonymous-namespace'::Baz");
static_assert(TypeInfo<Lambda>::name().starts_with("class `anonymous-namespace'::<lambda_"));
// underscore followed by some 32-bit hash that seems sensitive to the content of the current source code file
static_assert(TypeInfo<Lambda>::name().ends_with(">"));
#else
#error how did I ended up here?
#endif
// clang-format on

}  // namespace
