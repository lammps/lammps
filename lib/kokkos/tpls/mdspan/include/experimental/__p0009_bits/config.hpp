//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER
#pragma once

#ifndef MDSPAN_IMPL_HAS_INCLUDE
#  ifndef __has_include
#    define MDSPAN_IMPL_HAS_INCLUDE(x) 0
#  else
#    define MDSPAN_IMPL_HAS_INCLUDE(x) __has_include(x)
#  endif
#endif

#if MDSPAN_IMPL_HAS_INCLUDE(<version>)
#  include <version>
#else
#  include <type_traits>
#  include <utility>
#endif

#ifdef _MSVC_LANG
#define MDSPAN_IMPL_CPLUSPLUS _MSVC_LANG
#else
#define MDSPAN_IMPL_CPLUSPLUS __cplusplus
#endif

#define MDSPAN_CXX_STD_14 201402L
#define MDSPAN_CXX_STD_17 201703L
#define MDSPAN_CXX_STD_20 202002L
// Note GCC has not updated this in version 13
#ifdef __clang__
#define MDSPAN_CXX_STD_23 202302L
#else
#define MDSPAN_CXX_STD_23 202100L
#endif

#define MDSPAN_HAS_CXX_14 (MDSPAN_IMPL_CPLUSPLUS >= MDSPAN_CXX_STD_14)
#define MDSPAN_HAS_CXX_17 (MDSPAN_IMPL_CPLUSPLUS >= MDSPAN_CXX_STD_17)
#define MDSPAN_HAS_CXX_20 (MDSPAN_IMPL_CPLUSPLUS >= MDSPAN_CXX_STD_20)
#define MDSPAN_HAS_CXX_23 (MDSPAN_IMPL_CPLUSPLUS >= MDSPAN_CXX_STD_23)

static_assert(MDSPAN_IMPL_CPLUSPLUS >= MDSPAN_CXX_STD_14, "mdspan requires C++14 or later.");

#ifndef MDSPAN_IMPL_COMPILER_CLANG
#  if defined(__clang__)
#    define MDSPAN_IMPL_COMPILER_CLANG __clang__
#  endif
#endif

#if !defined(MDSPAN_IMPL_COMPILER_MSVC) && !defined(MDSPAN_IMPL_COMPILER_MSVC_CLANG)
#  if defined(_MSC_VER)
#    if !defined(MDSPAN_IMPL_COMPILER_CLANG)
#      define MDSPAN_IMPL_COMPILER_MSVC _MSC_VER
#    else
#      define MDSPAN_IMPL_COMPILER_MSVC_CLANG _MSC_VER
#    endif
#  endif
#endif

#ifndef MDSPAN_IMPL_COMPILER_INTEL
#  ifdef __INTEL_COMPILER
#    define MDSPAN_IMPL_COMPILER_INTEL __INTEL_COMPILER
#  endif
#endif

#ifndef MDSPAN_IMPL_COMPILER_APPLECLANG
#  ifdef __apple_build_version__
#    define MDSPAN_IMPL_COMPILER_APPLECLANG __apple_build_version__
#  endif
#endif

#ifndef MDSPAN_IMPL_HAS_CUDA
#  if defined(__CUDACC__)
#    define MDSPAN_IMPL_HAS_CUDA __CUDACC__
#  endif
#endif

#ifndef MDSPAN_IMPL_HAS_HIP
#  if defined(__HIPCC__)
#    define MDSPAN_IMPL_HAS_HIP __HIPCC__
#  endif
#endif

#ifndef MDSPAN_IMPL_HAS_SYCL
#  if defined(SYCL_LANGUAGE_VERSION)
#    define MDSPAN_IMPL_HAS_SYCL SYCL_LANGUAGE_VERSION
#  endif
#endif

#ifndef MDSPAN_IMPL_HAS_CPP_ATTRIBUTE
#  ifndef __has_cpp_attribute
#    define MDSPAN_IMPL_HAS_CPP_ATTRIBUTE(x) 0
#  else
#    define MDSPAN_IMPL_HAS_CPP_ATTRIBUTE(x) __has_cpp_attribute(x)
#  endif
#endif

#ifndef MDSPAN_IMPL_PRESERVE_STANDARD_LAYOUT
// Preserve standard layout by default, but we're not removing the old version
// that turns this off until we're sure this doesn't have an unreasonable cost
// to the compiler or optimizer.
#  define MDSPAN_IMPL_PRESERVE_STANDARD_LAYOUT 1
#endif

#if !defined(MDSPAN_IMPL_USE_ATTRIBUTE_NO_UNIQUE_ADDRESS)
#  if ((MDSPAN_IMPL_HAS_CPP_ATTRIBUTE(no_unique_address) >= 201803L) && \
       (!defined(__NVCC__) || ((__CUDACC_VER_MAJOR__ * 100 + __CUDACC_VER_MINOR__ * 10 < 1290) && MDSPAN_HAS_CXX_20)) && \
       (!defined(MDSPAN_IMPL_COMPILER_MSVC) || MDSPAN_HAS_CXX_20))
#    define MDSPAN_IMPL_USE_ATTRIBUTE_NO_UNIQUE_ADDRESS 1
#    define MDSPAN_IMPL_NO_UNIQUE_ADDRESS [[no_unique_address]]
#  else
#    define MDSPAN_IMPL_NO_UNIQUE_ADDRESS
#  endif
#endif

// NVCC older than 11.6 chokes on the no-unique-address-emulation
// so just pretend to use it (to avoid the full blown EBO workaround
// which NVCC also doesn't like ...), and leave the macro empty
#ifndef MDSPAN_IMPL_NO_UNIQUE_ADDRESS
#  if defined(__NVCC__)
#    define MDSPAN_IMPL_USE_ATTRIBUTE_NO_UNIQUE_ADDRESS 1
#    define MDSPAN_IMPL_USE_FAKE_ATTRIBUTE_NO_UNIQUE_ADDRESS
#  endif
#  define MDSPAN_IMPL_NO_UNIQUE_ADDRESS
#endif

// AMDs HIP compiler seems to have issues with concepts
// it pretends concepts exist, but doesn't ship <concept>
#ifndef __HIPCC__
#ifndef MDSPAN_IMPL_USE_CONCEPTS
#  if defined(__cpp_concepts) && __cpp_concepts >= 201507L
#    define MDSPAN_IMPL_USE_CONCEPTS 1
#  endif
#endif
#endif

#ifndef MDSPAN_IMPL_USE_FOLD_EXPRESSIONS
#  if (defined(__cpp_fold_expressions) && __cpp_fold_expressions >= 201603L) \
          || (!defined(__cpp_fold_expressions) && MDSPAN_HAS_CXX_17)
#    define MDSPAN_IMPL_USE_FOLD_EXPRESSIONS 1
#  endif
#endif

#ifndef MDSPAN_IMPL_USE_INLINE_VARIABLES
#  if defined(__cpp_inline_variables) && __cpp_inline_variables >= 201606L \
         || (!defined(__cpp_inline_variables) && MDSPAN_HAS_CXX_17)
#    define MDSPAN_IMPL_USE_INLINE_VARIABLES 1
#  endif
#endif

#ifndef MDSPAN_IMPL_NEEDS_TRAIT_VARIABLE_TEMPLATE_BACKPORTS
#  if (!(defined(__cpp_lib_type_trait_variable_templates) && __cpp_lib_type_trait_variable_templates >= 201510L) \
          || !MDSPAN_HAS_CXX_17)
#    if !(defined(MDSPAN_IMPL_COMPILER_APPLECLANG) && MDSPAN_HAS_CXX_17)
#      define MDSPAN_IMPL_NEEDS_TRAIT_VARIABLE_TEMPLATE_BACKPORTS 1
#    endif
#  endif
#endif

#ifndef MDSPAN_IMPL_USE_VARIABLE_TEMPLATES
#  if (defined(__cpp_variable_templates) && __cpp_variable_templates >= 201304 && MDSPAN_HAS_CXX_17) \
        || (!defined(__cpp_variable_templates) && MDSPAN_HAS_CXX_17)
#    define MDSPAN_IMPL_USE_VARIABLE_TEMPLATES 1
#  endif
#endif // MDSPAN_IMPL_USE_VARIABLE_TEMPLATES

#ifndef MDSPAN_IMPL_USE_CONSTEXPR_14
#  if (defined(__cpp_constexpr) && __cpp_constexpr >= 201304) \
        || (!defined(__cpp_constexpr) && MDSPAN_HAS_CXX_14) \
        && (!(defined(__INTEL_COMPILER) && __INTEL_COMPILER <= 1700))
#    define MDSPAN_IMPL_USE_CONSTEXPR_14 1
#  endif
#endif

#ifndef MDSPAN_IMPL_USE_IF_CONSTEXPR_17
#  if (defined(__cpp_if_constexpr) && __cpp_if_constexpr >= 201606) \
        || (!defined(__cpp_constexpr) && MDSPAN_HAS_CXX_17)
#    define MDSPAN_IMPL_USE_IF_CONSTEXPR_17 1
#  endif
#endif

#ifndef MDSPAN_IMPL_USE_INTEGER_SEQUENCE_14
#  if defined(MDSPAN_IMPL_COMPILER_MSVC)
#    if (defined(__cpp_lib_integer_sequence) && __cpp_lib_integer_sequence >= 201304)
#      define MDSPAN_IMPL_USE_INTEGER_SEQUENCE_14 1
#    endif
#  endif
#endif
#ifndef MDSPAN_IMPL_USE_INTEGER_SEQUENCE_14
#  if (defined(__cpp_lib_integer_sequence) && __cpp_lib_integer_sequence >= 201304) \
        || (!defined(__cpp_lib_integer_sequence) && MDSPAN_HAS_CXX_14) \
        /* as far as I can tell, libc++ seems to think this is a C++11 feature... */ \
        || (defined(__GLIBCXX__) && __GLIBCXX__ > 20150422 && __GNUC__ < 5 && !defined(__INTEL_CXX11_MODE__))
     // several compilers lie about integer_sequence working properly unless the C++14 standard is used
#    define MDSPAN_IMPL_USE_INTEGER_SEQUENCE_14 1
#  elif defined(MDSPAN_IMPL_COMPILER_APPLECLANG) && MDSPAN_HAS_CXX_14
     // appleclang seems to be missing the __cpp_lib_... macros, but doesn't seem to lie about C++14 making
     // integer_sequence work
#    define MDSPAN_IMPL_USE_INTEGER_SEQUENCE_14 1
#  endif
#endif

#ifndef MDSPAN_IMPL_USE_RETURN_TYPE_DEDUCTION
#  if (defined(__cpp_return_type_deduction) && __cpp_return_type_deduction >= 201304) \
          || (!defined(__cpp_return_type_deduction) && MDSPAN_HAS_CXX_14)
#    define MDSPAN_IMPL_USE_RETURN_TYPE_DEDUCTION 1
#  endif
#endif

#ifndef MDSPAN_IMPL_USE_CLASS_TEMPLATE_ARGUMENT_DEDUCTION
#  if (!defined(__NVCC__) || (__CUDACC_VER_MAJOR__ * 100 + __CUDACC_VER_MINOR__ * 10 >= 1170)) && \
      ((defined(__cpp_deduction_guides) && __cpp_deduction_guides >= 201703) || \
       (!defined(__cpp_deduction_guides) && MDSPAN_HAS_CXX_17))
#    define MDSPAN_IMPL_USE_CLASS_TEMPLATE_ARGUMENT_DEDUCTION 1
#  endif
#endif

#ifndef MDSPAN_IMPL_USE_STANDARD_TRAIT_ALIASES
#  if (defined(__cpp_lib_transformation_trait_aliases) && __cpp_lib_transformation_trait_aliases >= 201304) \
          || (!defined(__cpp_lib_transformation_trait_aliases) && MDSPAN_HAS_CXX_14)
#    define MDSPAN_IMPL_USE_STANDARD_TRAIT_ALIASES 1
#  elif defined(MDSPAN_IMPL_COMPILER_APPLECLANG) && MDSPAN_HAS_CXX_14
     // appleclang seems to be missing the __cpp_lib_... macros, but doesn't seem to lie about C++14
#    define MDSPAN_IMPL_USE_STANDARD_TRAIT_ALIASES 1
#  endif
#endif

#ifndef MDSPAN_IMPL_DEFAULTED_CONSTRUCTORS_INHERITANCE_WORKAROUND
#  ifdef __GNUC__
#    if __GNUC__ < 9
#      define MDSPAN_IMPL_DEFAULTED_CONSTRUCTORS_INHERITANCE_WORKAROUND 1
#    endif
#  endif
#endif

#ifndef MDSPAN_CONDITIONAL_EXPLICIT
#  if MDSPAN_HAS_CXX_20
#    define MDSPAN_CONDITIONAL_EXPLICIT(COND) explicit(COND)
#  else
#    define MDSPAN_CONDITIONAL_EXPLICIT(COND)
#  endif
#endif

#ifndef MDSPAN_USE_BRACKET_OPERATOR
#  if defined(__cpp_multidimensional_subscript)
// The following if/else is necessary to workaround a clang issue
// relative to using a parameter pack inside a bracket operator in C++2b/C++23 mode
#    if defined(MDSPAN_IMPL_COMPILER_CLANG) &&                                         \
        ((__clang_major__ < 17) ||                                                 \
         (__clang_major__ == 17 && __clang_minor__ == 0 &&                         \
          __clang_patchlevel__ == 0))
#      define MDSPAN_USE_BRACKET_OPERATOR 0
#    else
#      define MDSPAN_USE_BRACKET_OPERATOR 1
#    endif
#  else
#    define MDSPAN_USE_BRACKET_OPERATOR 0
#  endif
#endif

#ifndef MDSPAN_USE_PAREN_OPERATOR
#  if !MDSPAN_USE_BRACKET_OPERATOR
#    define MDSPAN_USE_PAREN_OPERATOR 1
#  else
#    define MDSPAN_USE_PAREN_OPERATOR 0
#  endif
#endif

#if MDSPAN_USE_BRACKET_OPERATOR
#  define MDSPAN_IMPL_OP(mds,...) mds[__VA_ARGS__]
// Corentins demo compiler for subscript chokes on empty [] call,
// though I believe the proposal supports it?
#ifdef MDSPAN_NO_EMPTY_BRACKET_OPERATOR
#  define MDSPAN_IMPL_OP0(mds) mds.accessor().access(mds.data_handle(),0)
#else
#  define MDSPAN_IMPL_OP0(mds) mds[]
#endif
#  define MDSPAN_IMPL_OP1(mds, a) mds[a]
#  define MDSPAN_IMPL_OP2(mds, a, b) mds[a,b]
#  define MDSPAN_IMPL_OP3(mds, a, b, c) mds[a,b,c]
#  define MDSPAN_IMPL_OP4(mds, a, b, c, d) mds[a,b,c,d]
#  define MDSPAN_IMPL_OP5(mds, a, b, c, d, e) mds[a,b,c,d,e]
#  define MDSPAN_IMPL_OP6(mds, a, b, c, d, e, f) mds[a,b,c,d,e,f]
#else
#  define MDSPAN_IMPL_OP(mds,...) mds(__VA_ARGS__)
#  define MDSPAN_IMPL_OP0(mds) mds()
#  define MDSPAN_IMPL_OP1(mds, a) mds(a)
#  define MDSPAN_IMPL_OP2(mds, a, b) mds(a,b)
#  define MDSPAN_IMPL_OP3(mds, a, b, c) mds(a,b,c)
#  define MDSPAN_IMPL_OP4(mds, a, b, c, d) mds(a,b,c,d)
#  define MDSPAN_IMPL_OP5(mds, a, b, c, d, e) mds(a,b,c,d,e)
#  define MDSPAN_IMPL_OP6(mds, a, b, c, d, e, f) mds(a,b,c,d,e,f)
#endif
