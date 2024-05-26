/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

// define integer data types used by LAMMPS and associated size limits

// smallint = variables for on-processor system (nlocal, nmax, etc)
// imageint = variables for atom image flags (image)
// tagint = variables for atom IDs and molecule IDs (tag,molecule)
// bigint = variables for total system (natoms, ntimestep, etc)

// smallint must be an int, as defined by C compiler
// imageint can be 32-bit or 64-bit int, must be >= smallint
// tagint can be 32-bit or 64-bit int, must be >= smallint
// bigint can be 32-bit or 64-bit int, must be >= imageint,tagint

// MPI_LMP_BIGINT = MPI data type corresponding to a bigint

#ifndef LMP_LMPTYPE_H
#define LMP_LMPTYPE_H

// C++11 check

#if __cplusplus < 201103L
#error LAMMPS requires a C++11 (or later) compliant compiler. Enable C++11 compatibility or upgrade the compiler.
#endif

#ifndef __STDC_LIMIT_MACROS
#define __STDC_LIMIT_MACROS
#endif

#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif

#include <cinttypes>    // IWYU pragma: export
#include <climits>      // IWYU pragma: export
#include <cstdlib>      // IWYU pragma: export

// grrr - IBM Power6 does not provide this def in their system header files

#ifndef PRId64
#define PRId64 "ld"
#endif

namespace LAMMPS_NS {

// reserve 2 highest bits in molecular system neigh list for special bonds flag
// reserve 3rd highest bit in neigh list for fix neigh/history flag
// max local + ghost atoms per processor = 2^29 - 1

#define SBBITS 30
#define HISTBITS 29
#define NEIGHMASK 0x1FFFFFFF
#define HISTMASK 0xDFFFFFFF
#define SPECIALMASK 0x3FFFFFFF

// default to 32-bit smallint and other ints, 64-bit bigint

#if !defined(LAMMPS_SMALLSMALL) && !defined(LAMMPS_BIGBIG) && !defined(LAMMPS_SMALLBIG)
#define LAMMPS_SMALLBIG
#endif

// allow user override of LONGLONG to LONG, necessary for some machines/MPI

#ifdef LAMMPS_LONGLONG_TO_LONG
#define MPI_LL MPI_LONG
#define ATOLL atoll
#else
#define MPI_LL MPI_LONG_LONG
#define ATOLL atol
#endif

// for atomic problems that exceed 2 billion (2^31) atoms
// 32-bit smallint/imageint/tagint, 64-bit bigint
// atom IDs and molecule IDs are limited to 32-bit

#ifdef LAMMPS_SMALLBIG

typedef int smallint;
typedef int imageint;
typedef int tagint;
typedef int64_t bigint;

#define MAXSMALLINT INT_MAX
#define MAXTAGINT INT_MAX
#define MAXBIGINT INT64_MAX
#define MAXDOUBLEINT 9007199254740992    // 2^53

#define MPI_LMP_TAGINT MPI_INT
#define MPI_LMP_IMAGEINT MPI_INT
#define MPI_LMP_BIGINT MPI_LL

#define TAGINT_FORMAT "%d"
#define BIGINT_FORMAT "%" PRId64

#define ATOTAGINT atoi
#define ATOBIGINT ATOLL

#define LAMMPS_TAGINT LAMMPS_INT
#define LAMMPS_TAGINT_2D LAMMPS_INT_2D
#define LAMMPS_BIGINT LAMMPS_INT64
#define LAMMPS_BIGINT_2D LAMMPS_INT64_2D

#define IMGMASK 1023
#define IMGMAX 512
#define IMGBITS 10
#define IMG2BITS 20

#endif

// for molecular problems that exceed 2 billion (2^31) atoms
// or problems where atoms wrap around the periodic box more than 512 times
// 32-bit smallint, 64-bit imageint/tagint/bigint

#ifdef LAMMPS_BIGBIG

typedef int smallint;
typedef int64_t imageint;
typedef int64_t tagint;
typedef int64_t bigint;

#define MAXSMALLINT INT_MAX
#define MAXTAGINT INT64_MAX
#define MAXBIGINT INT64_MAX
#define MAXDOUBLEINT 9007199254740992    // 2^53

#define MPI_LMP_TAGINT MPI_LL
#define MPI_LMP_IMAGEINT MPI_LL
#define MPI_LMP_BIGINT MPI_LL

#define TAGINT_FORMAT "%" PRId64
#define BIGINT_FORMAT "%" PRId64

#define ATOTAGINT ATOLL
#define ATOBIGINT ATOLL

#define LAMMPS_TAGINT LAMMPS_INT64
#define LAMMPS_TAGINT_2D LAMMPS_INT64_2D
#define LAMMPS_BIGINT LAMMPS_INT64
#define LAMMPS_BIGINT_2D LAMMPS_INT64_2D

#define IMGMASK 2097151
#define IMGMAX 1048576
#define IMGBITS 21
#define IMG2BITS 42

#endif

// for machines that do not support 64-bit ints
// 32-bit smallint/imageint/tagint/bigint

#ifdef LAMMPS_SMALLSMALL

typedef int smallint;
typedef int imageint;
typedef int tagint;
typedef int bigint;

#define MAXSMALLINT INT_MAX
#define MAXTAGINT INT_MAX
#define MAXBIGINT INT_MAX
#define MAXDOUBLEINT INT_MAX

#define MPI_LMP_TAGINT MPI_INT
#define MPI_LMP_IMAGEINT MPI_INT
#define MPI_LMP_BIGINT MPI_INT

#define TAGINT_FORMAT "%d"
#define BIGINT_FORMAT "%d"

#define ATOTAGINT atoi
#define ATOBIGINT atoi

#define LAMMPS_TAGINT LAMMPS_INT
#define LAMMPS_TAGINT_2D LAMMPS_INT_2D
#define LAMMPS_BIGINT LAMMPS_INT
#define LAMMPS_BIGINT_2D LAMMPS_INT_2D

#define IMGMASK 1023
#define IMGMAX 512
#define IMGBITS 10
#define IMG2BITS 20

#endif

/** Data structure for packing 32-bit and 64-bit integers
 * into double (communication) buffers
 *
 * Using this union avoids aliasing issues by having member types
 * (double, int) referencing the same buffer memory location.
 *
 * The explicit constructor for 32-bit integers prevents compilers
 * from (incorrectly) calling the double constructor when storing
 * an int into a double buffer.
\verbatim embed:rst

**Usage:**

.. code-block:: c++
   :caption: To copy an integer into a double buffer:

   double buf[2];
   int    foo =   1;
   tagint bar = 2<<40;
   buf[1] = ubuf(foo).d;
   buf[2] = ubuf(bar).d;

.. code-block:: c++
   :caption: To copy from a double buffer back to an int:

   foo = (int)    ubuf(buf[1]).i;
   bar = (tagint) ubuf(buf[2]).i;

The typecasts prevent compiler warnings about possible truncation issues.
\endverbatim
  */
union ubuf {
  double d;
  int64_t i;
  ubuf(const double &arg) : d(arg) {}
  ubuf(const int64_t &arg) : i(arg) {}
  ubuf(const int &arg) : i(arg) {}
};

/** Data structure for dynamic typing of int, bigint, and double
 *
 * Using this union allows to store any of the supported data types
 * in the same container and allows to "see" its current type.
\verbatim embed:rst

**Usage:**

.. code-block:: c++
   :caption: To store data in multitype array:

   multitype m[5];
   int    foo = 1;
   double bar = 2.5;
   bigint baz = 1<<40 - 1;
   m[0] = foo;
   m[1] = bar;
   m[2] = -1;
   m[3] = 2.0;
   m[4] = baz;

.. code-block:: c++
   :caption: To format data from multitype array into a space separated string:

   std::string str;
   for (int i = 0; i < 5; ++i) {
       switch (m[i].type) {
           case multitype::DOUBLE:
               str += std::to_string(m[i].data.d) + ' ';
               break;
           case multitype::INT:
               str += std::to_string(m[i].data.i) + ' ';
               break;
           case multitype::BIGINT:
               str += std::to_string(m[i].data.b) + ' ';
               break;
           default:
               break;
       }
   }
\endverbatim
  */
struct multitype {
  /** Data type constants for extracting data from atoms, computes and fixes
   *
   * This enum must be kept in sync with the corresponding enum or constants
   * in ``python/lammps/constants.py``, ``fortran/lammps.f90``, ``tools/swig/lammps.i``,
   * ``src/library.h``, and ``examples/COUPLE/plugin/liblammpsplugin.h`` */
  enum _LMP_DATATYPE_CONST {
    LAMMPS_NONE = -1,     /*!< no data type assigned (yet) */
    LAMMPS_INT = 0,       /*!< 32-bit integer (array) */
    LAMMPS_INT_2D = 1,    /*!< two-dimensional 32-bit integer array */
    LAMMPS_DOUBLE = 2,    /*!< 64-bit double (array) */
    LAMMPS_DOUBLE_2D = 3, /*!< two-dimensional 64-bit double array */
    LAMMPS_INT64 = 4,     /*!< 64-bit integer (array) */
    LAMMPS_INT64_2D = 5,  /*!< two-dimensional 64-bit integer array */
    LAMMPS_STRING = 6     /*!< C-String */
  };

  int type;
  union {
    double d;
    int i;
    int64_t b;
  } data;

  multitype() : type(LAMMPS_NONE) { data.d = 0.0; }
  multitype(const multitype &) = default;
  multitype(multitype &&) = default;
  ~multitype() = default;

  multitype &operator=(const double &_d)
  {
    type = LAMMPS_DOUBLE;
    data.d = _d;
    return *this;
  }
  multitype &operator=(const int &_i)
  {
    type = LAMMPS_INT;
    data.i = _i;
    return *this;
  }
  multitype &operator=(const int64_t &_b)
  {
    type = LAMMPS_INT64;
    data.b = _b;
    return *this;
  }
};

}    // namespace LAMMPS_NS

// preprocessor macros for compiler specific settings
// clear previous definitions to avoid redefinition warning

#ifdef _alignvar
#undef _alignvar
#endif
#ifdef _noalias
#undef _noalias
#endif
#ifdef _noopt
#undef _noopt
#endif

// define stack variable alignment

#if defined(__INTEL_COMPILER)
#define _alignvar(expr, val) __declspec(align(val)) expr
#elif defined(__GNUC__) || defined(__PGI) || defined(__INTEL_LLVM_COMPILER)
#define _alignvar(expr, val) expr __attribute((aligned(val)))
#else
#define _alignvar(expr, val) expr
#endif

// declaration to lift aliasing restrictions

#if defined(__INTEL_COMPILER) || (defined(__PGI) && !defined(__NVCOMPILER))
#define _noalias restrict
#elif defined(__GNUC__) || defined(__INTEL_LLVM_COMPILER) || defined(__NVCOMPILER)
#define _noalias __restrict
#else
#define _noalias
#endif

// Declaration to turn off optimization for specific noncritical
// functions and avoid compiler warnings about variable tracking.
// Disable for broken -D_FORTIFY_SOURCE feature.

#if defined(__clang__)
#define _noopt __attribute__((optnone))
#elif defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER)
#define _noopt
#elif defined(__PGI)
#define _noopt
#elif defined(__GNUC__)
#if (__GNUC__ > 4) || ((__GNUC__ == 4) && (__GNUC_MINOR__ >= 9))
#if defined(_FORTIFY_SOURCE) && (_FORTIFY_SOURCE > 0)
#define _noopt __attribute__((optimize("no-var-tracking-assignments")))
#else
#define _noopt __attribute__((optimize("O0", "no-var-tracking-assignments")))
#endif
#else
#if defined(_FORTIFY_SOURCE) && (_FORTIFY_SOURCE > 0)
#define _noopt
#else
#define _noopt __attribute__((optimize("O0")))
#endif
#endif
#else
#define _noopt
#endif

// suppress unused parameter warning

#define LMP_UNUSED_PARAM(x) (void) (x)

#endif
