/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

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

#ifndef __STDC_LIMIT_MACROS
#define __STDC_LIMIT_MACROS
#endif

#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif

#include <climits>
#include <stdint.h>   // <cstdint> requires C++-11
#include <inttypes.h> // <cinttypes> requires C++-11

// grrr - IBM Power6 does not provide this def in their system header files

#ifndef PRId64
#define PRId64 "ld"
#endif

namespace LAMMPS_NS {

// enum used for KOKKOS host/device flags

enum ExecutionSpace{Host,Device};

// reserve 2 hi bits in molecular system neigh list for special bonds flag
// max local + ghost atoms per processor = 2^30 - 1

#define SBBITS 30
#define NEIGHMASK 0x3FFFFFFF

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

#ifdef LAMMPS_SMALLBIG

typedef int smallint;
typedef int imageint;
typedef int tagint;
typedef int64_t bigint;

#define MAXSMALLINT INT_MAX
#define MAXTAGINT INT_MAX
#define MAXBIGINT INT64_MAX

#define MPI_LMP_TAGINT MPI_INT
#define MPI_LMP_IMAGEINT MPI_INT
#define MPI_LMP_BIGINT MPI_LL

#define TAGINT_FORMAT "%d"
#define BIGINT_FORMAT "%" PRId64

#define ATOTAGINT atoi
#define ATOBIGINT ATOLL

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

#define MPI_LMP_TAGINT MPI_LL
#define MPI_LMP_IMAGEINT MPI_LL
#define MPI_LMP_BIGINT MPI_LL

#define TAGINT_FORMAT "%" PRId64
#define BIGINT_FORMAT "%" PRId64

#define ATOTAGINT ATOLL
#define ATOBIGINT ATOLL

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

#define MPI_LMP_TAGINT MPI_INT
#define MPI_LMP_IMAGEINT MPI_INT
#define MPI_LMP_BIGINT MPI_INT

#define TAGINT_FORMAT "%d"
#define BIGINT_FORMAT "%d"

#define ATOTAGINT atoi
#define ATOBIGINT atoi

#define IMGMASK 1023
#define IMGMAX 512
#define IMGBITS 10
#define IMG2BITS 20

#endif

}

// preprocessor macros for compiler specific settings
// clear previous definitions to avoid redefinition warning

#ifdef _alignvar
#undef _alignvar
#endif
#ifdef _noalias
#undef _noalias
#endif

// define stack variable alignment

#if defined(__INTEL_COMPILER)
#define _alignvar(expr,val) __declspec(align(val)) expr
#elif defined(__GNUC__)
#define _alignvar(expr,val) expr __attribute((aligned(val)))
#else
#define _alignvar(expr,val) expr
#endif

// declaration to lift aliasing restrictions

#if defined(__INTEL_COMPILER)
#define _noalias restrict
#elif defined(__GNUC__)
#define _noalias __restrict
#else
#define _noalias
#endif

// settings to enable LAMMPS to build under Windows

#ifdef _WIN32
#include "lmpwindows.h"
#endif

// suppress unused parameter warning

#define LMP_UNUSED_PARAM(x) (void)(x)

#endif
