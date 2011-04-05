/* ----------------------------------------------------------------------
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

// smallint = variables for on-procesor system (nlocal, nmax, etc)
// tagint = variables for atom IDs (tag)
// bigint = variables for total system (natoms, ntimestep, etc)

// smallint must be an int, as defined by C compiler
// tagint can be 32-bit or 64-bit int, must be >= smallint
// NOTE: 64-bit tagint is not yet supported
// bigint can be 32-bit or 64-bit int, must be >= tagint

// MPI_LMP_TAGINT = MPI data type corresponding to a tagint
// MPI_LMP_BIGINT = MPI data type corresponding to a bigint

// NOTE: if your machine/MPI does not support "long long" ints,
//       but only "long" ints, then you will need to change
//       MPI_LONG_LONG to MPI_LONG, and atoll to atol

#ifndef LMP_LMPTYPE_H
#define LMP_LMPTYPE_H

#define __STDC_LIMIT_MACROS
#define __STDC_FORMAT_MACROS

#include "limits.h"
#include "stdint.h"
#include "inttypes.h"

namespace LAMMPS_NS {

// reserve 2 hi bits in molecular system neigh list for special bonds flag
// max local + ghost atoms per processor = 2^30 - 1

#define SBBITS 30
#define NEIGHMASK 0x3FFFFFFF

// default settings
// 32-bit smallint and tagint, 64-bit bigint

typedef int smallint;
typedef int tagint;
typedef int64_t bigint;

#define MAXSMALLINT INT_MAX
#define MAXTAGINT INT_MAX
#define MAXBIGINT INT64_MAX

#define MPI_LMP_TAGINT MPI_INT
#define MPI_LMP_BIGINT MPI_LONG_LONG

#define TAGINT_FORMAT "%d"
#define BIGINT_FORMAT "%" PRId64

#define ATOTAGINT atoi
#define ATOBIGINT atoll

// for molecular problems that exceed 2 billion (2^31) atoms
// 32-bit smallint, 64-bit tagint and bigint
// NOTE: 64-bit tagint is not yet supported

/*
typedef int smallint;
typedef int64_t tagint;
typedef int64_t bigint;

#define MAXSMALLINT INT_MAX
#define MAXTAGINT INT64_MAX
#define MAXBIGINT INT64_MAX

#define MPI_LMP_TAGINT MPI_LONG_LONG
#define MPI_LMP_BIGINT MPI_LONG_LONG

#define TAGINT_FORMAT "%" PRId64
#define BIGINT_FORMAT "%" PRId64

#define ATOTAGINT atoll
#define ATOBIGINT atoll
*/

// for machines that do not support 64-bit ints
// 32-bit smallint and tagint and bigint

/*
typedef int smallint;
typedef int tagint;
typedef int bigint;

#define MAXSMALLINT INT_MAX
#define MAXTAGINT INT_MAX
#define MAXBIGINT INT_MAX

#define MPI_LMP_TAGINT MPI_INT
#define MPI_LMP_BIGINT MPI_INT

#define TAGINT_FORMAT "%d"
#define BIGINT_FORMAT "%d"

#define ATOTAGINT atoi
#define ATOBIGINT atoi
*/

}

#endif
