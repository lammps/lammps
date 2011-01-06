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

// smallint = variables for system on 1 processor (nlocal, etc)
// tagint = variables for atom IDs (tag)
// bigint = variables for total system and timesteps (natoms, ntimestep, etc)

// smallint must be an int, as defined by C compiler
// tagint can be 32-bit or 64-bit int, must be >= smallint
// bigint can be 32-bit or 64-bit int, must be >= smallint and >= tagint

// MAXSMALLINT = max value of a smallint
// MAXTAGINT = max value of a tagint
// MAXBIGINT = max value of a bigint

// MPI_LMP_TAGINT = MPI data type corresponding to tagint
// MPI_LMP_BIGINT = MPI data type corresponding to bigint

#ifndef LMP_LMPTYPE_H
#define LMP_LMPTYPE_H

#include "stdint.h"

namespace LAMMPS_NS {

// default settings: 4-byte smallint, 4-byte tagint, 8-byte bigint

typedef int smallint;
typedef int tagint;
typedef int64_t bigint;

#define MAXSMALLINT 0x7FFFFFFF
#define MAXTAGINT 0x7FFFFFFF
#define MAXBIGINT 0x7FFFFFFFFFFFFFFFL

#define MPI_LMP_TAGINT MPI_INT
#define MPI_LMP_BIGINT MPI_LONG_LONG

}

#endif
